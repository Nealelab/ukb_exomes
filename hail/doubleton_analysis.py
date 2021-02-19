import argparse
import logging

import hail as hl

from gnomad.utils.slack import slack_notifications

from ukbb_qc.resources.basics import (
    array_sample_map_ht_path,
    geographical_ht_path,
    get_checkpoint_path,
    get_doubleton_ht_path,
    get_pair_ht_path,
    get_ukbb_data,
    logging_path,
)
from ukbb_qc.resources.resource_utils import CURRENT_FREEZE
from ukbb_qc.slack_creds import slack_token

logging.basicConfig(
    format="%(asctime)s (%(name)s %(lineno)s): %(message)s",
    datefmt="%m/%d/%Y %I:%M:%S %p",
)
logger = logging.getLogger("doubleton_analysis")
logger.setLevel(logging.INFO)


def get_samples_with_geo_data(data_source: str, freeze: int, overwrite: bool) -> None:
    """
    Reads in geographical data Table, maps country ints to names, and removes samples with missing data.

    Overwrites HT at `geographical_ht_path()`.

    Also maps UKBB IDs to exome IDs and converts geographical coordinates from strings into ints.

    Assumes geographical data HT has the following fields:
        - ukbb_app_26041_id (key)
        - north 
        - east
        - country
    Also assumes all fields in geographical data HT are strings.

    :param str data_source: One of 'regeneron' or 'broad'. Used to checkpoint geographical HT.
    :param int freeze: One of the data freezes. Used to checkpoint geographical HT.
    :param bool: Whether to overwrite existing geographical HT.
    :return: None
    """
    logger.info("Reading in location HT and array sample map HT...")
    map_ht = hl.read_table(array_sample_map_ht_path())
    geo_ht = hl.read_table(geographical_ht_path())
    map_ht = map_ht.key_by("ukbb_app_26041_id")
    geo_ht = geo_ht.annotate(s=map_ht[geo_ht.key].s).key_by("s")
    geo_ht = geo_ht.drop("ukbb_app_26041_id")
    logger.info(f"Found {geo_ht.count()} rows in geographical HT...")

    logger.info("Mapping country ints to names...")
    # Mapping here: https://biobank.ndph.ox.ac.uk/showcase/coding.cgi?id=100420
    # Skipping country codes 6, -1, and -3:
    # 6 = "Elsewhere", -1 = "Do not know", -3 = "Prefer not to answer"
    geo_ht = geo_ht.transmute(
        country=hl.case()
        .when(geo_ht.country == "1", "England")
        .when(geo_ht.country == "2", "Wales")
        .when(geo_ht.country == "3", "Scotland")
        .when(geo_ht.country == "4", "Northern Ireland")
        .when(geo_ht.country == "5", "Republic of Ireland")
        .or_missing()
    )

    logger.info("Removing samples with missing geographical data...")
    geo_ht = geo_ht.filter(
        (
            hl.is_defined(geo_ht.north)
            & hl.is_defined(geo_ht.east)
            & hl.is_defined(geo_ht.country)
        )
        # NOTE: Samples born in Northern Ireland have no birth coordinates
        # Most samples born in the Republic of Ireland also have no birth coordinates
        | (geo_ht.country == "Northern Ireland")
        | (geo_ht.country == "Republic of Ireland")
    )

    logger.info("Converting north and east coordinates from strings to ints...")
    geo_ht = geo_ht.transmute(north=hl.int(geo_ht.north), east=hl.int(geo_ht.east),)

    geo_ht = geo_ht.checkpoint(
        get_checkpoint_path(data_source, freeze, name="geographical_data.ht"),
        overwrite=True,
    )
    geo_ht.write(geographical_ht_path(), overwrite=overwrite)
    logger.info(
        f"Kept {geo_ht.count()} rows after removing samples with missing data..."
    )
    logger.info(f"Country counter: {geo_ht.aggregate(hl.agg.counter(geo_ht.country))}")


def get_doubletons(mt: hl.MatrixTable) -> hl.Table:
    """
    Filters input MatrixTable to doubletons and annotates each doubletons with relevant sample IDs.

    Assumes input MatrixTable has the following annotations:
        - freq (row annotation containing struct of frequency information)
        - meta (column annotation containing sample metadata struct)
    Also assumes `freq` annotation was calculated using `hl.agg.call_stats` and `meta` annotation contains the fields
    `sample_filters` and `related`.

    :param hl.MatrixTable mt: Input MatrixTable.
    :return: Table with doubletons and relevant sample IDs for each doubleton.
    :rtype: hl.Table
    """
    logger.info("Filtering to rows where AC == 2 and homozygote count == 0...")
    mt = mt.filter_rows(
        # Need to check AC[1] and homozygote_count[1] here because
        # call_stats returns one element for each allele (including reference)
        (mt.freq.AC[1] == 2)
        & (mt.freq.homozygote_count[1] == 0)
    )

    logger.info("Annotating each doubleton with sample IDs...")
    mt = mt.annotate_rows(pair=hl.agg.filter(mt.GT.is_het(), hl.agg.collect(mt.s)))
    return mt.annotate_rows(s_1=mt.pair[0], s_2=mt.pair[1]).rows()


def get_random_pairs(
    ht: hl.Table, n_pairs: int, data_source: str, freeze: int
) -> hl.Table:
    """
    Selects specified number of random sample pairs from input Table.

    :param hl.Table ht: Input MatrixTable.
    :param int n_pairs: Desired number of sample pairs.
    :param str data_source: One of 'regeneron' or 'broad'. Used to checkpoint temporary HT.
    :param int freeze: One of the data freezes. Used to checkpoint temporary HT.
    :return: Table of random sample pairs.
    :rtype: hl.Table
    """
    logger.info(f"Selecting {n_pairs} random samples...")
    ht = ht.add_index()
    indices = hl.range(ht.count())
    rand_indices = hl.shuffle(indices)[: n_pairs - 1]
    ht = ht.annotate(s_1=rand_indices.contains(hl.int(ht.idx)))
    ht = ht.checkpoint(
        get_checkpoint_path(data_source, freeze, name="random_pairs_temp.ht"),
        overwrite=True,
    )

    logger.info(f"Selecting a second set of {n_pairs} random samples...")
    new_indices = hl.filter(lambda x: ~rand_indices.contains(x), indices)
    rand_indices = hl.shuffle(new_indices)[: n_pairs - 1]
    ht = ht.annotate(s_2=rand_indices.contains(ht.s))
    return ht


def main(args):

    hl.init(log="/doubleton_analysis.log", default_reference="GRCh38")
    data_source = "broad"
    freeze = args.freeze
    tranche_data = (data_source, freeze)
    pops = set(args.pops_to_include.split(","))

    try:
        if args.filter_geographical_data:
            get_samples_with_geo_data(*tranche_data, args.overwrite)

        logger.info("Reading in adj-filtered hardcalls MT...")
        mt = get_ukbb_data(*tranche_data, adj=True, meta_root="meta")

        logger.info("Filtering to autosomes only...")
        mt = mt.filter_rows(mt.locus.in_autosome())

        logger.info(f"Only including these populations: {pops}")
        mt = mt.annotate_cols(
            meta=mt.meta.annotate(
                pop=mt.meta.gnomad_pc_project_pop_data.pop
                if args.use_gnomad_pop
                else mt.meta.pan_ancestry_meta.pop
            )
        )
        mt = mt.filter_cols(hl.literal(pops).contains(mt.meta.pop))

        logger.info("Filtering MT to samples with geographical information...")
        geo_ht = hl.read_table(geographical_ht_path())
        mt = mt.annotate_cols(**geo_ht[mt.col_key])
        mt = mt.filter_cols(hl.is_defined(mt.country))
        mt = mt.filter_rows(hl.agg.any(mt.GT.is_non_ref()))

        if args.get_doubletons:
            logger.info("Calculating frequency using call_stats...")
            if args.unrelated_only:
                logger.info("Removing related samples and their variants...")
                mt = mt.filter_cols(~mt.meta.sample_filters.related)
                mt = mt.filter_rows(hl.agg.any(mt.GT.is_non_ref()))
            mt = mt.annotate_rows(freq=hl.agg.call_stats(mt.GT, mt.alleles))
            ht = get_doubletons(mt)
            ht.write(
                get_doubleton_ht_path(*tranche_data, args.unrelated_only),
                overwrite=args.overwrite,
            )

        if args.get_random_pairs:
            ht = get_random_pairs(mt.cols(), args.n_pairs, *tranche_data)
            ht.write(get_pair_ht_path(*tranche_data), overwrite=args.overwrite)

    finally:
        logger.info("Copying hail log to logging bucket...")
        hl.copy_log(logging_path(*tranche_data))


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        "This script extracts pairs of individuals sharing rare doubletons and compares them to ~100k randomly selected sample pairs"
    )
    parser.add_argument(
        "--filter_geographical_data",
        help="Filter geographical data HT to remove missing data",
        action="store_true",
    )
    parser.add_argument(
        "--get_doubletons", help="Get samples with doubletons", action="store_true"
    )
    parser.add_argument(
        "--unrelated_only", help="Keep unrelated samples only", action="store_true"
    )
    parser.add_argument(
        "--get_random_pairs",
        help="Extract random sample pairs from MatrixTable",
        action="store_true",
    )
    parser.add_argument(
        "--n_pairs", help="Desired number of sample pairs", type=int,
    )
    parser.add_argument(
        "--pops_to_include", help="Populations to include in this analysis",
    )
    parser.add_argument(
        "--use_gnomad_pop",
        help="Use gnomAD ancestry assignments when generating LoF matrix summary Table. Will use pan-ancestry assignment if not set",
        action="store_true",
    )
    parser.add_argument(
        "-f", "--freeze", help="Data freeze to use", default=CURRENT_FREEZE, type=int,
    )
    parser.add_argument(
        "--overwrite", help="Overwrite existing data", action="store_true"
    )
    parser.add_argument(
        "--slack_channel", help="Send message to Slack channel/user", default="@kc"
    )
    args = parser.parse_args()

    if args.slack_channel:
        with slack_notifications(slack_token, args.slack_channel):
            main(args)
    else:
        main(args)
