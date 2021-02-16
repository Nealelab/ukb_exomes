import argparse
import logging

import hail as hl

from gnomad.utils.slack import slack_notifications

from ukbb_qc.resources.basics import (
    get_doubleton_ht_path,
    get_ukbb_data,
    get_pair_ht_path,
    logging_path,
    release_ht_path,
)
from ukbb_qc.resources.resource_utils import CURRENT_FREEZE
from ukbb_qc.slack_creds import slack_token

logging.basicConfig(
    format="%(asctime)s (%(name)s %(lineno)s): %(message)s",
    datefmt="%m/%d/%Y %I:%M:%S %p",
)
logger = logging.getLogger("doubleton_analysis")
logger.setLevel(logging.INFO)


def get_doubletons(
    mt: hl.MatrixTable, unrelated_only: bool, freq_index: int = 0
) -> hl.Table:
    """
    Filters input MatrixTable to doubletons and annotates each doubletons with relevant sample IDs.

    Assumes input MatrixTable has the following annotations:
        - freq (row annotation containing array of frequency structs)
        - meta (column annotation containing sample metadata struct)
    Also assumes `freq` annotation was calculated using `annotate_freq` and `meta` annotation contains the fields
    `sample_filters` and `related`.

    :param hl.MatrixTable mt: Input MatrixTable.
    :param bool unrelated_only: Whether to get doubletons from unrelated samples only.
    :param int freq_index: Which index of freq struct array to use. Default is 0.
    :return: Table with doubletons and relevant sample IDs for each doubleton.
    :rtype: hl.Table
    """
    if unrelated_only:
        logger.info("Removing related samples and their variants...")
        mt = mt.filter_cols(~mt.meta.sample_filters.related)
        mt = mt.filter_rows(hl.agg.any(mt.GT.is_non_ref()))

    logger.info("Filtering to rows where AC == 2 and homozygote count == 0...")
    mt = mt.filter_rows(
        (mt.freq[freq_index].AC == 2) & (mt.freq[freq_index].homozygote_count == 0)
    )

    logger.info("Annotating each doubleton with sample IDs...")
    mt = mt.annotate_rows(pair=hl.agg.filter(mt.GT.is_het(), hl.agg.collect(mt.s)),)
    return mt.annotate_rows(s_1=mt.pair[0], s_2=mt.pair[1],).rows()


def get_random_pairs(mt: hl.MatrixTable, n_pairs: int) -> hl.Table:
    """
    Selects specified number of random sample pairs from input MatrixTable.

    :param hl.MatrixTable mt: Input MatrixTable.
    :param int n_pairs: Desired number of sample pairs.
    :return: Table of random sample pairs.
    :rtype: hl.Table
    """
    logger.info("Selecting 100k random samples...")
    mt = mt.add_col_index()
    indices = hl.range(mt.count_cols())
    rand_indices = hl.shuffle(indices)[: n_pairs - 1]
    mt = mt.annotate_cols(s_1=rand_indices.contains(mt.s))

    logger.info("Selecting a second set of 100k random samples...")
    mt = mt.annotate_cols(col_idx=hl.or_missing(~mt.s_1, mt.col_idx))
    indices = mt.aggregate(
        hl.agg.filter(hl.is_defined(mt.col_idx), hl.agg.collect(mt.col_idx()))
    )
    rand_indices = hl.shuffle(indices)[: n_pairs - 1]
    mt = mt.annotate_cols(s_2=rand_indices.contains(mt.s))
    return mt.cols()


def main(args):

    data_source = "broad"
    freeze = args.freeze
    tranche_data = (data_source, freeze)
    pops = set(args.pops_to_include.split(","))

    try:
        logger.info("Reading in adj-filtered hardcalls MT...")
        mt = get_ukbb_data(*tranche_data, adj=True, meta_root="meta")

        logger.info(f"Only including these populations: {pops}")
        mt = mt.annotate_cols(
            meta=mt.meta.annotate(
                pop=mt.meta.hybrid_pop_data.pop
                if args.use_hybrid_pop
                else mt.meta.pan_ancestry_meta.pop
            )
        )
        logger.info(
            f"Removing variants from samples outside of included populations..."
        )
        mt = mt.filter_cols(hl.literal(pops).contains(mt.pop))
        mt = mt.filter_rows(hl.agg.any(mt.GT.is_non_ref()))

        if args.get_doubletons:
            logger.info("Reading in release HT and selecting frequency field...")
            release_ht = (
                hl.read_table(release_ht_path(*tranche_data))
                .select_globals()
                .select("freq")
            )
            mt = mt.annotate(freq=release_ht[mt.row_key].freq)
            ht = get_doubletons(mt, args.unrelated_only)
            ht.write(get_doubleton_ht_path(*tranche_data, args.unrelated_only))

        if args.get_random_pairs:
            ht = get_random_pairs(mt, args.n_pairs)
            ht.write(get_pair_ht_path(*tranche_data))

    finally:
        logger.info("Copying hail log to logging bucket...")
        hl.copy_log(logging_path(*tranche_data))


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        "This script extracts pairs of individuals sharing rare doubletons and compares them to ~100k randomly selected sample pairs"
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
        "--use_hybrid_pop",
        help="Use hybrid ancestry assignments when generating LoF matrix summary Table. Will use pan-ancestry assignment if not set",
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
