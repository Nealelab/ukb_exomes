#!/usr/bin/env python3

__author__ = 'konradk'

from ukb_exomes import *

MIN_CALL_RATE = 0.95
VARIANTS_PER_MAC_CATEGORY = 2000
VARIANTS_PER_MAF_CATEGORY = 10000


def create_gene_map_ht(ht, check_gene_contigs=False):
    ht = process_consequences(ht)
    ht = ht.explode(ht.vep.worst_csq_by_gene_canonical)
    ht = ht.annotate(
        variant_id=ht.locus.contig + ':' + hl.str(ht.locus.position) + '_' + ht.alleles[0] + '/' + ht.alleles[1],
        annotation=annotation_case_builder(ht.vep.worst_csq_by_gene_canonical))
    if check_gene_contigs:
        gene_contigs = ht.group_by(
            gene_id=ht.vep.worst_csq_by_gene_canonical.gene_id,
            gene_symbol=ht.vep.worst_csq_by_gene_canonical.gene_symbol,
        ).aggregate(
            contigs=hl.agg.collect_as_set(ht.locus.contig)
        )
        assert gene_contigs.all(hl.len(gene_contigs.contigs) == 1)

    gene_map_ht = ht.group_by(
        gene_id=ht.vep.worst_csq_by_gene_canonical.gene_id,
        gene_symbol=ht.vep.worst_csq_by_gene_canonical.gene_symbol,
    ).partition_hint(100).aggregate(
        interval=hl.interval(
            start=hl.locus(hl.agg.take(ht.locus.contig, 1)[0], hl.agg.min(ht.locus.position)),
            end=hl.locus(hl.agg.take(ht.locus.contig, 1)[0], hl.agg.max(ht.locus.position))
        ),
        variants=hl.agg.group_by(ht.annotation, hl.agg.collect(ht.variant_id)),
    )
    return gene_map_ht


def post_process_gene_map_ht(gene_ht):
    groups = ['pLoF', 'missense|LC', 'pLoF|missense|LC', 'synonymous', 'missense']
    variant_groups = hl.map(lambda group: group.split('\\|').flatmap(lambda csq: gene_ht.variants.get(csq)), groups)
    gene_ht = gene_ht.transmute(
        variant_groups=hl.zip(groups, variant_groups)
    ).explode('variant_groups')
    gene_ht = gene_ht.transmute(annotation=gene_ht.variant_groups[0], variants=hl.sorted(gene_ht.variant_groups[1]))
    gene_ht = gene_ht.key_by(start=gene_ht.interval.start)
    return gene_ht.filter(hl.len(gene_ht.variants) > 0)


def mac_category_case_builder(call_stats_expr):
    return (hl.case()
            .when(call_stats_expr.AC[1] <= 5, call_stats_expr.AC[1])
            .when(call_stats_expr.AC[1] <= 10, 10)
            .when(call_stats_expr.AC[1] <= 20, 20)
            .when(call_stats_expr.AF[1] <= 0.001, 0.001)
            .when(call_stats_expr.AF[1] <= 0.01, 0.01)
            .when(call_stats_expr.AF[1] <= 0.1, 0.1)
            .default(0.99))


def prepare_mt_for_plink(mt: hl.MatrixTable, min_call_rate: float = MIN_CALL_RATE,
                         variants_per_mac_category: int = VARIANTS_PER_MAC_CATEGORY,
                         variants_per_maf_category: int = VARIANTS_PER_MAF_CATEGORY):
    mt = filter_to_autosomes(mt)
    hq_samples = mt.aggregate_cols(hl.agg.count_where(mt.meta.high_quality))
    call_stats_ht = hl.read_table(ukb.var_annotations_ht_path(*TRANCHE_DATA[CURRENT_TRANCHE], 'call_stats'))
    mt = mt.annotate_rows(call_stats=call_stats_ht[mt.row_key].qc_callstats[0])
    mt = mt.filter_rows((mt.call_stats.AN >= hq_samples * 2 * min_call_rate) &
                        (mt.call_stats.AC[1] > 0))
    mt = mt.annotate_rows(mac_category=mac_category_case_builder(mt.call_stats))
    category_counter = mt.aggregate_rows(hl.agg.counter(mt.mac_category))
    print(category_counter)
    mt = mt.annotate_globals(category_counter=category_counter)
    mt = mt.filter_rows(hl.rand_unif(0, 1) <
                        hl.cond(mt.mac_category >= 1, variants_per_mac_category, variants_per_maf_category) / mt.category_counter[mt.mac_category]
                        )
    return mt


def main(args):
    hl.init(default_reference='GRCh38')

    if args.create_plink_file:
        mt = prepare_mt_for_plink(get_filtered_mt(adj=True, interval_filter=True))
        mt = mt.checkpoint(get_ukb_grm_mt_path(), _read_if_exists=not args.overwrite, overwrite=args.overwrite)
        mt = mt.naive_coalesce(1000).checkpoint(get_ukb_grm_mt_path(coalesced=True), _read_if_exists=not args.overwrite, overwrite=args.overwrite)
        mt = mt.unfilter_entries()
        ht = hl.ld_prune(mt.GT, r2=0.1)
        ht = ht.checkpoint(get_ukb_grm_pruned_ht_path(), _read_if_exists=not args.overwrite, overwrite=args.overwrite)
        mt = mt.filter_rows(hl.is_defined(ht[mt.row_key]))

        if args.overwrite or not hl.hadoop_exists(f'{get_ukb_grm_plink_path()}.bed'):
            hl.export_plink(mt, get_ukb_grm_plink_path())

        mt = get_filtered_mt()
        if args.overwrite or not hl.hadoop_exists(get_ukb_samples_file_path()):
            with hl.hadoop_open(get_ukb_samples_file_path(), 'w') as f:
                f.write('\n'.join(mt.ukbb_app_26041_id.collect()) + '\n')

    if args.create_gene_mapping_files:
        ht = hl.read_table(get_ukb_vep_path())
        gene_map_ht = create_gene_map_ht(ht)
        gene_map_ht.write(get_ukb_gene_map_ht_path(post_processed=False), args.overwrite)

        gene_map_ht = hl.read_table(get_ukb_gene_map_ht_path(post_processed=False))
        gene_map_ht = post_process_gene_map_ht(gene_map_ht)
        gene_map_ht.write(get_ukb_gene_map_ht_path(), args.overwrite)


if __name__ == '__main__':
    parser = argparse.ArgumentParser()

    parser.add_argument('--overwrite', help='Overwrite everything', action='store_true')
    parser.add_argument('--create_plink_file', help='Overwrite everything', action='store_true')
    parser.add_argument('--create_gene_mapping_files', help='Overwrite everything', action='store_true')
    parser.add_argument('--slack_channel', help='Send message to Slack channel/user', default='@konradjk')
    args = parser.parse_args()

    if args.slack_channel:
        try_slack(args.slack_channel, main, args)
    else:
        main(args)


