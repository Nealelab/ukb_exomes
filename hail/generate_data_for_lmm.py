#!/usr/bin/env python3

__author__ = 'konradk'

from gnomad_hail import *
import ukbb_qc.resources as ukb
from ukb_exomes.resources import *

DATA_SOURCE = 'broad'
MIN_CALL_RATE = 0.95
VARIANTS_PER_MAC_CATEGORY = 2000
VARIANTS_PER_MAF_CATEGORY = 10000


def create_gene_map_ht(ht, check_gene_contigs=False):
    ht = process_consequences(ht)
    ht = ht.explode(ht.vep.worst_csq_by_gene_canonical)
    ht = ht.annotate(
        variant_id=ht.locus.contig + ':' + hl.str(ht.locus.position) + '_' + ht.alleles[0] + '/' + ht.alleles[1],
        annotation=hl.case(missing_false=True)
            .when(ht.vep.worst_csq_by_gene_canonical.lof == 'HC', 'pLoF')
            .when(ht.vep.worst_csq_by_gene_canonical.lof == 'LC', 'LC')
            .when((ht.vep.worst_csq_by_gene_canonical.most_severe_consequence == 'missense_variant') |
                  (ht.vep.worst_csq_by_gene_canonical.most_severe_consequence == 'inframe_insertion') |
                  (ht.vep.worst_csq_by_gene_canonical.most_severe_consequence == 'inframe_deletion'),
                  'missense')
            .when(ht.vep.worst_csq_by_gene_canonical.most_severe_consequence == 'synonymous_variant',
                  'synonymous')
            .or_missing())
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


def prepare_mt_for_plink(mt: hl.MatrixTable, min_call_rate: float = MIN_CALL_RATE,
                         variants_per_mac_category: int = VARIANTS_PER_MAC_CATEGORY,
                         variants_per_maf_category: int = VARIANTS_PER_MAF_CATEGORY):
    mt = filter_to_autosomes(mt)
    mt = annotate_sample_data(mt)
    hq_samples = 100217  # mt.aggregate_cols(hl.agg.count_where(mt.meta.high_quality))
    call_stats_ht = hl.read_table(ukb.var_annotations_ht_path(DATA_SOURCE, ukb.CURRENT_FREEZE, 'call_stats'))
    mt = mt.annotate_rows(call_stats=call_stats_ht[mt.row_key].qc_callstats[0])
    mt = mt.filter_rows((mt.call_stats.AN >= hq_samples * 2 * min_call_rate) &
                        (mt.call_stats.AC[1] > 0))
    mt = mt.annotate_rows(
        mac_category=hl.case()
            .when(mt.call_stats.AC[1] <= 5, mt.call_stats.AC[1])
            .when(mt.call_stats.AC[1] <= 10, 10)
            .when(mt.call_stats.AC[1] <= 20, 20)
            .when(mt.call_stats.AF[1] <= 0.001, 0.001)
            .when(mt.call_stats.AF[1] <= 0.01, 0.01)
            .when(mt.call_stats.AF[1] <= 0.1, 0.1)
            .default(0.99)
    )
    # category_counter = mt.aggregate_rows(hl.agg.counter(mt.mac_category))
    # mt = mt.annotate_globals(category_counter=category_counter)
    mt = mt.annotate_globals(category_counter={5.0: 343544, 10.0: 792088, 20.0: 476789, 1.0: 5625273, 0.99: 74394,
                                               2.0: 1711929, 0.1: 65559, 0.01: 153733, 3.0: 868302, 4.0: 517775, 0.001: 583818})
    mt = mt.filter_rows(hl.rand_unif(0, 1) <
                        hl.cond(mt.mac_category >= 1, variants_per_mac_category, variants_per_maf_category) / mt.category_counter[mt.mac_category]
                        ).naive_coalesce(5000)
    return mt


def annotate_sample_data(mt):
    # meta_ht = hl.read_table(ukb.meta_ht_path(DATA_SOURCE, ukb.CURRENT_FREEZE))
    meta_ht = hl.read_table('gs://broad-ukbb/regeneron.freeze_4/sample_qc/meta_w_pop_adj.ht')
    mapping_ht = hl.import_table(sample_mapping_file, key='eid_sample', delimiter=',')
    mt = mt.annotate_cols(meta=meta_ht[mt.col_key],
                          exome_id=mapping_ht[mt.s.split('_')[1]].eid_26041)
    mt = mt.filter_cols(hl.is_defined(mt.exome_id) & ~mt.meta.is_filtered &
                        (mt.meta.hybrid_pop == '12')).key_cols_by('exome_id')  # TODO: fix populations
    return mt


def main(args):
    hl.init(default_reference='GRCh38')

    if args.vep:
        ht = hl.vep(get_ukb_exomes_mt().rows(), vep_config_path('GRCh38'))
        ht.write(get_ukb_vep_path(), args.overwrite)

        # gt_mt = hl.read_matrix_table()
        # ht = hl.vep(mt.rows(), vep_config_path('GRCh38'))
        # ht.write(genotypes_vep_path)

    if args.create_plink_file:
        mt = prepare_mt_for_plink(ukb.get_ukbb_data(DATA_SOURCE, ukb.CURRENT_FREEZE, adj=True))
        mt = mt.checkpoint(ukb_for_grm_mt_path, _read_if_exists=True)
        mt = mt.unfilter_entries()
        ht = hl.ld_prune(mt.GT, r2=0.1)
        ht = ht.annotate_globals(**mt.index_globals())
        ht = ht.checkpoint(ukb_for_grm_pruned_ht_path, _read_if_exists=False, overwrite=True)
        mt = mt.filter_rows(hl.is_defined(ht[mt.row_key]))

        if args.overwrite or not hl.hadoop_exists(f'{ukb_for_grm_plink_path}.bed'):
            hl.export_plink(mt, ukb_for_grm_plink_path)

        mt = ukb.get_ukbb_data(DATA_SOURCE, ukb.CURRENT_FREEZE, adj=True)
        mt = annotate_sample_data(mt)
        if args.overwrite or not hl.hadoop_exists(get_ukb_samples_file_path()):
            with hl.hadoop_open(get_ukb_samples_file_path(), 'w') as f:
                f.write('\n'.join(mt.exome_id.collect()) + '\n')

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
    parser.add_argument('--vep', help='Overwrite everything', action='store_true')
    parser.add_argument('--create_plink_file', help='Overwrite everything', action='store_true')
    parser.add_argument('--create_gene_mapping_files', help='Overwrite everything', action='store_true')
    parser.add_argument('--slack_channel', help='Send message to Slack channel/user', default='@konradjk')
    args = parser.parse_args()

    if args.slack_channel:
        try_slack(args.slack_channel, main, args)
    else:
        main(args)


