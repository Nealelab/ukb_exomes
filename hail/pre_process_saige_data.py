#!/usr/bin/env python3

__author__ = 'konradk'

import argparse
from gnomad.utils.vep import process_consequences
from gnomad.utils import slack
from ukb_common import *
from ukb_exomes import *


def count_variants():
    from gnomad.utils.vep import process_consequences
    ht = hl.read_table(get_ukb_vep_path())
    qual_ht = hl.read_table(get_ukb_exomes_qual_ht_path(CURRENT_TRANCHE))
    qual_ht = qual_ht.filter(hl.len(qual_ht.filters) == 0)
    ht = ht.filter(hl.is_defined(qual_ht[ht.key]))

    ht = process_consequences(ht)
    ht = ht.explode(ht.vep.worst_csq_by_gene_canonical)
    ht = ht.annotate(
        variant_id=ht.locus.contig + ':' + hl.str(ht.locus.position) + '_' + ht.alleles[0] + '/' + ht.alleles[1],
        annotation=annotation_case_builder(ht.vep.worst_csq_by_gene_canonical))
    ht = ht.filter(hl.literal({'pLoF', 'LC', 'missense', 'synonymous'}).contains(ht.annotation))
    print(ht.count())


def main(args):
    hl.init(default_reference='GRCh38')

    if args.create_plink_file:
        call_stats_ht = hl.read_table(ukb.var_annotations_ht_path('ukb_freq', *TRANCHE_DATA[CURRENT_TRANCHE]))
        freq_index = hl.zip_with_index(call_stats_ht.freq_meta).find(lambda f: (hl.len(f[1]) == 1) & (f[1].get('group') == 'cohort'))[0]
        call_stats_ht = call_stats_ht.annotate(call_stats=call_stats_ht.freq[freq_index])
        hq_samples = call_stats_ht.aggregate(hl.agg.max(call_stats_ht.call_stats.AN) / 2)
        call_stats_ht = filter_ht_for_plink(call_stats_ht, hq_samples)
        call_stats_ht = call_stats_ht.naive_coalesce(1000).checkpoint(get_ukb_sites_ht_path(), _read_if_exists=not args.overwrite, overwrite=args.overwrite)
        # This needed SSDs and highmems, and had to revert to 7bb57d65fa28 for some reason
        mt = get_filtered_mt(adj=True, interval_filter=True, densify_sites_ht=call_stats_ht, semi_join_rows=True)
        mt = mt.annotate_rows(call_stats=call_stats_ht[mt.row_key].call_stats)
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
        qual_ht = hl.read_table(get_ukb_exomes_qual_ht_path(CURRENT_TRANCHE))
        qual_ht = qual_ht.filter(hl.len(qual_ht.filters) == 0)
        ht = ht.filter(hl.is_defined(qual_ht[ht.key]))
        gene_map_ht = create_gene_map_ht(ht)
        gene_map_ht.write(get_ukb_gene_map_ht_path(post_processed=False), args.overwrite)

        gene_map_ht = hl.read_table(get_ukb_gene_map_ht_path(post_processed=False))
        gene_map_ht = post_process_gene_map_ht(gene_map_ht)
        gene_map_ht.write(get_ukb_gene_map_ht_path(), args.overwrite)

    count_variants()

    ht = hl.read_table(get_ukb_vep_path())
    ht = process_consequences(ht)
    ht = ht.explode(ht.vep.worst_csq_by_gene_canonical)
    ht = ht.group_by(
        gene=ht.vep.worst_csq_by_gene_canonical.gene_symbol,
        consequence=annotation_case_builder(ht.vep.worst_csq_by_gene_canonical)
    ).partition_hint(100).aggregate(
        n_variants=hl.agg.count()
    )
    ht.write(get_ukb_gene_summary_ht_path(), args.overwrite)


if __name__ == '__main__':
    parser = argparse.ArgumentParser()

    parser.add_argument('--overwrite', help='Overwrite everything', action='store_true')
    parser.add_argument('--create_plink_file', help='Overwrite everything', action='store_true')
    parser.add_argument('--create_gene_mapping_files', help='Overwrite everything', action='store_true')
    parser.add_argument('--slack_channel', help='Send message to Slack channel/user', default='@konradjk')
    args = parser.parse_args()

    if args.slack_channel:
        from slack_token_pkg.slack_creds import slack_token
        with slack.slack_notifications(slack_token, args.slack_channel):
            main(args)
    else:
        main(args)


