#!/usr/bin/env python3

__author__ = 'konradk'

from gnomad_hail import *
import ukbb_qc.resources as ukb
from ukb_exomes.resources import *


def main(args):
    if args.create_plink_file:
        # mt = ukb.get_ukbb_data('regeneron', raw=True, split=False, adj=True)
        # mt = hl.filter_intervals(mt, [hl.parse_locus_interval('chr1-chr22', 'GRCh38')])
        # mt = mt.filter_rows((hl.agg.fraction(hl.is_defined(mt.GT)) == 1) &
        #                     (hl.agg.call_stats(mt.GT, mt.alleles).AF[1] > 0.01))
        # mt = mt.checkpoint(common_high_callrate_path, _read_if_exists=True)
        # mt = hl.read_matrix_table(common_high_callrate_path)
        # mt = mt.filter_rows(hl.len(mt.alleles) == 2)
        # ht = hl.ld_prune(mt.GT, r2=0.1)
        # ht = ht.checkpoint(common_high_callrate_pruned_path)#, _read_if_exists=True)

        mt = hl.read_matrix_table(ukb.qc_mt_path('regeneron', ld_pruned=True))
        ht = hl.import_table(sample_mapping_file, key='eid_sample', delimiter=',')
        mt = mt.annotate_cols(exome_id=ht[mt.s.split('_')[1]].eid_26041)
        mt = mt.filter_cols(hl.is_defined(mt.exome_id)).key_cols_by('exome_id')

        if args.overwrite or not hl.hadoop_exists(f'{common_high_callrate_pruned_plink_path}.bed'):
            autosomes_mt = filter_to_autosomes(mt, reference='GRCh38')
            hl.export_plink(autosomes_mt, common_high_callrate_pruned_plink_path)

        if args.overwrite or not hl.hadoop_exists(f'{common_high_callrate_pruned_plink_path}.vcf.bgz'):
            export_mt = hl.filter_intervals(mt, [hl.parse_locus_interval('chr1', 'GRCh38')])
            hl.export_vcf(export_mt, f'{common_high_callrate_pruned_plink_path}.vcf.bgz')

            hl.export_vcf(mt, f'{common_high_callrate_pruned_plink_path}.vcf.bgz')

        if args.overwrite or not hl.hadoop_exists(f'{common_high_callrate_pruned_plink_path}.samples'):
            with hl.hadoop_open(f'{common_high_callrate_pruned_plink_path}.samples', 'w') as f:
                f.write('\n'.join(mt.exome_id.collect()) + '\n')


if __name__ == '__main__':
    parser = argparse.ArgumentParser()

    parser.add_argument('--overwrite', help='Overwrite everything', action='store_true')
    parser.add_argument('--create_plink_file', help='Overwrite everything', action='store_true')
    parser.add_argument('--slack_channel', help='Send message to Slack channel/user', default='@konradjk')
    args = parser.parse_args()

    if args.slack_channel:
        try_slack(args.slack_channel, main, args)
    else:
        main(args)


