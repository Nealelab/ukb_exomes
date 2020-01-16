#!/usr/bin/env python3

__author__ = 'konradk'

from ukb_common import *
from ukb_exomes import *


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


