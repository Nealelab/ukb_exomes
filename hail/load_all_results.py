#!/usr/bin/env python3

__author__ = 'konradk'

from ukb_common import *
from ukb_exomes import *

bucket = 'gs://ukbb-pharma-exome-analysis'
temp_bucket = 'gs://ukbb-pharma-exome-analysis-temp'
results_dir = 'gs://ukbb-pharma-exome-analysis/data/result/'
final_gene_results_ht = 'gs://ukbb-pharma-exome-analysis/data/results.ht'
final_pheno_info_ht = 'gs://ukbb-pharma-exome-analysis/data/pheno_info.ht'
final_variant_results_ht = 'gs://ukbb-pharma-exome-analysis/data/variant_results.ht'


def main(args):
    hl.init(default_reference='GRCh38')
    all_phenos = hl.hadoop_ls(results_dir)

    if args.load_gene_results:
        all_gene_outputs = []
        all_gene_mt_outputs = []
        to_run = 0
        for directory in all_phenos:
            # if 'E109' in directory['path'] or \
            #         'E119' in directory['path'] or \
            #         'K509' in directory['path'] or \
            #         'J459' in directory['path']: continue
            gene_results_ht_path = f'{directory["path"]}/gene_results.ht'
            gene_results_mt_path = f'{directory["path"]}/gene_results.mt'
            if args.skip_per_pheno:
                if hl.hadoop_exists(f'{gene_results_mt_path}/_SUCCESS'):
                    all_gene_outputs.append(gene_results_ht_path)
                    all_gene_mt_outputs.append(gene_results_mt_path)
                continue
            if args.force_overwrite_per_pheno or not hl.hadoop_exists(gene_results_mt_path):
                pheno = directory['path'].split('/')[-1].split('-')
                coding = '' if len(pheno) == 1 else pheno[1]
                pheno = pheno[0]
                to_run += 1
                if to_run % 10 == 0: print(to_run)
                if args.dry_run: continue
                load_gene_data(directory['path'], pheno, coding, get_ukb_gene_map_ht_path(False))
            all_gene_outputs.append(gene_results_ht_path)
            all_gene_mt_outputs.append(gene_results_mt_path)
        if args.dry_run:
            print(f'Would run {to_run} tasks...')
            sys.exit(0)

        print(f'Got {len(all_gene_outputs)} HT paths...')
        pheno_ht = hl.read_matrix_table(get_ukb_pheno_mt_path()).cols()
        all_hts = list(map(lambda x: hl.read_table(x).select_globals(), all_gene_outputs))
        print(f'Unioning {len(all_hts)} HTs...')
        ht = all_hts[0].union(*all_hts[1:], unify=True)
        ht = ht.annotate(**pheno_ht[hl.struct(pheno=ht.pheno, coding=ht.coding)])
        ht.write(final_gene_results_ht, overwrite=args.overwrite)

        all_mts = list(map(lambda x: hl.read_matrix_table(x), all_gene_mt_outputs))
        print(f'Read {len(all_mts)} MTs...')
        mt = union_mts_by_tree(all_mts, temp_bucket + '/gene')
        print(f'Unioned MTs...')
        mt = mt.annotate_cols(**pheno_ht[hl.struct(pheno=mt.pheno, coding=mt.coding)])
        mt = mt.checkpoint(final_gene_results_ht.replace('.ht', '.mt'), overwrite=args.overwrite, _read_if_exists=not args.overwrite)
        mt.cols().write(final_pheno_info_ht, overwrite=args.overwrite)

    if args.load_variant_results:
        all_variant_outputs = []
        all_variant_mt_outputs = []
        to_run = 0
        for directory in all_phenos:
            # if 'K519' not in directory['path']: continue
            variant_results_ht_path = f'{directory["path"]}/variant_results.ht'
            variant_results_mt_path = f'{directory["path"]}/variant_results.mt'
            if args.skip_per_pheno:
                if hl.hadoop_exists(f'{variant_results_mt_path}/_SUCCESS'):
                    all_variant_outputs.append(variant_results_ht_path)
                    all_variant_mt_outputs.append(variant_results_mt_path)
                continue
            if args.force_overwrite_per_pheno or not hl.hadoop_exists(variant_results_mt_path):
                pheno = directory['path'].split('/')[-1].split('-')
                coding = '' if len(pheno) == 1 else pheno[1]
                pheno = pheno[0]
                to_run += 1
                if to_run % 10 == 0: print(to_run)
                if args.dry_run: continue
                print(f'Loading {directory["path"]}')
                load_variant_data(directory["path"], pheno, coding, get_ukb_vep_path(), overwrite=args.overwrite)
            all_variant_outputs.append(variant_results_ht_path)
            all_variant_mt_outputs.append(variant_results_mt_path)
        if args.dry_run:
            print(f'Would run {to_run} tasks...')
            sys.exit(0)

        pheno_ht = hl.read_matrix_table(get_ukb_pheno_mt_path()).cols()
        all_hts = list(map(lambda x: hl.read_table(x), all_variant_outputs))
        ht = all_hts[0].union(*all_hts[1:], unify=True)
        ht = ht.annotate(**pheno_ht[hl.struct(pheno=ht.pheno, coding=ht.coding)])
        ht.write(final_variant_results_ht, overwrite=args.overwrite)

        all_mts = list(map(lambda x: hl.read_matrix_table(x), all_variant_mt_outputs))
        mt = union_mts_by_tree(all_mts, temp_bucket + '/variant')
        mt = mt.annotate_cols(**pheno_ht[hl.struct(pheno=mt.pheno, coding=mt.coding)])
        mt = mt.checkpoint(final_variant_results_ht.replace('.ht', '.mt').replace(bucket, temp_bucket),
                           overwrite=args.overwrite, _read_if_exists=not args.overwrite)
        mt.repartition(5000, shuffle=False).write(final_variant_results_ht.replace('.ht', '.mt'), overwrite=args.overwrite)

    if args.find_errors:
        for directory in all_phenos[:10]:
            all_files = hl.hadoop_ls(directory['path'])
            log_files = [x['path'] for x in all_files if x['path'].endswith('.log')]
            hl.grep('[Ee]rror', log_files)


if __name__ == '__main__':
    parser = argparse.ArgumentParser()

    parser.add_argument('--overwrite', help='Overwrite everything', action='store_true')
    parser.add_argument('--skip_per_pheno', help='Overwrite everything', action='store_true')
    parser.add_argument('--force_overwrite_per_pheno', help='Overwrite everything', action='store_true')
    parser.add_argument('--dry_run', help='Overwrite everything', action='store_true')
    parser.add_argument('--load_gene_results', help='Overwrite everything', action='store_true')
    parser.add_argument('--load_variant_results', help='Overwrite everything', action='store_true')
    parser.add_argument('--find_errors', help='Overwrite everything', action='store_true')
    parser.add_argument('--slack_channel', help='Send message to Slack channel/user', default='@konradjk')
    args = parser.parse_args()

    if args.slack_channel:
        try_slack(args.slack_channel, main, args)
    else:
        main(args)