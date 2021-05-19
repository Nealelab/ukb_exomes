#!/usr/bin/env python3

__author__ = 'konradk'

import argparse
from collections import defaultdict
from pprint import pprint
from gnomad.utils import slack
from ukbb_common import *
from ukb_exomes import *


def join_random_phenos(all_random_phenos, table_name, row_key, inner_mode):
    all_random_data = get_files_in_parent_directory(all_random_phenos, table_name)
    unify_func = unify_saige_ht_variant_schema if table_name == 'variant_results.ht' else unify_saige_burden_ht_schema
    all_hts = list(map(lambda x: unify_func(hl.read_table(x)), all_random_data))
    mt = join_pheno_hts_to_mt(all_hts, row_key, PHENO_KEY_FIELDS, f'{temp_bucket}/random_{table_name}', inner_mode=inner_mode)
    return mt.naive_coalesce(1000)


def main(args):
    hl.init(default_reference='GRCh38', log='/merge_results.log')
    all_phenos_dir = hl.hadoop_ls(get_results_dir())

    pheno_ht = hl.read_matrix_table(get_ukb_pheno_mt_path()).cols()
    pheno_dict = create_broadcast_dict(pheno_ht)

    inner_mode = 'overwrite' if args.overwrite else '_read_if_exists'

    if args.load_gene_results:
        all_gene_outputs = get_files_in_parent_directory(all_phenos_dir, 'gene_results.ht')

        print(f'Got {len(all_gene_outputs)} HT paths...')
        all_hts = list(map(lambda x: unify_saige_burden_ht_schema(hl.read_table(x)), all_gene_outputs))
        # union_ht(all_hts, col_keys, pheno_dict, temp_bucket + '/gene_ht', inner_mode=inner_mode).write(get_results_mt_path(extension='ht'), overwrite=args.overwrite)

        row_keys = ['gene_id', 'gene_symbol', 'annotation']
        mt = join_pheno_hts_to_mt(all_hts, row_keys, PHENO_KEY_FIELDS, temp_bucket + '/gene', inner_mode=inner_mode)
        mt = mt.annotate_cols(**pheno_dict[mt.col_key])
        mt = pull_out_fields_from_entries(mt, ['interval', 'markerIDs', 'markerAFs', 'total_variants'] +
                                          [f'Nmarker_MACCate_{i}' for i in range(1, 9)], 'rows')
        mt = mt.naive_coalesce(1000).checkpoint(get_results_mt_path(), overwrite=args.overwrite, _read_if_exists=not args.overwrite)
        mt.cols().write(get_final_pheno_info_ht_path(), overwrite=args.overwrite)

    if args.load_variant_results:
        all_variant_outputs = get_files_in_parent_directory(all_phenos_dir, 'variant_results.ht')
        all_hts = list(map(lambda x: unify_saige_ht_variant_schema(hl.read_table(x)), all_variant_outputs))
        print(f'Got {len(all_variant_outputs)} HT paths...')
        # union_ht(all_hts, col_keys, pheno_dict, temp_bucket + '/variant_ht', inner_mode='_read_if_exists'
        #          ).naive_coalesce(20000).write(get_results_mt_path('variant', extension='ht'), overwrite=args.overwrite)

        row_keys = ['locus', 'alleles']
        inner_mode = '_read_if_exists'
        mt = join_pheno_hts_to_mt(all_hts, row_keys, PHENO_KEY_FIELDS, temp_bucket + '/variant', repartition_final=20000, inner_mode=inner_mode)
        mt = mt.annotate_cols(**pheno_dict[mt.col_key])
        mt = pull_out_fields_from_entries(mt, ['markerID', 'AC', 'AF', 'gene', 'annotation'], 'rows')
        mt = pull_out_fields_from_entries(mt, ['N.Cases', 'N.Controls'], 'cols')
        mt.write(get_results_mt_path('variant'), overwrite=args.overwrite)

    if args.find_errors:
        pheno_dirs = [x for x in all_phenos_dir if x['is_dir']]
        all_errors = defaultdict(dict)
        for directory in pheno_dirs:
            _, pheno = directory['path'].rsplit('/', 1)
            all_files = hl.hadoop_ls(directory['path'])
            log_files = [x['path'] for x in all_files if x['path'].endswith('.log')]
            errors = hl.grep('[Ee]rror', log_files, show=False)
            for fname, errstrings in errors.items():
                for errstring in errstrings:
                    if pheno not in all_errors[errstring]:
                        all_errors[errstring][pheno] = []
                    all_errors[errstring][pheno].append(fname)
        pprint(dict(all_errors).keys())

    if args.find_unconverged:
        glmm_logs = hl.hadoop_ls(get_results_dir(location='null_glmm'))
        log_files = [x['path'] for x in glmm_logs if x['path'].endswith('.log')]
        errors = hl.grep('Large variance estimate observed in the iterations, model not converged...', log_files,
                         show=False, max_count=6000)
        pprint(errors)

    if args.load_random_phenos:
        all_random_phenos = [x for x in all_phenos_dir if 'random' in x['path']]
        join_random_phenos(all_random_phenos, 'gene_results.ht', ['gene_id', 'gene_symbol', 'annotation'], inner_mode
                           ).checkpoint(get_results_mt_path(random_phenos=True))
        join_random_phenos(all_random_phenos, 'variant_results.ht', ['locus', 'alleles'], inner_mode
                           ).checkpoint(get_results_mt_path('variant', random_phenos=True))


if __name__ == '__main__':
    parser = argparse.ArgumentParser()

    parser.add_argument('--overwrite', help='Overwrite everything', action='store_true')
    parser.add_argument('--skip_per_pheno', help='Overwrite everything', action='store_true')
    parser.add_argument('--force_overwrite_per_pheno', help='Overwrite everything', action='store_true')
    parser.add_argument('--dry_run', help='Overwrite everything', action='store_true')
    parser.add_argument('--load_gene_results', help='Overwrite everything', action='store_true')
    parser.add_argument('--load_variant_results', help='Overwrite everything', action='store_true')
    parser.add_argument('--find_unconverged', help='Overwrite everything', action='store_true')
    parser.add_argument('--load_random_phenos', help='Overwrite everything', action='store_true')
    parser.add_argument('--find_errors', help='Overwrite everything', action='store_true')
    parser.add_argument('--slack_channel', help='Send message to Slack channel/user', default='@konradjk')
    args = parser.parse_args()

    if args.slack_channel:
        from slack_token_pkg.slack_creds import slack_token
        with slack.slack_notifications(slack_token, args.slack_channel):
            main(args)
    else:
        main(args)