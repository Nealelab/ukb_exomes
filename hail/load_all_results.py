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


def get_cohort_index(call_stats_ht, field='freq_meta'):
    freq_index = hl.zip_with_index(call_stats_ht[field]).find(lambda f: (hl.len(f[1]) == 1) &
                                                                           (f[1].get('group') == 'cohort'))[0]
    return freq_index.collect()[0]



def main(args):
    hl.init(default_reference='GRCh38', log='/merge_results.log')
    all_phenos_dir = hl.hadoop_ls(get_results_dir())

    pheno_ht = hl.read_matrix_table(get_ukb_pheno_mt_path()).cols()
    pheno_dict = create_broadcast_dict(pheno_ht)

    inner_mode = 'overwrite' if args.overwrite else '_read_if_exists'
    take_func = lambda x: hl.agg.take(x, 1)[0]

    if args.load_gene_results:
        all_gene_outputs = get_files_in_parent_directory(all_phenos_dir, 'gene_results.ht')

        print(f'Got {len(all_gene_outputs)} HT paths...')
        all_hts = list(map(lambda x: unify_saige_burden_ht_schema(hl.read_table(x)), all_gene_outputs))
        # union_ht(all_hts, col_keys, pheno_dict, temp_bucket + '/gene_ht', inner_mode=inner_mode).write(get_results_mt_path(extension='ht'), overwrite=args.overwrite)

        row_keys = ['gene_id', 'gene_symbol', 'annotation']
        mt = join_pheno_hts_to_mt(all_hts, row_keys, PHENO_KEY_FIELDS, temp_bucket + '/gene', inner_mode=inner_mode)
        mt = mt.annotate_cols(**pheno_dict[mt.col_key])
        mt = pull_out_fields_from_entries(mt.annotate_entries(total_variants_pheno=mt.total_variants),
                                          ['interval', 'markerIDs', 'markerAFs', 'total_variants'] +
                                          [f'Nmarker_MACCate_{i}' for i in range(1, 9)], 'rows',
                                          agg_funcs=[take_func] * 3 + [hl.agg.max] * 9)
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
        mt = pull_out_fields_from_entries(mt, ['markerID', 'gene', 'annotation'], 'rows',
                                          agg_funcs=[take_func] * 3)
        mt = pull_out_fields_from_entries(mt, ['N.Cases', 'N.Controls'], 'cols', agg_funcs=hl.agg.max)
        call_stats_ht = hl.read_table('gs://broad-ukbb/broad.freeze_7/variant_qc/variant_annotations/ukb_cohort_freq.ht')
        # call_stats_ht = hl.read_table(ukb.var_annotations_ht_path('ukb_freq', *TRANCHE_DATA[CURRENT_TRANCHE]))
        freq_index = get_cohort_index(call_stats_ht, field='cohort_freq_meta')
        mt = mt.annotate_rows(call_stats=call_stats_ht[mt.row_key].freq[freq_index])
        mt.write(get_results_mt_path('variant'), overwrite=args.overwrite)

    mt = hl.read_matrix_table(get_results_mt_path('variant') + '.bak')
    call_stats_ht = hl.read_table('gs://broad-ukbb/broad.freeze_7/variant_qc/variant_annotations/ukb_freq_fix.ht')
    mt = mt.annotate_rows(call_stats=call_stats_ht[mt.row_key].freq[0])
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
                           ).checkpoint(get_results_mt_path(random_phenos=True), overwrite=args.overwrite)
        join_random_phenos(all_random_phenos, 'variant_results.ht', ['locus', 'alleles'], inner_mode
                           ).checkpoint(get_results_mt_path('variant', random_phenos=True), overwrite=args.overwrite)

    if args.write_vep_file:
        for tranche in ('300k', '500k'):
            ht = hl.read_matrix_table(get_results_mt_path('variant', tranche)).rows()
            vep_ht = hl.read_table(var_annotations_ht_path('vep', *TRANCHE_DATA[tranche]))
            ht = ht.annotate(vep=vep_ht[ht.key].vep)
            ht.write(f'{public_bucket}/{tranche}/results/vep.ht', args.overwrite)

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
    parser.add_argument('--write_vep_file', help='Overwrite everything', action='store_true')
    parser.add_argument('--slack_channel', help='Send message to Slack channel/user', default='@konradjk')
    args = parser.parse_args()

    if args.slack_channel:
        from slack_token_pkg.slack_creds import slack_token
        with slack.slack_notifications(slack_token, args.slack_channel):
            main(args)
    else:
        main(args)