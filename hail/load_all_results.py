#!/usr/bin/env python3

__author__ = 'konradk'

import argparse
from collections import defaultdict
from pprint import pprint
from ukb_common import *
from ukb_exomes import *

pheno_list = """
continuous-30600-irnt
continuous-30600-raw
continuous-30610-irnt
continuous-30610-raw
continuous-30620-irnt
continuous-30620-raw
continuous-30630-irnt
continuous-30630-raw
continuous-30640-irnt
continuous-30640-raw
continuous-30650-irnt
continuous-30650-raw
continuous-30660-irnt
continuous-30660-raw
continuous-30670-irnt
continuous-30670-raw
continuous-30680-irnt
continuous-30680-raw
continuous-30690-irnt
continuous-30690-raw
continuous-30700-irnt
continuous-30700-raw
continuous-30710-irnt
continuous-30710-raw
continuous-30720-irnt
continuous-30720-raw
continuous-30730-irnt
continuous-30730-raw
continuous-30740-irnt
continuous-30740-raw
continuous-30750-irnt
continuous-30750-raw
continuous-30760-irnt
continuous-30760-raw
continuous-30770-irnt
continuous-30770-raw
continuous-30780-irnt
continuous-30780-raw
continuous-30790-irnt
continuous-30790-raw
continuous-30800-irnt
continuous-30800-raw
continuous-30810-irnt
continuous-30810-raw
continuous-30820-irnt
continuous-30820-raw
continuous-30830-irnt
continuous-30830-raw
continuous-30840-irnt
continuous-30840-raw
continuous-30850-irnt
continuous-30850-raw
continuous-30860-irnt
continuous-30860-raw
continuous-30870-irnt
continuous-30870-raw
continuous-30880-irnt
continuous-30880-raw
continuous-30890-irnt
continuous-30890-raw""".strip().split('\n')


def main(args):
    hl.init(default_reference='GRCh38')
    col_keys = ['pheno', 'coding', 'trait_type']
    all_phenos_dir = hl.hadoop_ls(get_results_dir())
    if args.biomarkers_only:
        all_phenos_dir = [x for x in all_phenos_dir if any([y in x['path'] for y in pheno_list])]
        temp_directory = f'{temp_bucket}/biomarkers'
    else:
        temp_directory = temp_bucket

    pheno_ht = hl.read_matrix_table(get_ukb_pheno_mt_path()).cols()
    pheno_dict = hl.dict(pheno_ht.aggregate(hl.agg.collect(
        (hl.struct(pheno=pheno_ht.pheno, coding=pheno_ht.coding, trait_type=pheno_ht.data_type), pheno_ht.row_value)),
        _localize=False))

    if args.load_gene_results:
        all_gene_outputs = get_files_in_parent_directory(all_phenos_dir, 'gene_results.ht')

        print(f'Got {len(all_gene_outputs)} HT paths...')
        all_hts = list(map(lambda x: unify_saige_burden_ht_schema(hl.read_table(x)), all_gene_outputs))
        union_ht(all_hts, col_keys, pheno_dict, temp_directory + '/gene_ht', '_read_if_exists').write(
            get_results_mt_path(extension='ht', biomarkers=args.biomarkers_only), overwrite=args.overwrite)

        row_keys = ['gene_id', 'gene_symbol', 'annotation']
        mt = join_pheno_hts_to_mt(all_hts, row_keys, col_keys, temp_directory + '/gene', inner_mode=inner_mode)
        mt = mt.annotate_cols(**pheno_dict[mt.col_key])
        mt = pull_out_fields_from_entries(mt, ['interval', 'markerIDs', 'markerAFs', 'total_variants'] +
                                          [f'Nmarker_MACCate_{i}' for i in range(1, 9)], 'rows')
        mt = mt.naive_coalesce(1000).checkpoint(
            get_results_mt_path(biomarkers=args.biomarkers_only), overwrite=args.overwrite, _read_if_exists=not args.overwrite)
        mt.cols().write(get_final_pheno_info_ht_path(), overwrite=args.overwrite)

    if args.load_variant_results:
        all_variant_outputs = get_files_in_parent_directory(all_phenos_dir, 'variant_results.ht')

        all_hts = list(map(lambda x: unify_saige_ht_variant_schema(hl.read_table(x)), all_variant_outputs))
        print(f'Got {len(all_variant_outputs)} HT paths...')
        union_ht(all_hts, col_keys, pheno_dict, temp_directory + '/variant_ht', inner_mode='_read_if_exists'
                 ).naive_coalesce(20000).write(get_results_mt_path('variant', extension='ht', biomarkers=args.biomarkers_only), overwrite=args.overwrite)

        row_keys = ['locus', 'alleles']
        mt = join_pheno_hts_to_mt(all_hts, row_keys, col_keys, temp_directory + '/variant', inner_mode, repartition_final=20000)
        mt = mt.annotate_cols(**pheno_dict[mt.col_key])
        mt = pull_out_fields_from_entries(mt, ['markerID', 'AC', 'AF', 'N', 'gene', 'annotation'], 'rows')
        mt = pull_out_fields_from_entries(mt, ['N.Cases', 'N.Controls'], 'cols').drop('Tstat', 'varT', 'varTstar')
        mt.write(get_results_mt_path('variant', biomarkers=args.biomarkers_only), overwrite=args.overwrite)

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


if __name__ == '__main__':
    parser = argparse.ArgumentParser()

    parser.add_argument('--overwrite', help='Overwrite everything', action='store_true')
    parser.add_argument('--skip_per_pheno', help='Overwrite everything', action='store_true')
    parser.add_argument('--force_overwrite_per_pheno', help='Overwrite everything', action='store_true')
    parser.add_argument('--dry_run', help='Overwrite everything', action='store_true')
    parser.add_argument('--biomarkers_only', help='Overwrite everything', action='store_true')
    parser.add_argument('--load_gene_results', help='Overwrite everything', action='store_true')
    parser.add_argument('--load_variant_results', help='Overwrite everything', action='store_true')
    parser.add_argument('--find_errors', help='Overwrite everything', action='store_true')
    parser.add_argument('--slack_channel', help='Send message to Slack channel/user', default='@konradjk')
    args = parser.parse_args()

    if args.slack_channel:
        try_slack(args.slack_channel, main, args)
    else:
        main(args)