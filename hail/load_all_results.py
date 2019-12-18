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


def union_mts_by_tree(all_mts, temp_dir):
    chunk_size = int(len(all_mts) ** 0.5) + 1
    outer_mts = []
    for i in range(chunk_size):
        if i * chunk_size >= len(all_mts): break
        mt = all_mts[i * chunk_size]
        for j in range(1, chunk_size):
            if i * chunk_size + j >= len(all_mts): break
            # try:
            mt = mt.union_cols(all_mts[i * chunk_size + j], row_join_type='outer')
            # except:
            #     print(f'problem with {i * chunk_size} and {i * chunk_size + j}')
            #     mt.describe()
            #     all_mts[i * chunk_size + j].describe()
            #     sys.exit(1)
        outer_mts.append(mt.checkpoint(f'{temp_dir}/temp_output_{i}.mt', overwrite=True))
    mt = outer_mts[0]
    for next_mt in outer_mts[1:]:
        mt = mt.union_cols(next_mt, row_join_type='outer')
    return mt


def get_pheno_annotations():
    cont_ht = hl.read_matrix_table(get_ukb_pheno_mt_path('continuous', 'full')).cols()
    cont_ht = cont_ht.key_by(pheno=hl.str(cont_ht.pheno), coding=cont_ht.coding)
    cont_ht = cont_ht.select(meaning=hl.coalesce(cont_ht.both_sexes_pheno.Field,
                                                 cont_ht.females_pheno.Field,
                                                 cont_ht.males_pheno.Field),
                             path=hl.coalesce(cont_ht.both_sexes_pheno.Path,
                                              cont_ht.females_pheno.Path,
                                              cont_ht.males_pheno.Path)
                             )
    icd_ht = hl.read_table(icd_codings_ht_path)
    icd_ht = icd_ht.key_by(pheno=icd_ht.coding, coding='')
    icd_ht = icd_ht.select(meaning=icd_ht.short_meaning, path=icd_ht.meaning)
    pheno_ht = cont_ht.union(icd_ht)
    return pheno_ht


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
        pheno_ht = get_pheno_annotations()
        # raw_hts = list(map(lambda x: hl.read_table(x).select_globals(), all_gene_outputs))
        # print(f'Modifying HTs...')
        # all_hts = []
        # for ht in raw_hts:  # TODO: remove for 200k
        #     if ht.Pvalue.dtype == hl.tstr: ht = ht.annotate(Pvalue=hl.float(ht.Pvalue))
        #     if ht.Pvalue_Burden.dtype == hl.tstr: ht = ht.annotate(Pvalue_Burden=hl.float(ht.Pvalue_Burden))
        #     if ht.Pvalue_skato_NA.dtype == hl.tstr: ht = ht.annotate(Pvalue_skato_NA=hl.float(ht.Pvalue_skato_NA))
        #     if ht.Pvalue_burden_NA.dtype == hl.tstr: ht = ht.annotate(Pvalue_burden_NA=hl.float(ht.Pvalue_burden_NA))
        #     if ht.Pvalue_skat_NA.dtype == hl.tstr: ht = ht.annotate(Pvalue_skat_NA=hl.float(ht.Pvalue_skat_NA))
        #     all_hts.append(ht)
        # print(f'Unioning {len(raw_hts)} HTs...')
        # ht = all_hts[0].union(*all_hts[1:], unify=True)
        # ht = ht.annotate(**pheno_ht[hl.struct(pheno=ht.pheno, coding=ht.coding)])
        # ht = ht.checkpoint(final_gene_results_ht, overwrite=args.overwrite, _read_if_exists=not args.overwrite)

        raw_mts = list(map(lambda x: hl.read_matrix_table(x), all_gene_mt_outputs))
        print(f'Read {len(raw_mts)} MTs...')
        all_mts = []
        for ht in raw_mts:  # TODO: remove for 200k
            if ht.Pvalue.dtype == hl.tstr: ht = ht.annotate_entries(Pvalue=hl.float(ht.Pvalue))
            if ht.Pvalue_Burden.dtype == hl.tstr: ht = ht.annotate_entries(Pvalue_Burden=hl.float(ht.Pvalue_Burden))
            if ht.Pvalue_skato_NA.dtype == hl.tstr: ht = ht.annotate_entries(Pvalue_skato_NA=hl.float(ht.Pvalue_skato_NA))
            if ht.Pvalue_burden_NA.dtype == hl.tstr: ht = ht.annotate_entries(Pvalue_burden_NA=hl.float(ht.Pvalue_burden_NA))
            if ht.Pvalue_skat_NA.dtype == hl.tstr: ht = ht.annotate_entries(Pvalue_skat_NA=hl.float(ht.Pvalue_skat_NA))
            all_mts.append(ht)
        all_mts = [mt if 'n_cases' in list(mt.col) else  # TODO: remove for 200k
                   mt.annotate_cols(n_cases=hl.null(hl.tint32), n_controls=hl.null(hl.tint32))
                   for mt in all_mts]
        print(f'Modified MTs...')
        all_mts = [mt if 'pheno' in list(mt.col_key) else mt.key_cols_by('pheno', 'coding') for mt in all_mts]  # TODO: remove for 200k
        print(f'Modified MTs...')
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

        pheno_ht = get_pheno_annotations()
        raw_hts = list(map(lambda x: hl.read_table(x), all_variant_outputs))
        all_hts = []
        for ht in raw_hts:
            ht = ht.key_by('locus', 'alleles',
                           pheno=ht.pheno[0] if ht.pheno.dtype == hl.tarray(hl.tstr) else ht.pheno,  # TODO: remove for 200k
                           coding=ht.coding)
            if ht.SE.dtype == hl.tstr:  # TODO: remove for 200k
                ht = ht.annotate(SE=hl.null(hl.tfloat64))
            all_hts.append(ht)
        ht = all_hts[0].union(*all_hts[1:], unify=True)
        ht = ht.annotate(**pheno_ht[hl.struct(pheno=ht.pheno, coding=ht.coding)])
        ht = ht.checkpoint(final_variant_results_ht, overwrite=args.overwrite, _read_if_exists=not args.overwrite)

        raw_mts = list(map(lambda x: hl.read_matrix_table(x), all_variant_mt_outputs))
        all_mts = []
        for mt in raw_mts:
            mt = mt.key_cols_by(pheno=mt.pheno[0] if mt.pheno.dtype == hl.tarray(hl.tstr) else mt.pheno,  # TODO: remove for 200k
                                coding=mt.coding).key_rows_by('locus', 'alleles')
            if mt.SE.dtype == hl.tstr:  # TODO: remove for 200k
                mt = mt.annotate_entries(SE=hl.null(hl.tfloat64))
            if 'AF.Cases' not in list(mt.entry):  # TODO: remove for 200k
                fields = list(mt.entry)
                fields.remove('Pvalue')
                mt = mt.select_entries(*fields,
                    **{"AF.Cases": hl.null(hl.tfloat64), "AF.Controls": hl.null(hl.tfloat64),
                       "N.Cases": hl.null(hl.tint32), "N.Controls": hl.null(hl.tint32),
                       "Pvalue": mt.Pvalue})
            all_mts.append(mt)
        all_mts = [mt if 'n_cases' in list(mt.col) else  # TODO: remove for 200k
                   mt.annotate_cols(n_cases=hl.null(hl.tint32), n_controls=hl.null(hl.tint32))
                   for mt in all_mts]
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