#!/usr/bin/env python3

__author__ = 'konradk'

import csv
import subprocess

from gnomad_hail import *
from ukb_exomes import *
from scipy.cluster.hierarchy import dendrogram, linkage
from matplotlib import pyplot as plt


def pre_process_data_dictionary():
    local_pheno_description_path = '/tmp/Data_Dictionary_Showcase.csv'
    local_pheno_description_out_path = '/tmp/Data_Dictionary_Showcase.tsv'
    hl.hadoop_copy(pheno_description_raw_path, f'file://{local_pheno_description_path}')
    with open(local_pheno_description_path) as f, open(local_pheno_description_out_path, 'w') as g:
        reader = csv.reader(f)
        for line in reader:
            g.write('\t'.join(line) + '\n')
    hl.hadoop_copy(f'file://{local_pheno_description_out_path}', pheno_description_path)


def get_codings(overwrite: bool = False):
    root = '/tmp/PHESANT'
    # if subprocess.check_call(['git', 'clone', 'https://github.com/konradjk/PHESANT.git', '-b', 'more_codings', root]):
    if subprocess.check_call(['git', 'clone', 'https://github.com/astheeggeggs/PHESANT.git', root]):
        raise Exception('Could not clone repo')
    hts = []
    coding_dir = f'{root}/WAS/codings'
    for coding_file in os.listdir(f'{coding_dir}'):
        hl.hadoop_copy(f'file://{coding_dir}/{coding_file}', f'{coding_dir}/{coding_file}')
        ht = hl.import_table(f'{coding_dir}/{coding_file}')
        if 'node_id' not in ht.row:
            ht = ht.annotate(node_id=hl.null(hl.tstr), parent_id=hl.null(hl.tstr), selectable=hl.null(hl.tstr))
        ht = ht.annotate(coding_id=hl.int(coding_file.split('.')[0].replace('coding', '')))
        hts.append(ht)
    full_ht = hts[0].union(*hts[1:]).key_by('coding_id', 'coding')
    full_ht.repartition(10).write(coding_ht_path, overwrite=overwrite)


def pheno_ht_to_mt(pheno_ht: hl.Table, data_type: str, sex: str, overwrite: bool = False):
    if data_type == 'categorical':
        category_type = hl.int
        filter_type = {hl.tbool}
        value_type = hl.bool
    else:
        category_type = hl.str
        filter_type = {hl.tint, hl.tfloat}
        value_type = hl.float

    pheno_ht = pheno_ht.select(**{x: value_type(v) for x, v in pheno_ht.row_value.items() if v.dtype in filter_type})

    mt = pheno_ht.to_matrix_table_row_major(
        columns=list(pheno_ht.row_value), entry_field_name='value', col_field_name='phesant_pheno'
    )
    mt = mt.key_cols_by(
        pheno=hl.int(mt.phesant_pheno.split('_')[0]),
        coding=hl.cond(hl.len(mt.phesant_pheno.split('_')) > 1,
                       category_type(mt.phesant_pheno.split('_')[1]),
                       category_type(''))  # NB: This would error for categoricals if there were any missing a value
    )
    mt.write(get_ukb_pheno_mt_path(data_type, sex), overwrite)


def get_missing_codings(ht):
    data = ht.aggregate(hl.agg.filter(hl.is_missing(ht.meaning), hl.agg.collect_as_set(ht.coding_id)))
    if False and data != {1965}:
        import requests
        print(f'Missing: {data}')
        for coding in data:
            r = requests.post(url='http://biobank.ndph.ox.ac.uk/showcase/codown.cgi', data={'id': coding})
            with open(f'/tmp/coding{coding}.tsv', 'w') as f:
                f.write(r.text)

        data = ht.aggregate(hl.agg.filter(hl.is_missing(ht.meaning),
                                          hl.agg.collect_as_set((ht.coding_id, ht.coding))))


def get_phesant_reassignments(phesant_summary):
    phesant_summary = phesant_summary.annotate(
        reassign=phesant_summary['PHESANT.reassignments'].split(' ')[1].split('\|').map(lambda x: x.split('=')))
    ht = phesant_summary.explode('reassign')
    ht = ht.filter(ht.reassign[1] != 'NA')
    ht = ht.transmute(reassign_from=hl.int(ht.reassign[0]), reassign_to=hl.int(ht.reassign[1]))
    ht = ht.key_by(
        pheno=hl.int(ht.FieldID.split('_')[0]),
        coding=hl.or_missing(hl.len(ht.FieldID.split('_')) > 1, hl.int(ht.FieldID.split('_')[1]))
    )
    ht = ht.filter(ht.reassign_to == ht.coding)
    return ht


def add_coding_information(mt: hl.MatrixTable, coding_ht: hl.Table, sex: str) -> hl.MatrixTable:
    mt = mt.annotate_cols(**coding_ht[(mt.coding_id, hl.str(mt.coding))])
    get_missing_codings(mt.cols())
    phesant_summary = hl.import_table(get_ukb_phesant_summary_tsv_path(sex), impute=True, missing='', key='FieldID')
    phesant_reassign = get_phesant_reassignments(phesant_summary)
    mt = mt.annotate_cols(recoding=hl.or_missing(
        hl.is_missing(mt.meaning), phesant_reassign[mt.col_key].reassign_from
    ))
    mt = mt.annotate_cols(**hl.cond(hl.is_defined(mt.meaning),
                                    hl.struct(**{x: mt[x] for x in list(coding_ht.row_value)}),
                                    coding_ht[(mt.coding_id, hl.str(mt.recoding))]),
                          )
    return mt


def combine_datasets(data_type: str = 'categorical'):
    both_mt = hl.read_matrix_table(get_ukb_pheno_mt_path(data_type, 'both_sexes_no_sex_specific'))
    female_mt = hl.read_matrix_table(get_ukb_pheno_mt_path(data_type, 'females'))
    male_mt = hl.read_matrix_table(get_ukb_pheno_mt_path(data_type, 'males'))

    description_ht = hl.import_table(pheno_description_path, impute=True, missing='', key='FieldID')
    description_ht = description_ht.transmute(coding_id=description_ht.Coding)

    both_mt = both_mt.annotate_cols(**description_ht[both_mt.pheno])
    female_mt = female_mt.annotate_cols(**description_ht[female_mt.pheno])
    male_mt = male_mt.annotate_cols(**description_ht[male_mt.pheno])

    if data_type == 'categorical':
        coding_ht = hl.read_table(coding_ht_path)
        both_mt = add_coding_information(both_mt, coding_ht, 'both_sexes_no_sex_specific')
        female_mt = add_coding_information(female_mt, coding_ht, 'females')
        male_mt = add_coding_information(male_mt, coding_ht, 'males')

    mt = hl.experimental.full_outer_join_mt(both_mt, female_mt)
    mt = mt.select_entries(
        both_sexes=mt.left_entry.value,
        females=mt.right_entry.value,
    ).drop('left_row', 'right_row')
    # ht = mt.cols().persist()
    # assert ht.all(ht.left_col == ht.right_col)
    mt = mt.select_cols(both_sexes_pheno=mt.left_col.drop(*mt.col_key), females_pheno=mt.right_col.drop(*mt.col_key))

    mt = hl.experimental.full_outer_join_mt(mt, male_mt)
    mt = mt.select_entries(**mt.left_entry, males=mt.right_entry.value).drop('left_row', 'right_row')
    mt = mt.select_cols(**mt.left_col.drop(*mt.col_key), males_pheno=mt.right_col.drop(*mt.col_key))
    return mt


def load_icd_data(overwrite: bool = False):
    code_locations = {
        'primary_codes': '41202',
        'secondary_codes': '41204',
        'external_codes': '41201',
        'cause_of_death_codes': '40001'
    }
    ht = hl.import_table(pre_phesant_data_path, impute=True, min_partitions=100, missing='', key='userId')
    ht.write('gs://ukbb-pharma-exome-analysis/tmp/pre_phesant.ht', args.overwrite)
    ht = hl.read_table('gs://ukbb-pharma-exome-analysis/tmp/pre_phesant.ht')
    all_phenos = list(ht.row_value)
    ht = ht.select(
        **{code: [ht[x] for x in all_phenos if x.startswith(f'x{loc}')] for code, loc in code_locations.items()})
    ht = ht.annotate(**{code: ht[code].filter(lambda x: hl.is_defined(x)) for code in code_locations})
    all_codes = hl.sorted(hl.array(ht.aggregate(
        hl.agg.explode(lambda c: hl.agg.collect_as_set(c),
                       hl.flatmap(lambda x: x, [ht[code] for code in code_locations])),
        _localize=False)))
    ht = ht.select(
        bool_codes=all_codes.map(lambda x: hl.struct(**{code: ht[code].contains(x) for code in code_locations})))
    ht = ht.annotate_globals(all_codes=all_codes.map(lambda x: hl.struct(icd_code=x)))
    mt = ht._unlocalize_entries('bool_codes', 'all_codes', ['icd_code'])
    coding_ht = hl.import_table(icd_codings_path, impute=True, key='coding')
    mt = mt.annotate_cols(**coding_ht[mt.col_key]).annotate_globals(code_locations=code_locations)
    mt.write(get_ukb_pheno_mt_path('icd', 'full'), overwrite)


def export_pca_results(data_type):
    mt = hl.read_matrix_table(get_ukb_pheno_mt_path(data_type, 'full'))
    ht = hl.read_table(f'gs://ukbb-pharma-exome-analysis/tmp/{data_type}_pca.ht')
    if data_type == 'icd':
        mt = mt.annotate_cols(fraction_called=hl.agg.fraction(mt.primary_codes))
        ht = ht.annotate(fraction_called=mt.cols()[ht.key].fraction_called)
    else:
        mt = mt.annotate_cols(fraction_called=hl.agg.fraction(hl.is_defined(mt.both_sexes)))
        ht = ht.annotate(**mt.cols()[ht.key].both_sexes_pheno, fraction_called=mt.cols()[ht.key].fraction_called)
    ht = ht.transmute(**{'PC' + str(i + 1): ht.scores[i] for i in range(10)})
    ht.export(f'gs://ukbb-pharma-exome-analysis/tmp/{data_type}_pca.txt.bgz')


def main(args):
    sexes = ('both_sexes_no_sex_specific', 'females', 'males')
    data_types = ('categorical', 'continuous')
    if args.load_data:
        pre_process_data_dictionary()
        get_codings(args.overwrite)
        load_icd_data(args.overwrite)

        for sex in sexes:
            pheno_ht = hl.import_table(get_ukb_source_tsv_path(sex), impute=True, min_partitions=100, missing='', key='userID')
            pheno_ht.write(get_ukb_pheno_ht_path(sex), overwrite=args.overwrite)

            pheno_ht = hl.read_table(get_ukb_pheno_ht_path(sex))
            for data_type in data_types:
                pheno_ht_to_mt(pheno_ht, data_type, sex, args.overwrite)

            # All codings for categorical phenotypes are bools
            # pheno_types = [(int(x[0].split('_')[0]), x[1].dtype)
            #                for x in list(pheno_ht.row.items())
            #                if len(x[0].split('_')) > 1
            #                and x[0].split('_')[1] != 'raw' and x[0].split('_')[1] != 'irnt']
            # Counter(map(lambda x: x[1], set(pheno_types)))  # Counter({dtype('bool'): 74})
            # All codings for raw phenotypes are ints or floats
            # pheno_types = [(int(x[0].split('_')[0]), x[1].dtype)
            #                for x in list(pheno_ht.row.items())
            #                if len(x[0].split('_')) > 1 and x[0].split('_')[1] == 'raw']
            # Counter(map(lambda x: x[1], set(pheno_types)))  # Counter({dtype('float64'): 259, dtype('int32'): 67})
            # All codings for irnt phenotypes are floats
            # pheno_types = [(int(x[0].split('_')[0]), x[1].dtype)
            #                for x in list(pheno_ht.row.items())
            #                if len(x[0].split('_')) > 1 and x[0].split('_')[1] == 'irnt']
            # Counter(map(lambda x: x[1], set(pheno_types)))  # Counter({dtype('float64'): 326})

    if args.combine_data:
        for data_type in data_types:
            mt = combine_datasets(data_type)
            mt.write(get_ukb_pheno_mt_path(data_type, 'full'), args.overwrite)

    if args.pca_missingness:
        for data_type in data_types:
            mt = hl.read_matrix_table(get_ukb_pheno_mt_path(data_type, 'full'))
            mt = mt.filter_cols(hl.is_defined(mt.both_sexes_pheno))
            mt = mt.annotate_entries(pheno_defined=hl.is_defined(mt.both_sexes))
            mt = mt.annotate_rows(mean_pheno_defined=hl.agg.mean(mt.pheno_defined))
            mt = mt.annotate_entries(pheno_defined_norm=mt.pheno_defined - mt.mean_pheno_defined)
            eigenvalues, scores, _ = hl.pca(mt.pheno_defined_norm, compute_loadings=False)
            print('Eigenvalues: {}'.format(eigenvalues))
            scores.write(f'gs://ukbb-pharma-exome-analysis/tmp/{data_type}_pca.ht', args.overwrite)
            export_pca_results(data_type)

            mt = hl.read_matrix_table(get_ukb_pheno_mt_path(data_type, 'full'))
            mt = mt.filter_cols(hl.is_defined(mt.both_sexes_pheno))
            bm = hl.linalg.BlockMatrix.from_entry_expr(hl.is_defined(mt.both_sexes))
            # bm = hl.linalg.BlockMatrix.from_entry_expr(hl.is_defined(mt.both_sexes),
            #                                            mean_impute=True, center=True, normalize=True, axis='cols')
            corr = bm.T @ bm
            df = corr.to_numpy()
            # See gs://ukbb-pharma-exome-analysis/missingness_correlation.ipynb

        data_type = 'icd'
        mt = hl.read_matrix_table(get_ukb_pheno_mt_path(data_type, 'full'))
        mt = mt.annotate_rows(mean_primary_codes=hl.agg.mean(mt.primary_codes))
        mt = mt.annotate_entries(primary_codes_norm=hl.int(mt.primary_codes) - mt.mean_primary_codes)
        eigenvalues, scores, _ = hl.pca(mt.primary_codes_norm, compute_loadings=False)
        print('Eigenvalues: {}'.format(eigenvalues))
        scores.write(f'gs://ukbb-pharma-exome-analysis/tmp/{data_type}_pca.ht', args.overwrite)
        export_pca_results(data_type)


if __name__ == '__main__':
    parser = argparse.ArgumentParser()

    parser.add_argument('--overwrite', help='Overwrite everything', action='store_true')
    parser.add_argument('--load_data', help='Load data', action='store_true')
    parser.add_argument('--combine_data', help='Combine datasets', action='store_true')
    parser.add_argument('--pca_missingness', help='Run missingness PCA', action='store_true')
    parser.add_argument('--slack_channel', help='Send message to Slack channel/user', default='@konradjk')
    args = parser.parse_args()

    if args.slack_channel:
        try_slack(args.slack_channel, main, args)
    else:
        main(args)


