#!/usr/bin/env python3

__author__ = 'konradk'

import argparse
from pprint import pprint
from datetime import date
from gnomad.utils import slack
from ukb_common import *
from ukb_exomes import *

temp_bucket = 'gs://ukb-pharma-exome-analysis-temp'
pheno_directory = f'{bucket}/{CURRENT_TRANCHE}/pheno_combo_explore'
pairwise_correlation_ht_path = f'{pheno_directory}/pairwise_correlations.ht'

PRIORITIES = ('Abbvie_Priority', 'Biogen_Priority', 'Pfizer_Priority')


def read_covariate_data(pre_phesant_data_path):
    ht = hl.import_table(pre_phesant_data_path, impute=True, min_partitions=100, missing='', key='userId')
    columns = {
        'sex': 'x22001_0_0',
        'age': 'x21022_0_0'
    }
    columns.update(**{f'pc{i}': f'x22009_0_{i}' for i in range(1, 41)})
    ht = ht.select(*columns.values()).rename({v: k for k, v in columns.items()}).annotate_globals(coding_source=columns)
    meta_ht = hl.read_table(get_ukb_exomes_meta_ht_path(CURRENT_TRANCHE))
    batches = meta_ht.aggregate(hl.agg.counter(meta_ht.ukbb_meta.batch))
    meta_ht = meta_ht.key_by(userId=hl.int32(meta_ht.ukbb_meta.ukbb_app_26041_id))
    print(f'Found batches: {batches}')
    batch_field = meta_ht[ht.key].ukbb_meta.batch
    # One-hot encoding each batch (except 0)
    ht = ht.annotate(**{f'batch_{batch}': hl.int(batch_field == batch) for batch in batches if batch != '0' and batch is not None})
    return ht.annotate(age2=ht.age ** 2, age_sex=ht.age * ht.sex, age2_sex=ht.age ** 2 * ht.sex)


def extract_mt_by_type(ht: hl.Table, value_type = hl.expr.Int32Expression, trait_type = 'categorical', source = 'biogen'):
    columns = [x[0] for x in list(ht.row_value.items()) if isinstance(x[1], value_type)]
    ht2 = ht.select(values=[hl.struct(value=hl.float64(ht[x])) for x in columns]).annotate_globals(
        columns=hl.map(lambda x: hl.struct(trait_type=trait_type, phenocode=x, pheno_sex='both_sexes', coding=NULL_STR_KEY, modifier=source),
                       hl.literal(columns)))
    return ht2._unlocalize_entries('values', 'columns', PHENO_KEY_FIELDS)


def load_custom_data():
    custom_data_folder = 'gs://phenotype_pharma/custom_phenos'

    # Biogen
    ht = hl.import_table(f'{custom_data_folder}/new_phenos_BIOGEN.csv',
                         impute=True, delimiter=',', quote='"', missing="", key='UKBB_ID').rename({"UKBB_ID": 'userId'})
    ht = ht.annotate(Fluid_intelligence_score_BI=hl.float64(ht.Fluid_intelligence_score_BI),
                     Depressive_symptoms_BI=hl.float64(ht.Depressive_symptoms_BI),
                     Touchscreen_duration_BI=hl.float64(ht.Touchscreen_duration_BI))

    mt = extract_mt_by_type(ht, hl.expr.Int32Expression, trait_type='categorical')
    mt2 = extract_mt_by_type(ht, hl.expr.Float64Expression, trait_type='continuous')
    mt = mt.union_cols(mt2).annotate_cols(**{x: NULL_STR for x in PHENO_DESCRIPTION_FIELDS})

    # Abbvie
    ht = hl.import_table(f'{custom_data_folder}/Custom_Phenotypes_AbbVie_02042020_Merged.csv', impute=True, delimiter=' ', missing="NA", key='userId')
    ht2 = hl.import_table(f'{custom_data_folder}/Custom_Phenotypes_AbbVie_27042020_Merged.csv', impute=True, delimiter=' ', missing="NA", key='userId')
    ht = ht.annotate(**ht2[ht.key])
    mt_abbvie = extract_mt_by_type(ht, hl.expr.Int32Expression, trait_type='categorical', source='abbvie')
    abbvie_meta = hl.import_table(f'{custom_data_folder}/AbbVie_Custom_Phenotype_Category_Labels_Merged_05042020.csv', delimiter=',', impute=True, key='phenotype_id')
    abbvie_meta = abbvie_meta.select('description', description_more=abbvie_meta.secondary_category, coding_description=NULL_STR, category=abbvie_meta.primary_category)
    mt_abbvie = mt_abbvie.annotate_cols(**abbvie_meta[mt_abbvie.phenocode])

    # Pfizer
    ht = hl.import_table(f'{custom_data_folder}/pfe_phenotypes_07182020.txt', impute=True, delimiter='\t', missing="NA", key='f.eid', min_partitions=8)
    mt_pfizer = extract_mt_by_type(ht, hl.expr.Int32Expression, trait_type='categorical', source='pfizer')
    mt_pfizer2 = extract_mt_by_type(ht, hl.expr.Float64Expression, trait_type='continuous', source='pfizer')
    mt_pfizer = mt_pfizer.union_cols(mt_pfizer2).key_cols_by().select_cols('trait_type', 'phenocode')
    pfizer_meta = hl.import_table(f'{custom_data_folder}/Pfizer_custom_phenotypes_07182020.txt', delimiter='\t', impute=True)
    pfizer_meta = pfizer_meta.annotate(coding=NULL_STR_KEY, modifier=pfizer_meta.coding.lower()).select(*PHENO_KEY_FIELDS, *PHENO_DESCRIPTION_FIELDS).key_by('trait_type', 'phenocode')
    mt_pfizer = mt_pfizer.annotate_cols(**pfizer_meta[mt_pfizer.col.select('trait_type', 'phenocode')]).key_cols_by(*PHENO_KEY_FIELDS).rename({'f.eid': 'userId'})

    return mt.union_cols(mt_abbvie).union_cols(mt_pfizer)


def pre_process_first_occurrence():
    local_path = '/tmp/first_occurrence.csv'
    local_out_path = '/tmp/first_occurrence.tsv'
    hl.hadoop_copy(first_occcurrence_csv_path, f'file://{local_path}')
    with open(local_path, encoding='iso8859-1') as f, open(local_out_path, 'w') as g:
        reader = csv.reader(f)
        for line in reader:
            line = list(line)
            # if line[0].startswith('3246034') or line[0].startswith('4887441'):
            #     print(f'At: {line[0]}')
            #     print([(i, x) for i, x in enumerate(line) if '"' in x])
            #     print(line[8335])
            line.pop(8335)  # Offending field with misplaced quotes
            line.pop(8335)  # Offending field with misplaced quotes
            # if line[0].startswith('3246034') or line[0].startswith('4887441'):
            #     print(f'At: {line[0]}')
            #     print(line[8336])
            g.write('\t'.join(line) + '\n')
    hl.hadoop_copy(f'file://{local_out_path}', first_occcurrence_csv_path.replace('.csv', '.tsv'))


def read_random_phenos():
    random_mt = hl.import_matrix_table('gs://ukbb-pharma-exome-analysis/300k/random_phenos.tsv.bgz',
                                       row_fields={'userId': hl.tint32}, row_key='userId', entry_type=hl.tfloat64,
                                       min_partitions=100)
    return random_mt


def main(args):
    hl.init(log='/load_pheno.log')
    sexes = ('both_sexes_no_sex_specific', 'females', 'males')
    data_types = ('categorical', 'continuous') #, 'biomarkers')
    covid_data_wave = '01'
    curdate = date.today().strftime("%y%m%d")

    if args.load_data:
        # load_icd_data(get_pre_phesant_data_path(), icd_codings_ht_path, f'{temp_bucket}/{CURRENT_TRANCHE}').write(get_ukb_pheno_mt_path('icd', 'full'), args.overwrite)
        # read_covariate_data(get_pre_phesant_data_path()).write(get_ukb_covariates_ht_path(), args.overwrite)
        load_custom_data().write(get_ukb_pheno_mt_path('custom', 'full'), args.overwrite)

        # load_prescription_data(first_occcurrence_csv_path, prescription_mapping_path).write(
        #     get_ukb_pheno_mt_path('prescriptions'), args.overwrite)
        # pre_process_first_occurrence()
        # load_first_occurrence_data(first_occcurrence_tsv_path, first_occcurrence_tsv_path, '\t', quote=None).write(
        #     get_ukb_pheno_mt_path('icd_first_occurrence'), args.overwrite)

        # for sex in sexes:
        #     path = get_ukb_source_tsv_path(sex)
        #     pheno_ht = hl.import_table(path, impute=True, min_partitions=100, missing='', key='userId', force_bgz=path.endswith('.gz'))
        #     pheno_ht.write(get_ukb_pheno_ht_path(sex), overwrite=args.overwrite)
        #
        #     pheno_ht = hl.read_table(get_ukb_pheno_ht_path(sex))
        #     for data_type in data_types:
        #         pheno_ht_to_mt(pheno_ht, data_type).write(get_ukb_pheno_mt_path(data_type, sex), args.overwrite)
        #
        #     if CURRENT_TRANCHE == '100k':
        #         ht = hl.import_table(get_biomarker_source_tsv_path(sex), impute=True, min_partitions=100, missing='', key='userID')
        #         ht = ht.checkpoint(get_biomarker_ht_path(sex), overwrite=args.overwrite)
        #
        #         pheno_ht_to_mt(ht, 'continuous').write(get_ukb_pheno_mt_path('biomarkers', sex), overwrite=args.overwrite)

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
        # for data_type in data_types:
        #     mt = combine_datasets({sex: get_ukb_pheno_mt_path(data_type, sex) for sex in sexes},
        #                           {sex: get_ukb_phesant_summary_tsv_path(sex) for sex in sexes},
        #                           pheno_description_path, coding_ht_path, data_type)
        #     mt.write(get_ukb_pheno_mt_path(data_type, 'full'), args.overwrite)
        data_types = ('categorical', 'continuous', 'icd', 'custom', 'icd_first_occurrence')  # 'prescriptions',
        pheno_file_dict = {data_type: hl.read_matrix_table(get_ukb_pheno_mt_path(data_type, 'full')) for data_type in data_types}
        pheno_file_dict['random'] = read_random_phenos()
        cov_ht = hl.read_table(get_ukb_covariates_ht_path())
        mt = combine_pheno_files_multi_sex(pheno_file_dict, cov_ht)
        priority_ht = hl.import_table(pheno_description_with_priority_path, impute=True, missing='', key='FieldID'
                                      ).select(*PRIORITIES)
        mt = mt.annotate_cols(**priority_ht[mt.phenocode])
        mt = mt.annotate_cols(**{x: hl.if_else(mt.phenocode.startswith('random_'), 'h', mt[x]) for x in PRIORITIES})
        mt = mt.annotate_cols(
            score=2 * (hl.int(mt.Abbvie_Priority == 'h') + hl.int(mt.Biogen_Priority == 'h') + hl.int(mt.Pfizer_Priority == 'h')) +
                  hl.int(mt.Abbvie_Priority == 'm') + hl.int(mt.Biogen_Priority == 'm') + hl.int(mt.Pfizer_Priority == 'm'))
        mt.write(get_ukb_pheno_mt_path(), args.overwrite)

        mt = hl.read_matrix_table(get_ukb_pheno_mt_path())
        pprint(f'Defined scores: {mt.aggregate_cols(hl.agg.counter(hl.is_defined(mt.score)))}')
        mt.cols().export(f'{pheno_directory}/all_pheno_summary.txt.bgz')
        # make_correlation_ht(mt).write(pairwise_correlation_ht_path, args.overwrite)
        # hl.read_table(pairwise_correlation_ht_path).flatten().export(pairwise_correlation_ht_path.replace('.ht', '.txt.bgz'))

    if args.add_covid_wave:
        # hail load_phenotype_data --add_covid_wave 03 --overwrite
        # python saige_pan_ancestry.py --phenos .*COVID.*03.*
        ht = load_dob_ht(first_occcurrence_tsv_path, 'eid', year_field='34-0.0', month_field='52-0.0',
                         recruitment_center_field='54-0.0', quote=None)
        ht = ht.checkpoint(f'{bucket}/misc/covid_test/basic_dob.ht', _read_if_exists=True)
        mt = load_covid_data(ht, get_covid_data_path(args.add_covid_wave), wave=args.add_covid_wave).checkpoint(
            get_ukb_pheno_mt_path(f'covid_wave{args.add_covid_wave}'), args.overwrite)
        cov_ht = hl.read_table(get_ukb_covariates_ht_path())
        mt = combine_pheno_files_multi_sex({'custom': mt}, cov_ht)

        mt.group_rows_by(_x=True).aggregate(
            n_cases=hl.agg.count_where(mt.both_sexes == 1.0),
            n_controls=hl.agg.count_where(mt.both_sexes == 0.0)
        ).entries().drop(*[x for x in PHENO_COLUMN_FIELDS if x != 'description']).show(100, width=180)

        original_mt = hl.read_matrix_table(get_ukb_pheno_mt_path())
        original_mt = original_mt.checkpoint(get_ukb_pheno_mt_path(f'full_before_{curdate}', sex='full'),
                                             args.overwrite)
        original_mt.cols().export(f'{pheno_directory}/all_pheno_summary_before_{curdate}.txt.bgz')
        mt = mt.annotate_cols(Abbvie_Priority='h', Biogen_Priority='h', Pfizer_Priority='h', score=6)
        mt = original_mt.union_cols(mt, row_join_type='outer')
        mt.write(get_ukb_pheno_mt_path(), args.overwrite)
    # # TODO: pfizer male(s) bug
    # original_mt = hl.read_matrix_table(get_ukb_pheno_mt_path())
    # original_mt = original_mt.checkpoint(get_ukb_pheno_mt_path(f'full_before_{curdate}', sex='full'),
    #                                      True)
    # mt = original_mt
    # mt = mt.key_cols_by('trait_type', 'phenocode', pheno_sex=hl.if_else(mt.pheno_sex.endswith('male'),
    #                                                                     mt.pheno_sex + 's', mt.pheno_sex),
    #                     coding=mt.coding, modifier=mt.modifier)
    # mt.write(get_ukb_pheno_mt_path(), True)

    ht = hl.read_matrix_table(get_ukb_pheno_mt_path()).cols()
    exclusions = set()
    exclusions_tranche = '200k'
    with hl.hadoop_open(f'{bucket}/{exclusions_tranche}/results/exclusions_{exclusions_tranche}.txt') as f:
        for line in f:
            pheno = line.strip().split('\t')
            if len(pheno) == 2:
                pheno, coding = pheno
            else:
                pheno = pheno[0]
            exclusions.add((pheno, coding))

    ht = ht.filter(
        ~hl.set(exclusions).contains((ht.phenocode, ht.coding)) &
        (
                hl.set({'biogen', 'abbvie', 'pfizer'}).contains(ht.coding) |
                (
                        (ht[f'n_cases_both_sexes'] >= MIN_CASES) &
                        ((ht.score >= 1) | (ht.trait_type == 'icd'))
                )
        )
    )
    ht.export(f'{pheno_directory}/all_pheno_summary_to_run.txt.bgz')

if __name__ == '__main__':
    parser = argparse.ArgumentParser()

    parser.add_argument('--overwrite', help='Overwrite everything', action='store_true')
    parser.add_argument('--load_data', help='Load data', action='store_true')
    parser.add_argument('--combine_data', help='Combine datasets', action='store_true')
    parser.add_argument('--export_data', help='Export data for SAIGE', action='store_true')
    parser.add_argument('--add_covid_wave', help='Load data')
    parser.add_argument('--slack_channel', help='Send message to Slack channel/user', default='@konradjk')
    args = parser.parse_args()

    if args.slack_channel:
        from slack_token_pkg.slack_creds import slack_token
        with slack.slack_notifications(slack_token, args.slack_channel):
            main(args)
