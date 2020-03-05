#!/usr/bin/env python3

__author__ = 'konradk'

from ukb_common import *
from ukb_exomes import *

temp_bucket = 'gs://ukb-pharma-exome-analysis-temp'
pheno_directory = f'{bucket}/{CURRENT_TRANCHE}/pheno_combo_explore'
pairwise_correlation_ht_path = f'{pheno_directory}/pairwise_correlations.ht'


def read_covariate_data(pre_phesant_data_path):
    ht = hl.import_table(pre_phesant_data_path, impute=True, min_partitions=100, missing='', key='userId')
    columns = {
        'sex': 'x22001_0_0',
        'age': 'x21022_0_0'
    }
    columns.update(**{f'pc{i}': f'x22009_0_{i}' for i in range(1, 41)})
    ht = ht.select(*columns.values()).rename({v: k for k, v in columns.items()}).annotate_globals(coding_source=columns)
    meta_ht = hl.read_table(get_ukb_exomes_meta_ht_path(CURRENT_TRANCHE))
    batches = meta_ht.aggregate(hl.agg.counter(meta_ht.batch))
    meta_ht = meta_ht.key_by(userId=hl.int32(meta_ht.ukbb_app_26041_id))
    print(f'Found batches: {batches}')
    batch_field = meta_ht[ht.key].batch
    # One-hot encoding each batch (except 0)
    ht = ht.annotate(**{f'batch_{batch}': hl.int(batch_field == batch) for batch in batches if batch != '0' and batch is not None})
    return ht.annotate(age2=ht.age ** 2, age_sex=ht.age * ht.sex, age2_sex=ht.age ** 2 * ht.sex)


def extract_mt_by_type(ht: hl.Table, value_type = hl.expr.Int32Expression, data_type = 'categorical', source = 'biogen'):
    columns = [x[0] for x in list(ht.row_value.items()) if isinstance(x[1], value_type)]
    ht2 = ht.select(values=[hl.struct(value=hl.float64(ht[x])) for x in columns]).annotate_globals(
        columns=hl.map(lambda x: hl.struct(pheno=x, coding=source, data_type=data_type), hl.literal(columns)))
    return ht2._unlocalize_entries('values', 'columns', ['pheno', 'coding'])


def load_custom_data():
    ht = hl.import_table('gs://phenotype_pharma/custom_phenos/new_phenos_BIOGEN.csv',
                         impute=True, delimiter=',', quote='"', missing="", key='UKBB_ID').rename({"UKBB_ID": 'userId'})
    ht = ht.annotate(Fluid_intelligence_score_BI=hl.float64(ht.Fluid_intelligence_score_BI),
                     Depressive_symptoms_BI=hl.float64(ht.Depressive_symptoms_BI),
                     Touchscreen_duration_BI=hl.float64(ht.Touchscreen_duration_BI))

    mt = extract_mt_by_type(ht, hl.expr.Int32Expression, data_type='categorical')
    mt2 = extract_mt_by_type(ht, hl.expr.Float64Expression, data_type='continuous')
    mt = mt.union_cols(mt2)

    ht = hl.import_table('gs://phenotype_pharma/custom_phenos/Custom_Phenotypes_AbbVie_02042020_Merged.csv', impute=True, delimiter=' ', missing="NA", key='userId')
    mt_abbvie = extract_mt_by_type(ht, hl.expr.Int32Expression, data_type='categorical', source='abbvie')
    return mt.union_cols(mt_abbvie)


def main(args):
    hl.init(log='/load_pheno.log')
    sexes = ('both_sexes_no_sex_specific', 'females', 'males')
    data_types = ('categorical', 'continuous') #, 'biomarkers')

    if args.load_data:
        load_icd_data(get_pre_phesant_data_path(), icd_codings_ht_path, f'{temp_bucket}/{CURRENT_TRANCHE}').write(get_ukb_pheno_mt_path('icd', 'full'), args.overwrite)
        read_covariate_data(get_pre_phesant_data_path()).write(get_ukb_covariates_ht_path(), args.overwrite)
        # load_prescription_data(prescription_tsv_path, prescription_mapping_path).write(get_ukb_pheno_mt_path('prescriptions'), args.overwrite)
        load_custom_data().write(get_ukb_pheno_mt_path('custom', 'full'), args.overwrite)

        for sex in sexes:
            path = get_ukb_source_tsv_path(sex)
            pheno_ht = hl.import_table(path, impute=True, min_partitions=100, missing='', key='userId', force_bgz=path.endswith('.gz'))
            pheno_ht.write(get_ukb_pheno_ht_path(sex), overwrite=args.overwrite)

            pheno_ht = hl.read_table(get_ukb_pheno_ht_path(sex))
            for data_type in data_types:
                pheno_ht_to_mt(pheno_ht, data_type).write(get_ukb_pheno_mt_path(data_type, sex), args.overwrite)

            if CURRENT_TRANCHE == '100k':
                ht = hl.import_table(get_biomarker_source_tsv_path(sex), impute=True, min_partitions=100, missing='', key='userID')
                ht = ht.checkpoint(get_biomarker_ht_path(sex), overwrite=args.overwrite)

                pheno_ht_to_mt(ht, 'continuous').write(get_ukb_pheno_mt_path('biomarkers', sex), overwrite=args.overwrite)

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
            mt = combine_datasets({sex: get_ukb_pheno_mt_path(data_type, sex) for sex in sexes},
                                  {sex: get_ukb_phesant_summary_tsv_path(sex) for sex in sexes},
                                  pheno_description_path, coding_ht_path, data_type)
            mt.write(get_ukb_pheno_mt_path(data_type, 'full'), args.overwrite)
        data_types = ('categorical', 'continuous', 'icd', 'custom')  # 'prescriptions',
        pheno_file_dict = {data_type: hl.read_matrix_table(get_ukb_pheno_mt_path(data_type, 'full')) for data_type in data_types}
        cov_ht = hl.read_table(get_ukb_covariates_ht_path())
        mt = combine_pheno_files_multi_sex(pheno_file_dict, cov_ht)
        priority_ht = hl.import_table(pheno_description_with_priority_path, impute=True, missing='', key='FieldID'
                                      ).select('Abbvie_Priority', 'Biogen_Priority', 'Pfizer_Priority')
        mt = mt.annotate_cols(**priority_ht[mt.pheno])
        mt = mt.annotate_cols(
            score=2 * (hl.int(mt.Abbvie_Priority == 'h') + hl.int(mt.Biogen_Priority == 'h') + hl.int(mt.Pfizer_Priority == 'h')) +
                  hl.int(mt.Abbvie_Priority == 'm') + hl.int(mt.Biogen_Priority == 'm') + hl.int(mt.Pfizer_Priority == 'm'))
        mt.write(get_ukb_pheno_mt_path(), args.overwrite)

        mt = hl.read_matrix_table(get_ukb_pheno_mt_path())
        mt.cols().export(f'{pheno_directory}/all_pheno_summary.txt.bgz')
        # make_correlation_ht(mt).write(pairwise_correlation_ht_path, args.overwrite)
        # hl.read_table(pairwise_correlation_ht_path).flatten().export(pairwise_correlation_ht_path.replace('.ht', '.txt.bgz'))

    ht = hl.read_matrix_table(get_ukb_pheno_mt_path()).cols()
    exclusions = set()
    with hl.hadoop_open(f'{bucket}/{CURRENT_TRANCHE}/results/exclusions_200k.txt') as f:
        for line in f:
            pheno = line.strip().split('\t')
            if len(pheno) == 2:
                pheno, coding = pheno
            else:
                pheno = pheno[0]
            exclusions.add((pheno, coding))

    ht = ht.filter(
        ~hl.set(exclusions).contains((ht.pheno, ht.coding)) &
        (
                hl.set({'biogen', 'abbvie', 'pfizer'}).contains(ht.coding) |
                (
                        (ht[f'n_cases_both_sexes'] >= MIN_CASES) &
                        ((ht.score >= 1) | (ht.data_type == 'icd'))
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
    parser.add_argument('--slack_channel', help='Send message to Slack channel/user', default='@konradjk')
    args = parser.parse_args()

    if args.slack_channel:
        try_slack(args.slack_channel, main, args)
    else:
        main(args)
