#!/usr/bin/env python3

__author__ = 'konradk'

from gnomad_hail import *
from ukb_exomes import *
from ukb_phenotypes import *

temp_bucket = 'gs://ukbb-pharma-exome-analysis-temp'


def run_phenotype_missingness_pca(data_type, input_mt_path: str, overwrite: bool = False):
    # TODO: Move to notebook
    scores_ht_path = f'{temp_bucket}/{data_type}_pca.ht'

    mt = hl.read_matrix_table(input_mt_path)
    if data_type == 'icd':
        mt = mt.annotate_rows(mean_primary_codes=hl.agg.mean(mt.primary_codes))
        mt = mt.annotate_entries(pheno_norm=hl.int(mt.primary_codes) - mt.mean_primary_codes)
    else:
        mt = mt.filter_cols(hl.is_defined(mt.both_sexes_pheno))
        mt = mt.annotate_entries(pheno_defined=hl.is_defined(mt.both_sexes))
        mt = mt.annotate_rows(mean_pheno_defined=hl.agg.mean(mt.pheno_defined))
        mt = mt.annotate_entries(pheno_norm=mt.pheno_defined - mt.mean_pheno_defined)
    eigenvalues, scores, _ = hl.pca(mt.pheno_norm, compute_loadings=False)
    print('Eigenvalues: {}'.format(eigenvalues))
    ht = scores.checkpoint(scores_ht_path, overwrite)

    mt = hl.read_matrix_table(input_mt_path)
    if data_type == 'icd':
        mt = mt.annotate_cols(fraction_called=hl.agg.fraction(mt.primary_codes))
        ht = ht.annotate(fraction_called=mt.cols()[ht.key].fraction_called)
    else:
        mt = mt.annotate_cols(fraction_called=hl.agg.fraction(hl.is_defined(mt.both_sexes)))
        ht = ht.annotate(**mt.cols()[ht.key].both_sexes_pheno, fraction_called=mt.cols()[ht.key].fraction_called)
    ht = ht.transmute(**{'PC' + str(i + 1): ht.scores[i] for i in range(10)})
    ht.export(scores_ht_path.replace('.ht', '.txt.bgz'))


def main(args):
    sexes = ('both_sexes_no_sex_specific', 'females', 'males')
    data_types = ('categorical', 'continuous')
    if args.load_data:
        pre_process_data_dictionary(pheno_description_raw_path, pheno_description_path)
        get_codings().write(coding_ht_path, overwrite=args.overwrite)
        load_icd_data(pre_phesant_data_path, icd_codings_path, temp_bucket).write(get_ukb_pheno_mt_path('icd', 'full'), args.overwrite)
        read_covariate_data(pre_phesant_data_path).write(get_ukb_covariates_ht_path(), args.overwrite)

        for sex in sexes:
            pheno_ht = hl.import_table(get_ukb_source_tsv_path(sex), impute=True, min_partitions=100, missing='', key='userID')
            pheno_ht.write(get_ukb_pheno_ht_path(sex), overwrite=args.overwrite)

            pheno_ht = hl.read_table(get_ukb_pheno_ht_path(sex))
            for data_type in data_types:
                pheno_ht_to_mt(pheno_ht, data_type).write(get_ukb_pheno_mt_path(data_type, sex), args.overwrite)

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

    # run_phenotype_missingness_pca('categorical', get_ukb_pheno_mt_path('categorical', 'full'))
    # See gs://ukbb-pharma-exome-analysis/missingness_correlation.ipynb


if __name__ == '__main__':
    parser = argparse.ArgumentParser()

    parser.add_argument('--overwrite', help='Overwrite everything', action='store_true')
    parser.add_argument('--load_data', help='Load data', action='store_true')
    parser.add_argument('--combine_data', help='Combine datasets', action='store_true')
    parser.add_argument('--slack_channel', help='Send message to Slack channel/user', default='@konradjk')
    args = parser.parse_args()

    if args.slack_channel:
        try_slack(args.slack_channel, main, args)
    else:
        main(args)
