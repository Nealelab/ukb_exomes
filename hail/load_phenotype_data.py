#!/usr/bin/env python3

__author__ = 'konradk'

from ukb_exomes import *
from ukb_common import *

temp_bucket = 'gs://ukb-pharma-exome-analysis-temp'


def main(args):
    sexes = ('both_sexes_no_sex_specific', 'females', 'males')
    data_types = ('categorical', 'continuous')
    if args.load_data:
        pre_process_data_dictionary(pheno_description_raw_path, pheno_description_path)
        get_codings().write(coding_ht_path, overwrite=args.overwrite)
        load_icd_data(pre_phesant_data_path, icd_codings_path, temp_bucket).write(get_ukb_pheno_mt_path('icd', 'full'), args.overwrite)
        get_full_icd_data_description(icd_codings_path).write(icd_full_codings_ht_path, args.overwrite)
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
            # TODO: figure out what happened to 26413 (Health score (England)), for instance
            priority_ht = hl.import_table(pheno_description_with_priority_path, impute=True, missing='')
            priority_ht = priority_ht.filter(priority_ht.FieldID != '?')
            priority_ht = priority_ht.key_by(FieldID=hl.int32(priority_ht.FieldID)).select('Abbvie_Priority', 'Biogen_Priority', 'Pfizer_Priority')
            mt = mt.annotate_cols(**priority_ht[mt.pheno])
            mt.write(get_ukb_pheno_mt_path(data_type, 'full'), args.overwrite)


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
