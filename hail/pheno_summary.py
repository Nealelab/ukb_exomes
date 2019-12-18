from ukb_common import *
from ukb_exomes import *

temp_bucket = 'gs://ukb-pharma-exome-analysis-temp'


def main(args):
    hl.init(default_reference='GRCh38', log='/pheno_agg.log')

    ht = get_filtered_mt().cols()
    print(ht.count())

    mt = hl.read_matrix_table(get_ukb_pheno_mt_path(args.data_type, 'full'))
    mt = mt.filter_rows(hl.is_defined(ht[hl.str(mt.userId)]))

    if args.data_type == 'continuous':
        mt = mt.filter_cols(mt.coding != 'raw')

    extra_fields = compute_n_cases(mt, args.data_type)
    if args.data_type == 'icd':
        ht = mt.select_cols('truncated', 'meaning', **extra_fields).cols()
    else:
        ht = mt.select_cols('Abbvie_Priority', 'Biogen_Priority', 'Pfizer_Priority', **extra_fields).cols()
        ht = ht.annotate(
            score=2 * (hl.int(ht.Abbvie_Priority == 'h') + hl.int(ht.Biogen_Priority == 'h') + hl.int(ht.Pfizer_Priority == 'h')) +
                  hl.int(ht.Abbvie_Priority == 'm') + hl.int(ht.Biogen_Priority == 'm') + hl.int(ht.Pfizer_Priority == 'm'))
    ht.export(get_phenotype_summary_tsv_path(args.data_type))



if __name__ == '__main__':
    parser = argparse.ArgumentParser()

    parser.add_argument('--overwrite', help='Overwrite everything', action='store_true')
    parser.add_argument('--data_type', help='Data type', required=True, choices=('icd', 'continuous', 'categorical', 'biomarkers'))
    parser.add_argument('--slack_channel', help='Send message to Slack channel/user', default='@konradjk')
    args = parser.parse_args()

    if args.slack_channel:
        try_slack(args.slack_channel, main, args)
    else:
        main(args)