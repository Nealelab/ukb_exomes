

CURRENT_TRANCHE = '300k'
TRANCHE_PHENO_DATA = {
    '100k': 'January_2019',
    '200k': 'November_2019',
    '300k': 'April_2020',
    '500k': 'April_2020'
}

TRANCHE_DATA = {
    '100k': ('regeneron', 4),
    '200k': ('broad', 5),
    '300k': ('broad', 6)
}

BUCKETS = {'300k': 'gs://ukbb-pharma-exome-analysis-300k',
           '500k': 'gs://ukbb-pharma-exome-analysis',
           '100k': 'gs://ukbb-pharma-exome-analysis-backup',
           '200k': 'gs://ukbb-pharma-exome-analysis-backup'}

public_bucket = 'gs://ukbb-exome-public'
bucket = BUCKETS[CURRENT_TRANCHE]
temp_bucket = 'gs://ukbb-pharma-exome-analysis-temp'


TRANCHE_POPS = {
    '200k': {'nfe', '13'},
}


MIN_CASES = 200