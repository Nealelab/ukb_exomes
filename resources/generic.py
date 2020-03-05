

CURRENT_TRANCHE = '200k'
TRANCHE_PHENO_DATA = {
    '100k': 'January_2019',
    '200k': 'November_2019'
}

TRANCHE_DATA = {
    '100k': ('regeneron', 4),
    '200k': ('broad', 5)
}

public_bucket = 'gs://ukbb-exome-public/'
bucket = 'gs://ukbb-pharma-exome-analysis'
temp_bucket = 'gs://ukbb-pharma-exome-analysis-temp'


TRANCHE_POPS = {
    '200k': {'nfe', '13'}
}


MIN_CASES = 200