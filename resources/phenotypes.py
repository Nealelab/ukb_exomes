import hail as hl

CURRENT_TRANCHE = '200k'
TRANCHE_DATA = {
    '100k': 'January_2019',
    '200k': 'November_2019'
}

pheno_description_with_priority_path = 'gs://phenotype_pharma/misc_files/Data_Dictionary_Showcase_with_priority.tsv'
sample_mapping_file = 'gs://phenotype_pharma/PHESANT_input/linking_files/Project_26041_bridge.csv'


def get_ukb_phesant_summary_tsv_path(sex: str = 'both_sexes_no_sex_specific'):
    return f'gs://phenotype_pharma/misc_files/phesant_output_combined_{sex}_summary.modified.tsv'


def get_ukb_source_tsv_path(sex: str = 'both_sexes_no_sex_specific', tranche: str = CURRENT_TRANCHE):
    return f'gs://phenotype_pharma/PHESANT_output/{tranche}/{TRANCHE_DATA[tranche]}/phesant_output_combined_{sex}.tsv'


def get_biomarker_source_tsv_path(sex: str = 'both_sexes_no_sex_specific', tranche: str = CURRENT_TRANCHE, subdir: str = 'April_2019'):
    return f'gs://phenotype_pharma/PHESANT_output/{tranche}/{subdir}/biomarkers/phesant_output_combined_{sex}.tsv'


def get_biomarker_ht_path(sex: str = 'both_sexes_no_sex_specific', tranche: str = CURRENT_TRANCHE, subdir: str = 'April_2019'):
    return f'gs://phenotype_pharma/PHESANT_output/{tranche}/{subdir}/phenotype/biomarkers_phesant_output_combined_{sex}.ht'


def get_ukb_pheno_ht_path(sex: str = 'both_sexes_no_sex_specific', tranche: str = CURRENT_TRANCHE):
    return f'gs://ukbb-pharma-exome-analysis/{tranche}/{TRANCHE_DATA[tranche]}/phenotype/ht/{sex}.ht'


def get_ukb_pheno_mt_path(data_type: str, sex: str = 'full', tranche: str = CURRENT_TRANCHE):
    return f'gs://ukbb-pharma-exome-analysis/{tranche}/{TRANCHE_DATA[tranche]}/phenotype/{sex}/{data_type}.mt'


def get_pre_phesant_data_path(tranche: str = CURRENT_TRANCHE):
    return f'gs://phenotype_pharma/PHESANT_input/pharma_parsed_and_restricted_to_{tranche.upper()}_sample_subset.tsv'


def get_ukb_covariates_ht_path(tranche: str = CURRENT_TRANCHE):
    return f'gs://ukbb-pharma-exome-analysis/misc_files/{tranche}/covariates.ht'


def get_phenotype_summary_tsv_path(data_type: str, tranche: str = CURRENT_TRANCHE):
    return f'gs://ukbb-pharma-exome-analysis/misc/{tranche}/phenos_{data_type}.tsv'
