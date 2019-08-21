import hail as hl

CURRENT_TRANCHE = '100k'
TRANCHE_DATA = {
    '100k': 'January_2019'  # 'April_2019'
}

pheno_description_raw_path = 'gs://phenotype_pharma/misc_files/Data_Dictionary_Showcase.csv'
pheno_description_path = 'gs://phenotype_pharma/misc_files/Data_Dictionary_Showcase.tsv'
coding_ht_path = 'gs://phenotype_pharma/misc_files/categorical_codings.ht'
icd_codings_path = 'gs://phenotype_pharma/misc_files/coding19.tsv'
icd_full_codings_ht_path = 'gs://phenotype_pharma/misc_files/coding19.ht'
pre_phesant_data_path = 'gs://phenotype_pharma/PHESANT_input/pharma_parsed_and_restricted_to_100K_sample_subset.tsv'
sample_mapping_file = 'gs://phenotype_pharma/misc_files/Project_26041_bridge.csv'


def get_ukb_phesant_summary_tsv_path(sex: str = 'both_sexes_no_sex_specific'):
    return f'gs://phenotype_pharma/misc_files/phesant_output_combined_{sex}_summary.modified.tsv'


def get_ukb_source_tsv_path(sex: str = 'both_sexes_no_sex_specific', tranche: str = CURRENT_TRANCHE):
    return f'gs://phenotype_pharma/PHESANT_output/{tranche}/{TRANCHE_DATA[tranche]}/phesant_output_combined_{sex}.tsv'


def get_ukb_pheno_ht_path(sex: str = 'both_sexes_no_sex_specific', tranche: str = CURRENT_TRANCHE):
    return f'gs://ukbb-pharma-exome-analysis/{tranche}/{TRANCHE_DATA[tranche]}/phenotype/ht/{sex}.ht'


def get_ukb_pheno_mt_path(data_type: str, sex: str = 'full', tranche: str = CURRENT_TRANCHE):
    return f'gs://ukbb-pharma-exome-analysis/{tranche}/{TRANCHE_DATA[tranche]}/phenotype/{sex}/{data_type}.mt'


def get_ukb_covariates_ht_path():
    return f'gs://ukbb-pharma-exome-analysis/misc_files/covariates.ht'


