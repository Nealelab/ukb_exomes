from .generic import *

pheno_description_with_priority_path = 'gs://phenotype_pharma/misc_files/Data_Dictionary_Showcase_with_priority.tsv'
prescription_mapping_path = 'gs://phenotype_pharma/misc_files/ukb_prescription_mapping.tsv'
# prescription_tsv_path = 'gs://ukb31063/ukb31063.gp_scripts.20191008.txt'  # TODO: get gp scripts file


def get_ukb_phesant_summary_tsv_path(sex: str = 'both_sexes_no_sex_specific'):
    return f'gs://phenotype_pharma/misc_files/phesant_output_combined_{sex}_summary.modified.tsv'


def get_ukb_source_tsv_path(sex: str = 'both_sexes_no_sex_specific', tranche: str = CURRENT_TRANCHE):
    ext = '' if CURRENT_TRANCHE == '100k' else '.gz'
    return f'gs://phenotype_pharma/PHESANT_output/{tranche}/{TRANCHE_PHENO_DATA[tranche]}/phesant_output_combined_{sex}.tsv{ext}'


def get_pre_phesant_data_path(tranche: str = CURRENT_TRANCHE):
    return f'gs://phenotype_pharma/PHESANT_input/pharma_parsed_and_restricted_to_{tranche.upper()}_sample_subset.tsv'


def get_biomarker_source_tsv_path(sex: str = 'both_sexes_no_sex_specific'):
    # Note: biomarkers are integrated into continuous variables in future tranches
    return f'gs://phenotype_pharma/PHESANT_output/100k/April_2019/biomarkers/phesant_output_combined_{sex}.tsv'


def get_biomarker_ht_path(sex: str = 'both_sexes_no_sex_specific'):
    # Note: biomarkers are integrated into continuous variables in future tranches
    return f'gs://phenotype_pharma/PHESANT_output/100k/April_2019/phenotype/biomarkers_phesant_output_combined_{sex}.ht'


def get_ukb_pheno_ht_path(sex: str = 'both_sexes_no_sex_specific', tranche: str = CURRENT_TRANCHE):
    return f'{bucket}/{tranche}/{TRANCHE_PHENO_DATA[tranche]}/phenotype/ht/{sex}.ht'


def get_ukb_pheno_mt_path(data_type: str = 'full', sex: str = 'full', tranche: str = CURRENT_TRANCHE):
    return f'{bucket}/{tranche}/{TRANCHE_PHENO_DATA[tranche]}/phenotype/{sex}/{data_type}.mt'


def get_ukb_covariates_ht_path(tranche: str = CURRENT_TRANCHE):
    return f'{bucket}/{tranche}/covariates.ht'


def get_phenotype_summary_tsv_path(data_type: str, tranche: str = CURRENT_TRANCHE):
    return f'{bucket}/{tranche}/pheno_summary/phenos_{data_type}.tsv'
