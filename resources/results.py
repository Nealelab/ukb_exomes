from .generic import *

finngen_results_mt_path = f'{bucket}/finngen/variant_results.mt'
finngen_results_vep_ht_path = f'{bucket}/finngen/variant_results_vep.ht'


def get_results_dir(tranche: str = CURRENT_TRANCHE, location: str = 'result'):
    return f'{bucket}/{tranche}/results/{location}/'


def get_final_pheno_info_ht_path(tranche: str = CURRENT_TRANCHE):
    return f'{bucket}/{tranche}/results/pheno_info.ht'


def get_results_mt_path(result_type: str = 'gene', tranche: str = CURRENT_TRANCHE, extension: str = 'mt'):
    result_type = '' if result_type == 'gene' else 'variant_'
    return f'{bucket}/{tranche}/results/{result_type}results.{extension}'

