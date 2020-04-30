from .generic import *
from ukb_common.resources.results import *


def get_finngen_results_mt_path(release: str = 'r4', by_gene: bool = False):
    release = f'{release}/' if release != 'r2' else ''
    return f'{bucket}/finngen/{release}variant_results{"_by_gene" if by_gene else ""}.mt'


def get_finngen_results_vep_ht_path(release: str = 'r4'):
    release = f'{release}/' if release != 'r2' else ''
    return f'{bucket}/finngen/{release}variant_results_vep.ht'


def get_enriched_hits_mt_path(release: str = 'r4'):
    release = f'{release}/' if release != 'r2' else ''
    return f'{bucket}/finngen/{release}enriched_hits_by_gene.mt'


def get_ukb_finngen_joined_mt_path(release: str = 'r4', fin_by_gene: bool = False):
    release = f'{release}/' if release != 'r2' else ''
    return f'{bucket}/finngen/{release}ukb_finngen_joined{"_fin_by_gene" if fin_by_gene else ""}.mt'


def get_results_dir(tranche: str = CURRENT_TRANCHE, location: str = 'result'):
    return f'{bucket}/{tranche}/results/{location}/'


def get_final_pheno_info_ht_path(tranche: str = CURRENT_TRANCHE):
    return f'{bucket}/{tranche}/results/pheno_info.ht'


def get_results_mt_path(result_type: str = 'gene', tranche: str = CURRENT_TRANCHE, extension: str = 'mt'):
    result_type = '' if result_type == 'gene' else 'variant_'
    return f'{bucket}/{tranche}/results/{result_type}results.{extension}'

