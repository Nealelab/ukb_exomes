from .generic import *
from ukb_common.resources.results import *


CURRENT_FG_RELEASE = 'r5'


def get_finngen_results_mt_path(release: str = CURRENT_FG_RELEASE, by_gene: bool = False):
    release = f'{release}/' if release != 'r2' else ''
    return f'{bucket}/finngen/{release}variant_results{"_by_gene" if by_gene else ""}.mt'


def get_finngen_results_vep_ht_path(release: str = CURRENT_FG_RELEASE):
    release = f'{release}/' if release != 'r2' else ''
    return f'{bucket}/finngen/{release}variant_results_vep.ht'


def get_enriched_hits_mt_path(release: str = CURRENT_FG_RELEASE):
    release = f'{release}/' if release != 'r2' else ''
    return f'{bucket}/finngen/{release}enriched_hits_by_gene.mt'


def get_ukb_finngen_joined_mt_path(release: str = CURRENT_FG_RELEASE, fin_by_gene: bool = False):
    release = f'{release}/' if release != 'r2' else ''
    return f'{bucket}/finngen/{release}ukb_finngen_joined{"_fin_by_gene" if fin_by_gene else ""}.mt'


def get_results_dir(tranche: str = CURRENT_TRANCHE, location: str = 'result'):
    return f'{bucket}/{tranche}/results/{location}/'


def get_final_pheno_info_ht_path(tranche: str = CURRENT_TRANCHE):
    return f'{bucket}/{tranche}/results/pheno_info.ht'


def get_results_mt_path(result_type: str = 'gene', tranche: str = CURRENT_TRANCHE, extension: str = 'mt',
                        biomarkers: bool = False, random_phenos: bool = False):
    if random_phenos:
        result_type = f'random_phenos_{result_type}_'
    elif biomarkers:
        result_type = f'biomarkers_{result_type}'
    else:
        result_type = '' if result_type == 'gene' else 'variant_'
    return f'{bucket}/{tranche}/results/{result_type}results.{extension}'


def get_loeuf_results_path(analysis_type: str, extension: str = 'mt', tranche: str = CURRENT_TRANCHE):
    return f'{bucket}/{tranche}/loeuf_analysis/{analysis_type}.{extension}'


def get_results_timing_tsv_path(timing_type: str, trait_type: str = 'all', tranche: str = CURRENT_TRANCHE):
    check_timing_type(timing_type)

    if timing_type == 'saige':
        check_trait_types(trait_type)
        trait_type = f'_{trait_type}'
    else:
        trait_type = ''

    return f'{bucket}/{tranche}/results/misc/timings_{timing_type}{trait_type}.txt'


def get_results_timing_ht_path(timing_type: str, tranche: str = CURRENT_TRANCHE):
    check_timing_type(timing_type)
    return f'{bucket}/{tranche}/results/misc/timings_{timing_type}.ht'

