import hail as hl
import ukbb_qc.resources as ukb
import ukbb_qc.utils as ukb_utils
from .generic import *


def get_ukb_grm_mt_path(tranche: str = CURRENT_TRANCHE, coalesced: bool = False):
    return f'{bucket}/{tranche}/misc/ukb.for_grm{".coalesced" if coalesced else ""}.mt'


def get_ukb_grm_pruned_ht_path(tranche: str = CURRENT_TRANCHE):
    return f'{bucket}/{tranche}/misc/ukb.for_grm.pruned.ht'


def get_ukb_grm_plink_path(tranche: str = CURRENT_TRANCHE):
    return f'{bucket}/{tranche}/misc/ukb.for_grm.pruned.plink'


def get_ukb_samples_file_path(tranche: str = CURRENT_TRANCHE):
    return f'{bucket}/{tranche}/misc/ukb.exomes.samples'


def get_ukb_gene_map_ht_path(post_processed=True, tranche: str = CURRENT_TRANCHE):
    return f'{bucket}/{tranche}/misc/ukb.exomes.gene_map{"" if post_processed else ".raw"}.ht'


def get_ukb_exomes_mt_path(tranche: str = CURRENT_TRANCHE):
    return ukb.get_ukbb_data_path(*TRANCHE_DATA[tranche], hardcalls=True, split=True)


def get_ukb_exomes_mt(tranche: str = CURRENT_TRANCHE, adj: bool = True):
    return ukb.get_ukbb_data(*TRANCHE_DATA[tranche], split=True, adj=adj, meta_root='meta')


def get_ukb_exomes_meta_ht_path(tranche: str = CURRENT_TRANCHE):
    return ukb.meta_ht_path(*TRANCHE_DATA[tranche])


def get_ukb_exomes_qual_ht_path(tranche: str = CURRENT_TRANCHE):
    qual_source = 'vqsr' if tranche == '100k' else 'rf'
    return ukb.var_annotations_ht_path(*TRANCHE_DATA[tranche], qual_source)


def get_processed_ukb_exomes_mt(adj=False, key='ukbb_app_26041_id'):
    mt = get_ukb_exomes_mt(adj=adj)
    mt = mt.key_cols_by(**{key: mt.meta[key]})
    qual_ht = hl.read_table(get_ukb_exomes_qual_ht_path())
    mt = mt.annotate_rows(filters=qual_ht[mt.row_key].filters)
    return mt.filter_cols(hl.is_defined(mt.meta) & hl.is_defined(mt[key]))


def get_filtered_mt(adj=False, interval_filter=False, tranche: str = CURRENT_TRANCHE):
    mt = get_processed_ukb_exomes_mt(adj=adj)
    if tranche == '100k':
        mt = mt.filter_cols(~mt.meta.is_filtered & (mt.meta.hybrid_pop == '12'))
    else:
        mt = mt.filter_cols(mt.meta.high_quality & ~mt.meta.duplicate &
                            hl.set(TRANCHE_POPS[tranche]).contains(mt.meta.hybrid_pop))
    if interval_filter:
        mt = ukb_utils.interval_qc_filter(*TRANCHE_DATA[tranche], t=mt)
    return mt.filter_rows(hl.len(mt.filters) == 0)


def get_ukb_vep_path():
    return ukb.var_annotations_ht_path(*TRANCHE_DATA[CURRENT_TRANCHE], 'vep')
