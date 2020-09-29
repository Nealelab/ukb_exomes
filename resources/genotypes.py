import ukbb_qc.resources as ukb
import ukbb_qc.utils.utils as ukb_utils
from gnomad.utils.sparse_mt import *
from .generic import *


def get_ukb_grm_mt_path(tranche: str = CURRENT_TRANCHE, coalesced: bool = False):
    return f'{bucket}/{tranche}/misc/ukb.for_grm{".coalesced" if coalesced else ""}.mt'


def get_ukb_sites_ht_path(tranche: str = CURRENT_TRANCHE):
    return f'{bucket}/{tranche}/misc/ukb.sites_for_grm.ht'


def get_ukb_grm_pruned_ht_path(tranche: str = CURRENT_TRANCHE):
    return f'{bucket}/{tranche}/misc/ukb.for_grm.pruned.ht'


def get_ukb_grm_plink_path(tranche: str = CURRENT_TRANCHE):
    return f'{bucket}/{tranche}/misc/ukb.for_grm.pruned.plink'


def get_ukb_samples_file_path(tranche: str = CURRENT_TRANCHE):
    return f'{bucket}/{tranche}/misc/ukb.exomes.samples'


def get_ukb_gene_map_ht_path(post_processed=True, tranche: str = CURRENT_TRANCHE):
    return f'{bucket}/{tranche}/misc/ukb.exomes.gene_map{"" if post_processed else ".raw"}.ht'


def get_ukb_gene_summary_ht_path(tranche: str = CURRENT_TRANCHE):
    return f'{bucket}/{tranche}/misc/ukb.exomes.gene_summary.ht'


def get_ukb_exomes_mt_path(tranche: str = CURRENT_TRANCHE):
    return ukb.get_ukbb_data_path(*TRANCHE_DATA[tranche], hardcalls=True, split=True)


def get_ukb_exomes_mt(tranche: str = CURRENT_TRANCHE, adj: bool = True):
    return ukb.get_ukbb_data(*TRANCHE_DATA[tranche], split=True, adj=adj, meta_root='meta', key_by_locus_and_alleles=True)


def get_ukb_exomes_meta_ht_path(tranche: str = CURRENT_TRANCHE):
    return ukb.meta_ht_path(*TRANCHE_DATA[tranche])


def get_ukb_exomes_qual_ht_path(tranche: str = CURRENT_TRANCHE):
    qual_source = 'vqsr' if tranche == '100k' else 'rf'
    return ukb.var_annotations_ht_path(qual_source, *TRANCHE_DATA[tranche])


def get_processed_ukb_exomes_mt(adj=False, key='ukbb_app_26041_id', tranche: str = CURRENT_TRANCHE,
                                densify_sites_ht: hl.Table = None, interval: str = None, semi_join_rows: bool = False):
    mt = get_ukb_exomes_mt(tranche=tranche, adj=adj)
    mt = mt.key_cols_by(**{key: mt.meta.ukbb_meta[key]})
    qual_ht = hl.read_table(get_ukb_exomes_qual_ht_path(tranche))

    if tranche not in ('100k', '200k'):
        last_ref_block_end_ht = hl.read_table(ukb.last_END_positions_ht_path(TRANCHE_DATA[tranche][1]))
        if densify_sites_ht is None:
            densify_sites_ht = qual_ht
        if interval is not None:
            densify_sites_ht = hl.read_table(ukb.var_annotations_ht_path('ukb_freq', *TRANCHE_DATA[CURRENT_TRANCHE]))
            densify_sites_ht = hl.filter_intervals(densify_sites_ht, [hl.parse_locus_interval(interval)])
        mt = densify_sites(mt, densify_sites_ht, last_ref_block_end_ht, semi_join_rows=semi_join_rows)

    mt = mt.annotate_rows(filters=qual_ht[mt.row_key].filters)
    return mt.filter_cols(hl.is_defined(mt.meta) & hl.is_defined(mt[key]))


def get_filtered_mt(interval: str = None, adj=False, interval_filter=False, tranche: str = CURRENT_TRANCHE,
                    densify_sites_ht: hl.Table = None, semi_join_rows: bool = False):
    mt = get_processed_ukb_exomes_mt(adj=adj, tranche=tranche, densify_sites_ht=densify_sites_ht, interval=interval, semi_join_rows=semi_join_rows)
    if tranche == '100k':
        mt = mt.filter_cols(~mt.meta.is_filtered & (mt.meta.hybrid_pop == '12'))
    elif tranche == '200k':
        mt = mt.filter_cols(mt.meta.high_quality & ~mt.meta.duplicate &
                            hl.set(TRANCHE_POPS[tranche]).contains(mt.meta.hybrid_pop))
    else:
        mt = mt.filter_cols(mt.meta.sample_filters.high_quality & ~mt.meta.sample_filters.duplicate &
                            (mt.meta.gnomad_pc_project_pop_data.pop == 'nfe'))
    if interval_filter:
        mt = ukb_utils.annotate_interval_qc_filter(*TRANCHE_DATA[tranche], t=mt)
        mt = mt.filter_rows(mt.interval_qc_pass)
    return mt.filter_rows(hl.len(mt.filters) == 0)


def get_ukb_vep_path():
    return ukb.var_annotations_ht_path('vep', *TRANCHE_DATA[CURRENT_TRANCHE])
