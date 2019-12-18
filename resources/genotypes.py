import hail as hl
import ukbb_qc.resources as ukb
from .generic import *

bucket = 'gs://ukbb-pharma-exome-analysis'
root = f'{bucket}/mt'
ukb_for_grm_mt_path = f'{root}/ukb.for_grm.mt'
ukb_for_grm_pruned_ht_path = f'{root}/ukb.for_grm.pruned.ht'
ukb_for_grm_plink_path = f'{root}/ukb.for_grm.pruned.plink'


def get_ukb_exomes_mt_path(tranche: str = CURRENT_TRANCHE):
    return ukb.get_ukbb_data_path(*TRANCHE_DATA[tranche], hardcalls=True, split=True)


def get_ukb_exomes_mt(tranche: str = CURRENT_TRANCHE, adj: bool = True):
    return ukb.get_ukbb_data(*TRANCHE_DATA[tranche], split=True, adj=adj, meta_root='meta')


def get_ukb_exomes_meta_ht_path(tranche: str = CURRENT_TRANCHE):
    return ukb.meta_ht_path(*TRANCHE_DATA[tranche])


def get_ukb_exomes_qual_ht_path(tranche: str = CURRENT_TRANCHE):
    return ukb.var_annotations_ht_path(*TRANCHE_DATA[tranche], 'vqsr')


def get_processed_ukb_exomes_mt(adj=False):
    mt = get_ukb_exomes_mt(adj=adj)
    mt = mt.key_cols_by(exome_id=mt.array_id)
    qual_ht = hl.read_table(get_ukb_exomes_qual_ht_path())
    mt = mt.annotate_rows(filters=qual_ht[mt.row_key].filters)
    return mt.filter_cols(hl.is_defined(mt.meta) & hl.is_defined(mt.exome_id))


def get_filtered_mt(adj=False):
    mt = get_processed_ukb_exomes_mt(adj=adj)
    mt = mt.filter_rows(~mt.meta.is_filtered & (mt.meta.hybrid_pop == '12'))
    return mt.filter_rows(hl.len(mt.filters) == 0)


def get_ukb_vep_path():
    return f'{root}/ukb.exomes.vep.ht'


def get_ukb_samples_file_path():
    return f'{root}/ukb.exomes.samples'


def get_ukb_gene_map_ht_path(post_processed=True):
    return f'{root}/ukb.exomes.gene_map{"" if post_processed else ".raw"}.ht'


genotypes_vep_path = f'{root}/ukb.genotypes.vep.ht'
