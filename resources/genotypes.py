import hail as hl
import ukbb_qc.resources as ukb

bucket = 'gs://ukbb-pharma-exome-analysis'
root = f'{bucket}/mt'
ukb_for_grm_mt_path = f'{root}/ukb.for_grm.mt'
ukb_for_grm_pruned_ht_path = f'{root}/ukb.for_grm.pruned.ht'
ukb_for_grm_plink_path = f'{root}/ukb.for_grm.pruned.plink'


def get_ukb_exomes_mt_path():
    return ukb.get_ukbb_data_path('broad', 4, hardcalls=True, split=True)


def get_ukb_exomes_meta_ht_path():
    return 'gs://broad-ukbb/regeneron.freeze_4/sample_qc/meta.ht'


def get_ukb_exomes_mt():
    from .phenotypes import sample_mapping_file
    mt = hl.read_matrix_table(get_ukb_exomes_mt_path())
    meta = hl.read_table(get_ukb_exomes_meta_ht_path())
    mapping = hl.import_table(sample_mapping_file, key='eid_sample', delimiter=',')
    mt = mt.annotate_cols(meta=meta[mt.col_key],
                          exome_id=mapping[mt.s.split('_')[1]].eid_26041)
    return mt.filter_cols(hl.is_defined(mt.exome_id)).key_cols_by('exome_id')


def get_ukb_vep_path():
    return f'{root}/ukb.exomes.vep.ht'


def get_ukb_samples_file_path():
    return f'{root}/ukb.exomes.samples'


def get_ukb_gene_map_ht_path(post_processed=True):
    return f'{root}/ukb.exomes.gene_map{"" if post_processed else ".raw"}.ht'


genotypes_vep_path = f'{root}/ukb.genotypes.vep.ht'
public_bucket = 'gs://ukbb-exome-public/'
