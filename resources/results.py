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
                        biomarkers: bool = False):
    result_type = '' if result_type == 'gene' else 'variant_'
    if biomarkers: result_type = f'biomarkers_{result_type}'
    return f'{bucket}/{tranche}/results/{result_type}results.{extension}'


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


def get_lambda_gc_ht(result_type: str = 'gene', tranche: str = CURRENT_TRANCHE, cutoff: float = 0.0001):

	vep = hl.read_table(var_annotations_ht_path('vep', *TRANCHE_DATA[tranche]))
	vep = process_consequences(vep)
	vep = vep.explode(vep.vep.worst_csq_by_gene_canonical)
	annotation = (hl.case(missing_false=True)
	.when(vep.vep.worst_csq_by_gene_canonical.most_severe_consequence == 'missense_variant', 'missense|LC')
	.when(vep.vep.worst_csq_by_gene_canonical.lof == 'HC', 'pLoF')
	.when(vep.vep.worst_csq_by_gene_canonical.most_severe_consequence == 'synonymous_variant', 'synonymous')
	.or_missing())  # Create ids for different annotations.
	vep = vep.annotate(annotation=annotation)
	vep = vep.filter((hl.is_defined(vep.annotation)) & (vep.vep.worst_csq_by_gene_canonical.biotype == 'protein_coding'))
	vep = vep.select(vep.vep.worst_csq_by_gene_canonical.gene_id, vep.annotation)

	freq = hl.read_table(var_annotations_ht_path('ukb_freq', *TRANCHE_DATA[tranche]))
	var_af = vep.annotate(AF=freq[vep.key].freq[0].AF)
	sum_af = var_af.group_by(var_af.annotation, var_af.gene_id).aggregate(caf=hl.agg.sum(var_af.AF))

	result_mt = hl.read_matrix_table(get_results_mt_path(result_type, tranche))

	if result_type == 'gene':  # Gene Info
		sub_af = sum_af.filter(sum_af.caf > cutoff)
		output = result_mt.filter_rows(hl.is_defined(sub_af.index(result_mt['annotation'], result_mt['gene_id'])))
		lambda_gc = hl.agg.filter(hl.is_defined(output.Pvalue), hl.methods.statgen._lambda_gc_agg(output.Pvalue))
		output = output.annotate_cols(lambda_gc_skato=lambda_gc)
		output = output.select_cols('n_cases', 'lambda_gc_skato')

	else:  # Variant Info
		sub_af = var_af.filter(var_af.AF > cutoff)
		output = result_mt.filter_rows(hl.is_defined(sub_af.index(result_mt['locus'], result_mt['alleles']))) 
		lambda_gc = hl.agg.filter(hl.is_defined(output.Pvalue), hl.methods.statgen._lambda_gc_agg(output.Pvalue))
		output = output.annotate_cols(lambda_gc=lambda_gc)
		output = output.select_cols('n_cases', 'lambda_gc')

	output = output.cols()

	return output
