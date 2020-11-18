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


def get_lambda_gc_ht(result_type: str = 'gene', tranche: str = CURRENT_TRANCHE, 
	vep_ht_path: str = 'gs://broad-ukbb/broad.freeze_6/variant_qc/variant_annotations/vep.ht', 
	freq_ht_path: str = 'gs://broad-ukbb/broad.freeze_6/variant_qc/variant_annotations/ukb_freq.ht', 
	output_path: str = None, extension: str = 'txt.bgz', cutoff: float = 0.0001):

	vep = hl.read_table(vep_ht_path)
	vep = process_consequences(vep)
	vep = vep.explode(vep.vep.worst_csq_by_gene_canonical)
	annt_id = (hl.case(missing_false=True)
	.when(vep.vep.worst_csq_by_gene_canonical.most_severe_consequence == "missense_variant", 1)
	.when(vep.vep.worst_csq_by_gene_canonical.lof == "HC", 2)
	.when(vep.vep.worst_csq_by_gene_canonical.most_severe_consequence == "synonymous_variant", 4)
	.or_missing())  # Create ids for different annotations.
	vep = vep.annotate(annt_id=annt_id)
	vep = vep.filter(hl.is_defined(vep.annt_id))
	vep = vep.filter(vep.vep.worst_csq_by_gene_canonical.biotype == "protein_coding")
	vep = vep.select(vep.vep.worst_csq_by_gene_canonical.gene_id, vep.annt_id)

	freq = hl.read_table(freq_ht_path)
	var_af = vep.key_by('locus', 'alleles').join(freq.key_by('locus', 'alleles'))
	var_af = var_af.select(var_af.gene_id, var_af.annt_id, AF = var_af.freq[0].AF)
	var_af = var_af.filter(hl.is_defined(var_af.AF))
	sum_af = var_af.group_by(var_af.annt_id, var_af.gene_id).aggregate(caf=hl.agg.sum(var_af.AF))


	output_file = f'LambdaGC_AF{cutoff}_{result_type}.{extension}'

	result_type = '' if result_type == 'gene' else 'variant_'
	result_mt = hl.read_matrix_table(get_results_mt_path(result_type, tranche))

	if result_type == '':  # Gene Info
		annt_id = (hl.case(missing_false=True)
			.when(result_mt.annotation == "missense|LC", 1)
			.when(result_mt.annotation == "pLoF", 2)
			.when(result_mt.annotation == "synonymous", 4)
			.or_missing())  # Create ids for different annotations.
		result_mt = result_mt.annotate_rows(annt_id=annt_id)

		sub_af = sum_af.filter(sum_af.caf > cutoff)
		output = result_mt.filter_rows(~hl.is_missing(sub_af.index(result_mt['annt_id'], result_mt['gene_id'])))
		Lambda_GC = hl.agg.filter(hl.is_defined(output.Pvalue), hl.methods.statgen._lambda_gc_agg(output.Pvalue))
		output = output.annotate_cols(Lambda_GC_SKATO=Lambda_GC)
		output = output.key_cols_by('trait_type', 'phenocode')
		output = output.select_cols("n_cases", "Lambda_GC_SKATO")

	else:  # Variant Info
		sub_af = var_af.filter(var_af.AF > cutoff)
		output = result_mt.filter_rows(~hl.is_missing(sub_af.index(result_mt['locus'], result_mt['alleles']))) 
		Lambda_GC = hl.agg.filter(hl.is_defined(output.Pvalue), hl.methods.statgen._lambda_gc_agg(output.Pvalue))
		output = output.annotate_cols(Lambda_GC_Burden=Lambda_GC)
		output = output.key_cols_by('trait_type', 'phenocode')
		output = output.select_cols("n_cases", "Lambda_GC_Burden")



	output = output.cols()

	if output_path is not None:
		output.export(f'{output_path}/{output_file}')

	return output
