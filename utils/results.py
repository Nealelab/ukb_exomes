import hail as hl
from ukbb_qc.resources.variant_qc import var_annotations_ht_path
from gnomad.utils.vep import process_consequences


def compute_lambda_gc_ht(result_type: str = 'gene', tranche: str = CURRENT_TRANCHE, cutoff: float = 0.0001):

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
		output = output.select_cols('n_cases', lambda_gc_skato=lambda_gc)

	else:  # Variant Info
		sub_af = var_af.filter(var_af.AF > cutoff)
		output = result_mt.filter_rows(hl.is_defined(sub_af.index(result_mt['locus'], result_mt['alleles']))) 
		lambda_gc = hl.agg.filter(hl.is_defined(output.Pvalue), hl.methods.statgen._lambda_gc_agg(output.Pvalue))
		output = output.select_cols('n_cases', lambda_gc=lambda_gc)



	output = output.cols()


	return output
