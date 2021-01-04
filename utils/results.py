from ukbb_qc.resources.variant_qc import var_annotations_ht_path
from ukbb_qc.resources.basics import release_ht_path
from ukb_common.utils.annotations import annotation_case_builder
from gnomad.utils.vep import process_consequences
from ukb_exomes.resources import *


def compute_lambda_gc_ht(result_type: str = 'gene', tranche: str = CURRENT_TRANCHE, cutoff: float = 0.0001):
    af_info = get_af_info_ht(result_type, tranche)
    result_mt = hl.read_matrix_table(get_results_mt_path(result_type, tranche))

    if result_type == 'gene':  # Gene Info
        sub_af = af_info.filter(af_info.caf > cutoff)
        output = result_mt.filter_rows(hl.is_defined(sub_af.index(result_mt['annotation'], result_mt['gene_id'])))
        lambda_gc = hl.agg.filter(hl.is_defined(output.Pvalue), hl.methods.statgen._lambda_gc_agg(output.Pvalue))
        output = output.select_cols('n_cases', lambda_gc_skato=lambda_gc)
    else:  # Variant Info
        sub_af = af_info.filter(af_info.AF > cutoff)
        output = result_mt.filter_rows(hl.is_defined(sub_af.index(result_mt['locus'], result_mt['alleles'])))
        lambda_gc = hl.agg.filter(hl.is_defined(output.Pvalue), hl.methods.statgen._lambda_gc_agg(output.Pvalue))
        output = output.select_cols('n_cases', lambda_gc=lambda_gc)

    output = output.cols()
    return output


def compute_ukb_pheno_moments_ht(pheno_sex='both_sexes', phenocode: list = None):
    pheno = get_ukb_pheno_mt()
    if phenocode is not None:
        pheno = pheno.filter_cols(hl.literal(phenocode).contains(pheno.phenocode))
    pheno = pheno.annotate_cols(mean=hl.agg.mean(pheno.both_sexes))
    pheno = pheno.annotate_cols(variance=hl.agg.mean((pheno.both_sexes - pheno.mean) ** 2))
    pheno = pheno.annotate_cols(skewness=hl.agg.mean((pheno.both_sexes - pheno.mean) ** 3) / (pheno.variance ** 1.5),
                                kurtosis=hl.agg.mean((pheno.both_sexes - pheno.mean) ** 4) / (pheno.variance ** 2))
    output = pheno.select_cols('mean', 'variance', 'skewness', 'kurtosis').cols()
    return output


def get_af_info_ht(result_type: str = 'gene', tranche: str = CURRENT_TRANCHE):
    vep = hl.read_table(var_annotations_ht_path('vep', *TRANCHE_DATA[tranche]))
    vep = process_consequences(vep)
    vep = vep.explode(vep.vep.worst_csq_by_gene_canonical)
    annotation = annotation_case_builder(vep.vep.worst_csq_by_gene_canonical)
    vep = vep.annotate(
        annotation=hl.if_else(hl.literal({'missense', 'LC'}).contains(annotation), 'missense|LC', annotation))
    vep = vep.filter((hl.is_defined(vep.annotation)) & (vep.annotation != 'non-coding'))
    vep = vep.select(vep.vep.worst_csq_by_gene_canonical.gene_id, vep.annotation)

    freq = hl.read_table(release_ht_path(*TRANCHE_DATA[tranche]))
    var_af = vep.annotate(AF=freq[vep.key].freq[0].AF)

    if result_type == 'gene':
        sum_af = var_af.group_by(var_af.annotation, var_af.gene_id).aggregate(caf=hl.agg.sum(var_af.AF))
        return sum_af
    else:
        return var_af


def get_sig_cnt_mt(tranche: str = CURRENT_TRANCHE, phenocode: list = None, gene_symbol: list = None,
                   level: float = 1e-6):
    # Load and Combine Data
    af = get_af_info_ht(result_type='gene')
    mt = hl.read_matrix_table(get_results_mt_path())
    sub = mt.annotate_rows(caf=af[mt.gene_id, mt.annotation]['caf'])
    # Count Significant Hits per Gene for each Test
    sub = sub.annotate_rows(sig_pheno_cnt_skato=hl.agg.sum(sub.Pvalue < level),
                            sig_pheno_cnt_skat=hl.agg.sum(sub.Pvalue_SKAT < level),
                            sig_pheno_cnt_burden=hl.agg.sum(sub.Pvalue_Burden < level))

    # Count Significant Hits per Phenotype for each Test
    sub = sub.annotate_cols(sig_gene_cnt_skato=hl.agg.sum(sub.Pvalue < level),
                            sig_gene_cnt_skat=hl.agg.sum(sub.Pvalue_SKAT < level),
                            sig_gene_cnt_burden=hl.agg.sum(sub.Pvalue_Burden < level))
    if phenocode is not None:
        sub = sub.filter_cols(hl.literal(phenocode).contains(sub.phenocode))
    if gene_symbol is not None:
        sub = sub.filter_rows(hl.literal(gene_symbol).contains(sub.gene_symbol))

    output = sub.select_rows('sig_pheno_cnt_skato', 'sig_pheno_cnt_skat', 'sig_pheno_cnt_burden')
    output = output.select_cols('sig_gene_cnt_skato', 'sig_gene_cnt_skat', 'sig_gene_cnt_burden')
    return output
