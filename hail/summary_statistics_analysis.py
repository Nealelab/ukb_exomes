import hail as hl
from ukb_exomes import *

# --------Preliminary Lambda GC Check----------
# Overall Lambda GC
lambda_gene_full = compute_lambda_gc_ht()
lambda_gene_full = lambda_gene_full.annotate(trait_type2=hl.if_else(lambda_gene_full.trait_type == 'icd_first_occurrence', 'icd10', lambda_gene_full.trait_type))
lambda_gene_full.export('gs://ukbb-exome-public/summary_statistics_analysis/lambda_gene_full300k.txt.bgz')

lambda_var_full = compute_lambda_gc_ht(result_type='variant')
lambda_var_full = lambda_var_full.annotate(trait_type2=hl.if_else(lambda_var_full.trait_type == 'icd_first_occurrence', 'icd10', lambda_var_full.trait_type))
lambda_var_full.export('gs://ukbb-exome-public/summary_statistics_analysis/lambda_var_full300k.txt.bgz')

# Lambda GC with AF filters
g1e5 = compute_lambda_gc_ht(af_upper=0.0001)
g1e4 = compute_lambda_gc_ht(af_lower=0.0001, af_upper=0.001)
g1e3 = compute_lambda_gc_ht(af_lower=0.001, af_upper=0.01)
g1e2 = compute_lambda_gc_ht(af_lower=0.01, af_upper=0.1)
g1e1 = compute_lambda_gc_ht(af_lower=0.1)
lambda_gene_af = g1e5.union(g1e4, g1e3, g1e2, g1e1)
lambda_gene_af = lambda_gene_af.annotate(trait_type2=hl.if_else(lambda_gene_af.trait_type == 'icd_first_occurrence', 'icd10', lambda_gene_af.trait_type),)
lambda_gene_af.export('gs://ukbb-exome-public/summary_statistics_analysis/lambda_gene_af300k.txt.bgz')

v1e5 = compute_lambda_gc_ht(result_type='variant', af_upper=0.0001)
v1e4 = compute_lambda_gc_ht(result_type='variant', af_lower=0.0001, af_upper=0.001)
v1e3 = compute_lambda_gc_ht(result_type='variant', af_lower=0.001, af_upper=0.01)
v1e2 = compute_lambda_gc_ht(result_type='variant', af_lower=0.01, af_upper=0.1)
v1e1 = compute_lambda_gc_ht(result_type='variant', af_lower=0.1)
lambda_var_af = v1e5.union(v1e4, v1e3, v1e2, v1e1)
lambda_var_af = lambda_var_af.annotate(trait_type2=hl.if_else(lambda_var_af.trait_type == 'icd_first_occurrence', 'icd10', lambda_var_af.trait_type),)
lambda_var_af.export('gs://ukbb-exome-public/summary_statistics_analysis/lambda_var_af300k.txt.bgz')

# Lambda GC by Annotation
lambda_gene_annt = compute_lambda_gc_ht(by_annotation=True)
lambda_gene_annt = lambda_gene_annt.annotate(trait_type2=hl.if_else(lambda_gene_annt.trait_type == 'icd_first_occurrence', 'icd10', lambda_gene_annt.trait_type))
lambda_gene_annt.export('gs://ukbb-exome-public/summary_statistics_analysis/lambda_gene_annt300k.txt.bgz')

lambda_var_annt = compute_lambda_gc_ht(result_type='variant', by_annotation=True)
lambda_var_annt = lambda_var_annt.annotate(trait_type2=hl.if_else(lambda_var_annt.trait_type == 'icd_first_occurrence', 'icd10', lambda_var_annt.trait_type))
lambda_var_annt.export('gs://ukbb-exome-public/summary_statistics_analysis/lambda_var_annt300k.txt.bgz')

# Lambda GC with AF filters by Annotation
g1e5 = compute_lambda_gc_ht(by_annotation=True, af_upper=0.0001)
g1e4 = compute_lambda_gc_ht(by_annotation=True, af_lower=0.0001, af_upper=0.001)
g1e3 = compute_lambda_gc_ht(by_annotation=True, af_lower=0.001, af_upper=0.01)
g1e2 = compute_lambda_gc_ht(by_annotation=True, af_lower=0.01, af_upper=0.1)
g1e1 = compute_lambda_gc_ht(by_annotation=True, af_lower=0.1)
lambda_gene_af_annt = g1e5.union(g1e4, g1e3, g1e2, g1e1)
lambda_gene_af_annt = lambda_gene_af_annt.annotate(trait_type2=hl.if_else(lambda_gene_af_annt.trait_type == 'icd_first_occurrence', 'icd10', lambda_gene_af_annt.trait_type),)
lambda_gene_af_annt.export('gs://ukbb-exome-public/summary_statistics_analysis/lambda_gene_af300k.txt.bgz')

v1e5 = compute_lambda_gc_ht(result_type='variant', by_annotation=True, af_upper=0.0001)
v1e4 = compute_lambda_gc_ht(result_type='variant', by_annotation=True, af_lower=0.0001, af_upper=0.001)
v1e3 = compute_lambda_gc_ht(result_type='variant', by_annotation=True, af_lower=0.001, af_upper=0.01)
v1e2 = compute_lambda_gc_ht(result_type='variant', by_annotation=True, af_lower=0.01, af_upper=0.1)
v1e1 = compute_lambda_gc_ht(result_type='variant', by_annotation=True, af_lower=0.1)
lambda_var_af_annt = v1e5.union(v1e4, v1e3, v1e2, v1e1)
lambda_var_af_annt = lambda_var_af_annt.annotate(trait_type2=hl.if_else(lambda_var_af_annt.trait_type == 'icd_first_occurrence', 'icd10', lambda_var_af_annt.trait_type),)
lambda_var_af_annt.export('gs://ukbb-exome-public/summary_statistics_analysis/lambda_var_af300k.txt.bgz')

# Gene-Level Lambda GC
lambda_by_gene = compute_lambda_gc_ht(by_gene = True)
lambda_by_gene.export('gs://ukbb-exome-public/summary_statistics_analysis/lambda_by_gene300k.txt.bgz')

# --------QCed Lambda GC----------
lambda_gene = compute_lambda_gc_ht(result_type='gene', af_lower=1e-4, n_var_min=2)
lambda_by_gene = compute_lambda_gc_ht(result_type='gene', af_lower=1e-4, by_gene=True, n_var_min=2)
lambda_var = compute_lambda_gc_ht(result_type='variant', af_lower=2e-5)

#  --------Genes and Phenotypes to Keep----------
phenos_to_keep = get_lambda_filtered_ht(lambda_gene, lambda_name='lambda_gc_skato',lower=0.75, upper=1.5)
genes_to_keep = get_lambda_filtered_ht(lambda_by_gene, lambda_name='lambda_gc_skato',lower=0.75, upper=1.5)
phenos_to_remove = get_corr_phenos_ht(r_2=0.5, tie_breaker=more_cases_tie_breaker)
phenos_to_keep = phenos_to_keep.filter(hl.is_defined(phenos_to_remove.key_by(trait_type = phenos_to_remove.node.trait_type, phenocode = phenos_to_remove.node.phenocode,
                                                                             pheno_sex = phenos_to_remove.node.pheno_sex, coding = phenos_to_remove.node.coding,
                                                                             modifier = phenos_to_remove.node.modifier, )[phenos_to_keep.key]), keep=False)

# --------Significant Hits----------
# Load and Combine Data
var = get_sig_cnt_mt('variant', phenos_to_keep=phenos_to_keep)
gene = get_sig_cnt_mt(phenos_to_keep=phenos_to_keep, genes_to_keep=genes_to_keep)

# Gene Association Count
gene_sig = gene.rows()
gene_sig = gene_sig.select('caf', data=[(gene_sig.sig_pheno_cnt_skato, 'skato'), (gene_sig.sig_pheno_cnt_skat, 'skat'), (gene_sig.sig_pheno_cnt_burden, 'burden')]).explode('data')
gene_sig = gene_sig.transmute(sig_cnt=gene_sig.data[0], result_type=gene_sig.data[1])
gene_sig.export('gs://ukbb-exome-public/summary_statistics_analysis/gene_sig300k.txt.bgz')

# Phenotype Association Count (Gene-Level)
pheno_sig = gene.cols()
pheno_sig = pheno_sig.select('n_cases', 'description', data=[(pheno_sig.sig_gene_cnt_skato, 'skato'), (pheno_sig.sig_gene_cnt_skat, 'skat'), (pheno_sig.sig_gene_cnt_burden, 'burden')]).explode('data')
pheno_sig = pheno_sig.transmute(sig_cnt=pheno_sig.data[0], result_type=pheno_sig.data[1])
pheno_sig = pheno_sig.annotate(trait_type2=hl.if_else(pheno_sig.trait_type == 'icd_first_occurrence', 'icd10', pheno_sig.trait_type))
pheno_sig.export('gs://ukbb-exome-public/summary_statistics_analysis/pheno_sig300k.txt.bgz')

# Variant Association Count
var_sig = var.rows()
var_sig.export('gs://ukbb-exome-public/summary_statistics_analysis/var_sig300k.txt.bgz')

# Phenotype Association Count  (Variant-Level)
pheno_var_sig = var.cols().select('n_cases', 'description',  'sig_var_cnt')
pheno_var_sig = pheno_var_sig.annotate(trait_type2=hl.if_else(pheno_var_sig.trait_type == 'icd_first_occurrence', 'icd10', pheno_var_sig.trait_type))
pheno_var_sig.export('gs://ukbb-exome-public/summary_statistics_analysis/pheno_var_sig300k.txt.bgz')

# Phenotype Association by Annotation - Gene
pheno_annt = get_sig_cnt_annt_ht(phenos_to_keep=phenos_to_keep)
pheno_annt = pheno_annt.annotate(trait_type2=hl.if_else(pheno_annt.trait_type == 'icd_first_occurrence', 'icd10', pheno_annt.trait_type),)
pheno_annt.export('gs://ukbb-exome-public/summary_statistics_analysis/pheno_annt_gene300k.txt.bgz')

# Phenotype Association by Annotation - Variant
pheno_var_annt = get_sig_cnt_annt_ht(result_type='variant', phenos_to_keep=phenos_to_keep)
pheno_var_annt = pheno_var_annt.annotate(trait_type2=hl.if_else(pheno_var_annt.trait_type == 'icd_first_occurrence', 'icd10', pheno_var_annt.trait_type))
pheno_var_annt.export('gs://ukbb-exome-public/summary_statistics_analysis/pheno_annt_var300k.txt.bgz')

# Variant-Gene Association Count Comparison
var_gene = compare_gene_var_sig_cnt_ht(phenos_to_keep=phenos_to_keep, genes_to_keep=genes_to_keep)
var_gene = var_gene.select('n_cases', 'description', data=[(var_gene.gene_var_sig_cnt, 'Gene and Variant'), (var_gene.var_sig_cnt, 'Variant Only'), (var_gene.var_sig_cnt, 'Gene Only')]).explode('data')
var_gene = var_gene.transmute(sig_cnt=var_gene.data[0], cnt_type=var_gene.data[1])
var_gene = var_gene.annotate(trait_type2=hl.if_else(var_gene.trait_type == 'icd_first_occurrence', 'icd10', var_gene.trait_type))
var_gene.export('gs://ukbb-exome-public/summary_statistics_analysis/var_gene_comparison300k.txt.bgz')

# --------BETAs for significant hits (Singel-Variant Test)----------
var = hl.read_matrix_table(get_results_mt_path('variant'))
var = var.annotate_entries(BETA_sig = hl.or_missing(var.Pvalue<1e-6, var.BETA))
var = var.select_cols('n_cases')
var = var.annotate_cols(rare_beta = hl.agg.filter(var.AF < 0.001, hl.agg.mean(var.BETA_sig)),
                        common_beta = hl.agg.filter(var.AF >= 0.001, hl.agg.mean(var.BETA_sig)),)
var.cols().export('gs://ukbb-exome-pharma-wlu/subset_data/mean_beta_sig_var300k.txt.bgz')