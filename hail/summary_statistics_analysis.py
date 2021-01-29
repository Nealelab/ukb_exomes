import hail as hl
from ukb_exomes import *

# Load and Combine Data
var = hl.read_matrix_table(get_results_mt_path('variant'))
var = var.annotate_rows(sig_pheno_cnt=hl.agg.sum(var.Pvalue < 1e-6))
var = var.annotate_cols(sig_var_cnt=hl.agg.sum(var.Pvalue < 1e-6))
gene = get_sig_cnt_mt()

# Gene Association Count
gene_sig = gene.rows()
gene_sig = gene_sig.select('caf', data=[(gene_sig.sig_pheno_cnt_skato, 'SKATO'), (gene_sig.sig_pheno_cnt_skat, 'SKAT'), (gene_sig.sig_pheno_cnt_burden, 'Burden Test')]).explode('data')
gene_sig = gene_sig.transmute(sig_cnt=gene_sig.data[0], result_type=gene_sig.data[1])
gene_sig.export('gs://ukbb-exome-public/summary_statistics_analysis/gene_sig300k.txt.bgz')

# Phenotype Association Count (Gene-Level)
pheno_sig = gene.cols()
pheno_sig = pheno_sig.select('n_cases', 'description', data=[(pheno_sig.sig_gene_cnt_skato, 'SKATO'), (pheno_sig.sig_gene_cnt_skat, 'SKAT'), (pheno_sig.sig_gene_cnt_burden, 'Burden Test')]).explode('data')
pheno_sig = pheno_sig.transmute(sig_cnt=pheno_sig.data[0], result_type=pheno_sig.data[1])
pheno_sig = pheno_sig.annotate(trait_type2=hl.if_else(pheno_sig.trait_type == 'icd_first_occurrence', 'icd10', pheno_sig.trait_type))
pheno_sig.export('gs://ukbb-exome-public/summary_statistics_analysis/pheno_sig300k.txt.bgz')

# Variant Association Count
var_sig = var.rows()
var_sig = var_sig.annotate(annotation=hl.if_else(var_sig.annotation == 'missense', 'missense|LC', var_sig.annotation))
var_sig.export('gs://ukbb-exome-public/summary_statistics_analysis/var_sig300k.txt.bgz')

# Phenotype Association Count  (Variant-Level)
pheno_var_sig = var.cols().select('n_cases', 'description',  'sig_var_cnt')
pheno_var_sig = pheno_var_sig.annotate(trait_type2=hl.if_else(pheno_var_sig.trait_type == 'icd_first_occurrence', 'icd10', pheno_var_sig.trait_type))
pheno_var_sig.export('gs://ukbb-exome-public/summary_statistics_analysis/pheno_var_sig300k.txt.bgz')

# Phenotype Association by Annotation - Gene
pheno_annt = get_sig_cnt_annt_ht()
pheno_annt = pheno_annt.annotate(trait_type2=hl.if_else(pheno_annt.trait_type == 'icd_first_occurrence', 'icd10', pheno_annt.trait_type),)
pheno_annt.export('gs://ukbb-exome-public/summary_statistics_analysis/pheno_annt_gene300k.txt.bgz')

# Phenotype Association by Annotation - Variant
pheno_var_annt = get_sig_cnt_annt_ht(result_type='variant')
pheno_var_annt = pheno_var_annt.annotate(trait_type2=hl.if_else(pheno_var_annt.trait_type == 'icd_first_occurrence', 'icd10', pheno_var_annt.trait_type))
pheno_var_annt.export('gs://ukbb-exome-public/summary_statistics_analysis/pheno_annt_var300k.txt.bgz')

# Variant-Gene Association Count Comparison
var_gene = compare_gene_var_sig_cnt_ht()
var_gene = var_gene.select('n_cases', 'description', data=[(var_gene.gene_var_sig_cnt, 'Gene and Variant'), (var_gene.var_sig_cnt, 'Variant Only'), (var_gene.var_sig_cnt, 'Gene Only')]).explode('data')
var_gene = var_gene.transmute(sig_cnt=var_gene.data[0], cnt_type=var_gene.data[1])
var_gene = var_gene.annotate(trait_type2=hl.if_else(var_gene.trait_type == 'icd_first_occurrence', 'icd10', var_gene.trait_type))
var_gene.export('gs://ukbb-exome-public/summary_statistics_analysis/var_gene_comparison300k.txt.bgz')
