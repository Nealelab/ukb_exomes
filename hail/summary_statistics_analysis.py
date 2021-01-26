import hail as hl
from ukb_exomes import *

# Load and Combine Data
var = hl.read_matrix_table(get_results_mt_path('variant'))
var = var.annotate_rows(sig_pheno_cnt=hl.agg.sum(var.Pvalue < 1e-6))
var = var.annotate_cols(sig_var_cnt=hl.agg.sum(var.Pvalue < 1e-6))
gene = get_sig_cnt_mt()

# Gene Association Count
gene_sig = gene.rows()
gene_skato = gene_sig.select('caf', sig_cnt=gene_sig.sig_pheno_cnt_skato, result_type='SKATO')
gene_skat = gene_sig.select('caf', sig_cnt=gene_sig.sig_pheno_cnt_skat, result_type='SKAT')
gene_burden = gene_sig.select('caf', sig_cnt=gene_sig.sig_pheno_cnt_burden, result_type='Burden Test')
gene_sig = gene_skato.union(gene_skat, gene_burden)
gene_sig.export('gs://ukbb-exome-public/summary_statistics_analysis/gene_sig300k.txt.bgz')

# Phenotype Association Count (Gene-Level)
pheno_sig = gene.cols()
pheno_skato = pheno_sig.select('n_cases', 'description', sig_cnt=pheno_sig.sig_gene_cnt_skato, result_type='SKATO')
pheno_skat = pheno_sig.select('n_cases', 'description', sig_cnt=pheno_sig.sig_gene_cnt_skat, result_type='SKAT')
pheno_burden = pheno_sig.select('n_cases', 'description', sig_cnt=pheno_sig.sig_gene_cnt_burden, result_type='Burden Test')
pheno_sig = pheno_skato.union(pheno_skat, pheno_burden)
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
p_skato = get_sig_cnt_annt_ht()
p_skat = get_sig_cnt_annt_ht(test_type='skat')
p_burden = get_sig_cnt_annt_ht(test_type='burden')
s_skato = get_sig_cnt_annt_ht(annotation='synonymous')
s_skat = get_sig_cnt_annt_ht(test_type='skat', annotation='synonymous')
s_burden = get_sig_cnt_annt_ht(test_type='burden', annotation='synonymous')
m_skato = get_sig_cnt_annt_ht(annotation='missense|LC')
m_skat = get_sig_cnt_annt_ht(test_type='skat', annotation='missense|LC')
m_burden = get_sig_cnt_annt_ht(test_type='burden', annotation='missense|LC')
pheno_annt = p_skato.union(p_skat, p_burden, s_skato, s_skat, s_burden, m_skato, m_skat, m_burden)
pheno_annt = pheno_annt.annotate(trait_type2=hl.if_else(pheno_annt.trait_type == 'icd_first_occurrence', 'icd10', pheno_annt.trait_type),
                                 result_type=hl.if_else(pheno_annt.result_type == 'BURDEN', 'Burden Test', pheno_annt.result_type))
pheno_annt.export('gs://ukbb-exome-public/summary_statistics_analysis/pheno_annt_gene300k.txt.bgz')

# Phenotype Association by Annotation - Variant
p_var = get_sig_cnt_annt_ht(result_type='variant', test_type='variant')
s_var = get_sig_cnt_annt_ht(result_type='variant', test_type='variant', annotation='synonymous')
m_var = get_sig_cnt_annt_ht(result_type='variant', test_type='variant', annotation='missense')
pheno_var_annt = p_var.union(s_var, m_var)
pheno_var_annt = pheno_var_annt.annotate(trait_type2=hl.if_else(pheno_var_annt.trait_type == 'icd_first_occurrence', 'icd10', pheno_var_annt.trait_type))
pheno_var_annt = pheno_var_annt.annotate(annotation=hl.if_else(pheno_var_annt.annotation_type == 'missense', 'missense|LC', pheno_var_annt.annotation_type))
pheno_var_annt.export('gs://ukbb-exome-public/summary_statistics_analysis/pheno_annt_var300k.txt.bgz')

# Variant-Gene Association Count Comparison
var_gene = compare_gene_var_sig_cnt_ht()
var_gene1 = var_gene.select('n_cases', 'description', sig_cnt=var_gene.gene_var_sig_cnt, cnt_type='Gene and Variant')
var_gene2 = var_gene.select('n_cases', 'description', sig_cnt=var_gene.var_sig_cnt, cnt_type='Variant Only')
var_gene3 = var_gene.select('n_cases', 'description', sig_cnt=var_gene.var_sig_cnt, cnt_type='Gene Only')
var_gene = var_gene1.union(var_gene2, var_gene3)
var_gene = var_gene.annotate(trait_type2=hl.if_else(var_gene.trait_type == 'icd_first_occurrence', 'icd10', var_gene.trait_type))
var_gene.export('gs://ukbb-exome-public/summary_statistics_analysis/var_gene_comparison300k.txt.bgz')
