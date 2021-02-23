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
# Gene-based Tests
lambda_gene_full = compute_lambda_gc_ht(result_type='gene', af_lower=1e-4, n_var_min=2, coverage_min=20)
lambda_gene_full = lambda_gene_full.annotate(trait_type2=hl.if_else(lambda_gene_full.trait_type == 'icd_first_occurrence', 'icd10', lambda_gene_full.trait_type))
lambda_gene_full.export('gs://ukbb-exome-public/summary_statistics_analysis/lambda_gene_full300k_filtered_coverage.txt.bgz')

lambda_gene_annt = compute_lambda_gc_ht(by_annotation=True, af_lower=1e-4, n_var_min=2, coverage_min=20)
lambda_gene_annt = lambda_gene_annt.annotate(trait_type2=hl.if_else(lambda_gene_annt.trait_type == 'icd_first_occurrence', 'icd10', lambda_gene_annt.trait_type))
lambda_gene_annt.export('gs://ukbb-exome-public/summary_statistics_analysis/lambda_gene_annt300k_filtered_coverage.txt.bgz')

g1e5 = compute_lambda_gc_ht(af_upper=0.0001, coverage_min=20)
g1e4 = compute_lambda_gc_ht(af_lower=0.0001, af_upper=0.001, coverage_min=20)
g1e3 = compute_lambda_gc_ht(af_lower=0.001, af_upper=0.01, coverage_min=20)
g1e2 = compute_lambda_gc_ht(af_lower=0.01, af_upper=0.1, coverage_min=20)
g1e1 = compute_lambda_gc_ht(af_lower=0.1, coverage_min=20)
lambda_gene_af = g1e5.union(g1e4, g1e3, g1e2, g1e1)
lambda_gene_af = lambda_gene_af.annotate(trait_type2=hl.if_else(lambda_gene_af.trait_type == 'icd_first_occurrence', 'icd10', lambda_gene_af.trait_type),)
lambda_gene_af.export('gs://ukbb-exome-public/summary_statistics_analysis/lambda_gene_af300k_filtered_coverage.txt.bgz')

g1e4 = compute_lambda_gc_ht(by_annotation=True, af_lower=0.0001, af_upper=0.001, n_var_min=2)
g1e3 = compute_lambda_gc_ht(by_annotation=True, af_lower=0.001, af_upper=0.01, n_var_min=2)
g1e2 = compute_lambda_gc_ht(by_annotation=True, af_lower=0.01, af_upper=0.1, n_var_min=2)
g1e1 = compute_lambda_gc_ht(by_annotation=True, af_lower=0.1, n_var_min=2)
lambda_gene_af_annt = g1e4.union(g1e3, g1e2, g1e1)
lambda_gene_af_annt = lambda_gene_af_annt.annotate(trait_type2=hl.if_else(lambda_gene_af_annt.trait_type == 'icd_first_occurrence', 'icd10', lambda_gene_af_annt.trait_type),)
lambda_gene_af_annt.export('gs://ukbb-exome-public/summary_statistics_analysis/lambda_gene_af_annt300k_filtered.txt.bgz')

lambda_by_gene = compute_lambda_gc_ht(by_gene = True, af_lower=1e-4, n_var_min=2, coverage_min=20)
lambda_by_gene.export('gs://ukbb-exome-public/summary_statistics_analysis/lambda_by_gene300k_filtered_coverage.txt.bgz')

# Single-Variant Results
lambda_var_full = compute_lambda_gc_ht(result_type='variant', af_lower=2e-5)
lambda_var_full = lambda_var_full.annotate(trait_type2=hl.if_else(lambda_var_full.trait_type == 'icd_first_occurrence', 'icd10', lambda_var_full.trait_type))
lambda_var_full.export('gs://ukbb-exome-public/summary_statistics_analysis/lambda_var_full300k_filtered.txt.bgz')

lambda_var_annt = compute_lambda_gc_ht(result_type='variant', by_annotation=True, af_lower=2e-5)
lambda_var_annt = lambda_var_annt.annotate(trait_type2=hl.if_else(lambda_var_annt.trait_type == 'icd_first_occurrence', 'icd10', lambda_var_annt.trait_type))
lambda_var_annt.export('gs://ukbb-exome-public/summary_statistics_analysis/lambda_var_annt300k_filtered.txt.bgz')

v1e5 = compute_lambda_gc_ht(result_type='variant', af_lower=2e-5, af_upper=0.0001)
v1e4 = compute_lambda_gc_ht(result_type='variant', af_lower=0.0001, af_upper=0.001)
v1e3 = compute_lambda_gc_ht(result_type='variant', af_lower=0.001, af_upper=0.01)
v1e2 = compute_lambda_gc_ht(result_type='variant', af_lower=0.01, af_upper=0.1)
v1e1 = compute_lambda_gc_ht(result_type='variant', af_lower=0.1)
lambda_var_af = v1e5.union(v1e4, v1e3, v1e2, v1e1)
lambda_var_af = lambda_var_af.annotate(trait_type2=hl.if_else(lambda_var_af.trait_type == 'icd_first_occurrence', 'icd10', lambda_var_af.trait_type),)
lambda_var_af.export('gs://ukbb-exome-public/summary_statistics_analysis/lambda_var_af300k_filtered.txt.bgz')

v1e5 = compute_lambda_gc_ht(result_type='variant', by_annotation=True, af_lower=2e-5, af_upper=0.0001)
v1e4 = compute_lambda_gc_ht(result_type='variant', by_annotation=True, af_lower=0.0001, af_upper=0.001)
v1e3 = compute_lambda_gc_ht(result_type='variant', by_annotation=True, af_lower=0.001, af_upper=0.01)
v1e2 = compute_lambda_gc_ht(result_type='variant', by_annotation=True, af_lower=0.01, af_upper=0.1)
v1e1 = compute_lambda_gc_ht(result_type='variant', by_annotation=True, af_lower=0.1)
lambda_var_af_annt = v1e5.union(v1e4, v1e3, v1e2, v1e1)
lambda_var_af_annt = lambda_var_af_annt.annotate(trait_type2=hl.if_else(lambda_var_af_annt.trait_type == 'icd_first_occurrence', 'icd10', lambda_var_af_annt.trait_type),)
lambda_var_af_annt.export('gs://ukbb-exome-public/summary_statistics_analysis/lambda_var_af_annt300k_filtered.txt.bgz')

# Filter by Expected AC
mt = hl.read_matrix_table(get_results_mt_path('variant'))
mt = mt.annotate_entries(expected_AC = mt.AF*mt.n_cases)
mt = mt.select_cols('n_cases',
                    lambda_gc1 = hl.agg.filter(mt.expected_AC <= 1, hl.methods.statgen._lambda_gc_agg(mt.Pvalue)),
                    lambda_gc2 = hl.agg.filter((mt.expected_AC > 1) & (mt.expected_AC <= 10), hl.methods.statgen._lambda_gc_agg(mt.Pvalue)),
                    lambda_gc3 = hl.agg.filter((mt.expected_AC > 10) & (mt.expected_AC <= 100), hl.methods.statgen._lambda_gc_agg(mt.Pvalue)),
                    lambda_gc4 = hl.agg.filter((mt.expected_AC > 100) & (mt.expected_AC <= 1000), hl.methods.statgen._lambda_gc_agg(mt.Pvalue)),
                    lambda_gc5 = hl.agg.filter((mt.expected_AC > 1000) & (mt.expected_AC <= 10000), hl.methods.statgen._lambda_gc_agg(mt.Pvalue)),
                    lambda_gc6 = hl.agg.filter((mt.expected_AC > 10000) & (mt.expected_AC <= 100000), hl.methods.statgen._lambda_gc_agg(mt.Pvalue)),
                    lambda_gc7 = hl.agg.filter(mt.expected_AC > 100000, hl.methods.statgen._lambda_gc_agg(mt.Pvalue)),
                    trait_type2=hl.if_else(mt.trait_type == 'icd_first_occurrence', 'icd10', mt.trait_type))
mt.cols().export('gs://ukbb-exome-pharma-wlu/subset_data/lambda_var_expectedAC300k.txt.bgz')

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

# --------Annotate BLOSUM, polyphen, pext information to single-variant table----------
rg37 = hl.get_reference('GRCh37')
rg38 = hl.get_reference('GRCh38')
rg37.add_liftover('gs://hail-common/references/grch37_to_grch38.over.chain.gz', rg38)
ht = hl.read_matrix_table('gs://gcp-public-data--gnomad/papers/2019-tx-annotation/pre_computed/all.possible.snvs.tx_annotated.021520.ht').rows()
ht = ht.explode(ht.tx_annotation)
ht = ht.annotate(new_locus=hl.liftover(ht.locus, 'GRCh38'))
ht = ht.filter(hl.is_defined(ht.new_locus))
ht = ht.key_by(locus=ht.new_locus, alleles = ht.alleles)

vep = hl.read_table(var_annotations_ht_path('vep', *TRANCHE_DATA['300k']))
vep = process_consequences(vep)
vep = vep.explode(vep.vep.worst_csq_by_gene_canonical)

var = get_sig_cnt_mt('variant').rows()
var = var.annotate(annotation=hl.if_else(hl.literal({'missense', 'LC'}).contains(var.annotation), 'missense|LC', var.annotation),
                   polyphen2 = vep[var.key].vep.worst_csq_by_gene_canonical.polyphen_prediction,
                   amino_acids = vep[var.key].vep.worst_csq_by_gene_canonical.amino_acids,
                   lof = vep[var.key].vep.worst_csq_by_gene_canonical.lof,
                   most_severe_consequence = vep[var.key].vep.worst_csq_by_gene_canonical.most_severe_consequence,
                   tx_annotation_csq = vep[var.key].tx_annotation_csq,
                   mean_proportion = vep[var.key].mean_proportion)
var.export('gs://ukbb-exome-public/summary_statistics_analysis/var_info300k.txt.bgz')
