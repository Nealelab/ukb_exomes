# ukb_exomes

The results of this analysis are released in Google Cloud bucket `gs://ukbb-exome-public/`:
 - Main summary statistics MatrixTable:
 	+ Variant-level results: `gs://ukbb-exome-public/500k/results/variant_results.mt`
 	+ Gene-level results: `gs://ukbb-exome-public/500k/results/results.mt`
 - QC information annotated MatrixTable or Hail Table:
 	+ Variant-level: `gs://ukbb-exome-public/500k/qc/variant_qc_metrics_ukb_exomes_500k{.mt, .ht}` 
 	+ Gene-level: `gs://ukbb-exome-public/500k/qc/gene_qc_metrics_ukb_exomes_500k{.mt, .ht}` 
 	+ Phenotype: `gs://ukbb-exome-public/500k/qc/pheno_qc_metrics_ukb_exomes_500k.ht`
 	
 We also provide the following derived datasets for convenience:
 - Gene-annotation group cumulative allele frequency table: `gs://ukbb-exome-public/500k/qc/gene_caf_500k.ht`
 	

These files can be accessed by cloning this and the https://github.com/broadinstitute/ukbb_qc repo, import the `ukbb_common` python module and accessing them programmatically. We recommend using these functions, as they apply our QC metrics and include convenience metrics such as lambda GC.

```
%%bash
git clone https://github.com/broadinstitute/ukbb_qc
git clone https://github.com/Nealelab/ukb_exomes
```

```
from ukb_exomes import *
from ukbb_common import *
```

To read the original MatrixTables with 4529 phenotypes:
```
## Gene-level results
gene_mt = hl.read_matrix_table(get_results_mt_path(result_type='gene'))
## Variant-level results
var_mt = hl.read_matrix_table(get_results_mt_path(result_type='variant'))
```

To read the full MatrixTables with QC information annotated:
```
## Gene-level results
gene_mt = load_final_sumstats_table(result_type='gene', extension="mt")
## Variant-level results
var_mt = load_final_sumstats_table(result_type='variant', extension="mt")
```

To get the final QCed MatrixTables (Note that we have two options for `test_type`: `skato` and `burden`, which indicates which test the lambda GC used here were computed from):
```
## Gene-level results
gene_mt = get_qc_result_mt(result_type="gene", test_type="skato")
## Variant-level results
var_mt = get_qc_result_mt(result_type="variant", test_type="skato")
```

The basic summary statistics of the gene-based tests have the following schema:
```
----------------------------------------
Global fields:
    'coverage_min': int32
    'expected_AC_min': int32
    'n_var_min': int32
    'gene_syn_lambda_min': float64
    'pheno_lambda_min': float64
----------------------------------------
Column fields:
    'n_cases': int32
    'n_controls': int32
    'heritability': float64
    'saige_version': str
    'inv_normalized': str
    'trait_type': str
    'phenocode': str
    'pheno_sex': str
    'coding': str
    'modifier': str
    'n_cases_defined': int64
    'n_cases_both_sexes': int64
    'n_cases_females': int64
    'n_cases_males': int64
    'description': str
    'description_more': str
    'coding_description': str
    'category': str
    'expected_ac_col_filter': int64
    'lambda_gc_skat': float64
    'lambda_gc_burden': float64
    'lambda_gc_skato': float64
    'keep_pheno_skato': bool
    'keep_pheno_skat': bool
    'keep_pheno_burden': bool
    'keep_pheno_unrelated': bool
----------------------------------------
Row fields:
    'gene_id': str
    'gene_symbol': str
    'annotation': str
    'interval': interval<locus<GRCh38>>
    'markerIDs': str
    'markerAFs': str
    'total_variants': int32
    'Nmarker_MACCate_1': int32
    'Nmarker_MACCate_2': int32
    'Nmarker_MACCate_3': int32
    'Nmarker_MACCate_4': int32
    'Nmarker_MACCate_5': int32
    'Nmarker_MACCate_6': int32
    'Nmarker_MACCate_7': int32
    'Nmarker_MACCate_8': int32
    'CAF': float64
    'mean_coverage': float64
    'expected_ac_row_filter': int64
    'continuous_lambda_gc_skato': float64
    'continuous_lambda_gc_skat': float64
    'continuous_lambda_gc_burden': float64
    'categorical_lambda_gc_skato': float64
    'categorical_lambda_gc_skat': float64
    'categorical_lambda_gc_burden': float64
    'icd10_lambda_gc_skato': float64
    'icd10_lambda_gc_skat': float64
    'icd10_lambda_gc_burden': float64
    'annotation_lambda_gc_skato': float64
    'annotation_lambda_gc_skat': float64
    'annotation_lambda_gc_burden': float64
    'synonymous_lambda_gc_skato': float64
    'synonymous_lambda_gc_skat': float64
    'synonymous_lambda_gc_burden': float64
    'keep_gene_skato': bool
    'keep_gene_skat': bool
    'keep_gene_burden': bool
    'keep_gene_coverage': bool
    'keep_gene_expected_ac': bool
    'keep_gene_n_var': bool
----------------------------------------
Entry fields:
    'Pvalue': float64
    'Pvalue_Burden': float64
    'Pvalue_SKAT': float64
    'BETA_Burden': float64
    'SE_Burden': float64
    'Pvalue.NA': float64
    'Pvalue_Burden.NA': float64
    'Pvalue_SKAT.NA': float64
    'BETA_Burden.NA': float64
    'SE_Burden.NA': float64
    'total_variants_pheno': int32
    'expected_AC': float64
    'keep_entry_expected_ac': bool
----------------------------------------
Column key: ['trait_type', 'phenocode', 'pheno_sex', 'coding', 'modifier']
Row key: ['gene_id', 'gene_symbol', 'annotation']
----------------------------------------
```

The basic summary statistics of the variant-level tests have the following schema:
```
----------------------------------------
Global fields:
    'expected_AC_min': int32
    'pheno_lambda_min': float64
----------------------------------------
Column fields:
    'n_cases': int32
    'n_controls': int32
    'heritability': float64
    'saige_version': str
    'inv_normalized': str
    'trait_type': str
    'phenocode': str
    'pheno_sex': str
    'coding': str
    'modifier': str
    'n_cases_defined': int64
    'n_cases_both_sexes': int64
    'n_cases_females': int64
    'n_cases_males': int64
    'description': str
    'description_more': str
    'coding_description': str
    'category': str
    'expected_ac_col_filter': int64
    'lambda_gc_skat': float64
    'lambda_gc_burden': float64
    'lambda_gc_skato': float64
    'keep_pheno_skato': bool
    'keep_pheno_skat': bool
    'keep_pheno_burden': bool
    'keep_pheno_unrelated': bool
----------------------------------------
Row fields:
    'locus': locus<GRCh38>
    'alleles': array<str>
    'markerID': str
    'gene': str
    'annotation': str
    'call_stats': struct {
        AC: int32,
        AF: float64,
        AN: int32,
        homozygote_count: int32
    }
    'expected_ac_row_filter': int64
    'keep_var_expected_ac': bool
    'keep_var_annt': bool
----------------------------------------
Entry fields:
    'AC': int32
    'AF': float64
    'BETA': float64
    'SE': float64
    'AF.Cases': float64
    'AF.Controls': float64
    'Pvalue': float64
    'expected_AC': float64
    'keep_entry_expected_ac': bool
----------------------------------------
Column key: ['trait_type', 'phenocode', 'pheno_sex', 'coding', 'modifier']
Row key: ['locus', 'alleles']
----------------------------------------
```
