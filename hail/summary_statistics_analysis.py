import hail as hl
import argparse
from ukb_exomes import *


def main(args):
    if args.compute_lambdas:
        if args.lambda_without_filters:
            if args.gene_results:
                write_lambda_hts(result_type='gene', extension=args.extension, random_phenos=args.random_phenos, overwrite=args.overwrite)

            if args.variant_results:
                write_lambda_hts(result_type='var', extension=args.extension, random_phenos=args.random_phenos, overwrite=args.overwrite)

        else:
            if args.gene_results:
                write_lambda_hts(result_type='gene', freq_lower=args.caf_lower, n_var_min=args.n_var_min, coverage_min=args.coverage_min, extension=args.extension, overwrite=args.overwrite)

            if args.variant_results:
                write_lambda_hts(result_type='var', freq_lower=args.af_lower, var_filter=True, extension=args.extension, overwrite=args.overwrite)

    if args.get_sig_cnts:
        filter_flag = '_filtered'
        if args.sig_without_filters:
            phenos_to_keep = None
            genes_to_keep = None
            filter_flag = ''
        else:
            lambda_gene = hl.read_table(get_ukb_exomes_sumstat_path(subdir=lambda_folder, dataset='lambda_full_filtered'))
            lambda_by_gene = hl.read_table(get_ukb_exomes_sumstat_path(subdir=lambda_folder, dataset='lambda_by_gene_filtered', result_type=''))

            lambda_by_gene_filtered = get_lambda_filtered_ht(lambda_by_gene, lambda_name=f'all_lambda_gc_{args.result_type.lower()}', lower=args.lambda_lower, upper=args.lambda_upper)
            gene_ids_to_keep_syn = lambda_by_gene_filtered.filter(lambda_by_gene_filtered.annotation == 'synonymous')
            lambda_by_gene = lambda_by_gene.key_by('gene_id', 'gene_symbol')
            gene_ids_to_keep_syn = gene_ids_to_keep_syn.key_by('gene_id', 'gene_symbol')
            genes_to_keep = lambda_by_gene.filter(hl.is_defined(gene_ids_to_keep_syn.index(lambda_by_gene.key)))
            genes_to_keep = genes_to_keep.key_by('gene_id', 'gene_symbol', 'annotation')

            phenos_to_keep = get_lambda_filtered_ht(lambda_gene, lambda_name=f'lambda_gc_{args.result_type.lower()}', lower=args.lambda_lower, upper=args.lambda_upper)
            phenos_to_remove = get_corr_phenos_ht(r_2=args.r2_cut, tie_breaker=more_cases_tie_breaker)
            phenos_to_keep = phenos_to_keep.filter(hl.is_defined(phenos_to_remove.key_by(trait_type=phenos_to_remove.node.trait_type, phenocode=phenos_to_remove.node.phenocode,
                                                                                         pheno_sex=phenos_to_remove.node.pheno_sex, coding=phenos_to_remove.node.coding,
                                                                                         modifier=phenos_to_remove.node.modifier, )[phenos_to_keep.key]), keep=False)
            genes_to_keep.write(get_ukb_exomes_sumstat_path(subdir=lambda_folder, dataset='genes_to_keep', result_type=''), overwrite=args.overwrite)
            phenos_to_keep.write(get_ukb_exomes_sumstat_path(subdir=lambda_folder, dataset='phenos_to_keep', result_type=''), overwrite=args.overwrite)

        if args.get_p_values:
            if args.random_phenos:
                get_pvalue_by_freq_interval_ht(result_type='gene', test_type='skato', random_phenos=args.random_phenos
                                               ).write(get_ukb_exomes_sumstat_path(subdir=pvalue_folder, dataset='rp_p_value_skato', result_type='gene',
                                                                                   extension=args.extension), overwrite=args.overwrite)
                get_pvalue_by_freq_interval_ht(result_type='gene', test_type='burden', random_phenos=args.random_phenos
                                               ).write(get_ukb_exomes_sumstat_path(subdir=pvalue_folder, dataset='rp_p_value_burden', result_type='gene',
                                                                                   extension=args.extension), overwrite=args.overwrite)
                get_pvalue_by_freq_interval_ht(result_type='variant', random_phenos=args.random_phenos
                                               ).write(get_ukb_exomes_sumstat_path(subdir=pvalue_folder, dataset='rp_p_value', result_type='var',
                                                                                   extension=args.extension), overwrite=args.overwrite)
            else:
                get_pvalue_by_freq_interval_ht(result_type='gene', test_type=args.result_type, phenos_to_keep=phenos_to_keep, n_var_min=args.n_var_min, coverage_min=args.coverage_min
                                               ).write(get_ukb_exomes_sumstat_path(subdir=pvalue_folder, dataset=f'pvalue_{args.result_type}_{filter_flag}',
                                                                                   result_type='gene', extension=args.extension), overwrite=args.overwrite)

        if args.compare_var_gene:
            var_gene = compare_gene_var_sig_cnt_mt(test_type = args.result_type, phenos_to_keep=phenos_to_keep, genes_to_keep=genes_to_keep, var_min_freq=2e-5)
            var_gene.rows().write(get_ukb_exomes_sumstat_path(subdir=sig_folder, dataset=f'var_gene_comparison_by_gene{filter_flag}', result_type='', extension=args.extension))

            var_gene_by_pheno = var_gene.cols()
            var_gene_by_pheno = var_gene_pheno.select('n_cases', 'description', 'trait_type2', 
                                                      data=[(var_gene.gene_var_sig_cnt, 'Gene and Variant'), 
                                                            (var_gene.var_sig_cnt, 'Variant Only'), 
                                                            (var_gene.var_sig_cnt, 'Gene Only')]).explode('data')
            var_gene_by_pheno = var_gene_pheno.transmute(sig_cnt=var_gene_pheno.data[0], cnt_type=var_gene_pheno.data[1])
            var_gene_pheno.write(get_ukb_exomes_sumstat_path(subdir=sig_folder, dataset=f'var_gene_comparison_by_pheno{filter_flag}', result_type='', extension=args.extension))

        if args.get_icd_pvalue:
            get_icd_min_p_ht('var', phenos_to_keep=phenos_to_keep, var_min_freq=args.af_lower
                             ).write(get_ukb_exomes_sumstat_path(subdir=pvalue_folder, dataset=f'icd_min_pvalue{filter_flag}', result_type='var',
                                                                 extension=args.extension), overwrite=args.overwrite)
            get_icd_min_p_ht('gene', test_type = args.result_type.lower(), phenos_to_keep = phenos_to_keep, genes_to_keep = genes_to_keep
                             ).write(get_ukb_exomes_sumstat_path(subdir=pvalue_folder, dataset=f'icd_min_pvalue{filter_flag}', result_type='gene',
                                        extension=args.extension), overwrite=args.overwrite)

        if args.gene_results:
            gene = get_sig_cnt_mt(test_type=args.result_type.lower(), phenos_to_keep=phenos_to_keep, genes_to_keep=genes_to_keep)
            gene_sig = gene.rows()
            gene_sig.write(get_ukb_exomes_sumstat_path(subdir=sig_folder, dataset=f'sig_cnt{filter_flag}', result_type='gene', extension=args.extension), overwrite=args.overwrite)

            pheno_sig = gene.cols()
            pheno_sig.write(get_ukb_exomes_sumstat_path(subdir=sig_folder, dataset=f'pheno_sig_cnt{filter_flag}', result_type='gene', extension=args.extension), overwrite=args.overwrite)

        if args.variant_results:
            var = get_sig_cnt_mt('variant', phenos_to_keep=phenos_to_keep, var_min_freq=args.af_lower)
            if args.add_variant_info:
                var = annotate_additional_info_ht(result_mt=var, result_type='variant')
            var_sig = var.rows()
            var_sig.write(get_ukb_exomes_sumstat_path(subdir=sig_folder, dataset=f'sig_cnt{filter_flag}', result_type='var', extension=args.extension), overwrite=args.overwrite)

            pheno_var_sig = var.cols()
            pheno_var_sig.write(get_ukb_exomes_sumstat_path(subdir=sig_folder, dataset=f'pheno_sig_cnt{filter_flag}', result_type='var', extension=args.extension), overwrite=args.overwrite)

    if args.get_sig_betas:
        if args.variant_results:
            var = hl.read_matrix_table(get_results_mt_path('variant'))
            var = var.annotate_entries(BETA_sig=hl.or_missing(var.Pvalue < LEVELS['variant'], var.BETA))
            var = var.select_cols('n_cases')
            var = var.annotate_cols(rare_beta=hl.agg.filter(var.AF < 0.001, hl.agg.mean(var.BETA_sig)), 
                                    common_beta=hl.agg.filter(var.AF >= 0.001, hl.agg.mean(var.BETA_sig)), )
            var.cols().write(get_ukb_exomes_sumstat_path(subdir=beta_folder, dataset='var_beta', result_type='var', extension=args.extension), overwrite=args.overwrite)
        if args.gene_results:
            gene = compare_gene_var_sig_cnt_mt(test_type='burden', tranche=CURRENT_TRANCHE)
            gene = annotate_additional_info_ht(result_mt=gene, result_type='gene')
            gene = gene.annotate_rows(lambda_gc_gene=hl.agg.filter(hl.is_defined(gene.Pvalue_Burden), hl.methods.statgen._lambda_gc_agg(gene.Pvalue_Burden)))
            gene = gene.annotate_cols(lambda_gc_pheno=hl.agg.filter(hl.is_defined(gene.Pvalue_Burden), hl.methods.statgen._lambda_gc_agg(gene.Pvalue_Burden)))
            gene = gene.annotate_entries(BETA_sig=hl.or_missing(gene.Pvalue_Burden < LEVELS['burden'], gene.BETA_Burden))
            gene = gene.entries()
            gene = gene.filter(hl.is_defined(gene.BETA_sig))
            gene = gene.select('CAF', 'mean_coverage', 'total_variants', 'description', 'description_more', 'coding_description',
                               'lambda_gc_gene', 'lambda_gc_pheno', 'BETA_sig', 'Pvalue_Burden', 'var_cnt')

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--overwrite', help='Overwrite everything', action='store_true')
    parser.add_argument('--compute_lambdas', help='Compute lambda GC', action='store_true')
    parser.add_argument('--random_phenos', help='Compute lambda GC for random phenotypes', action='store_true')
    parser.add_argument('--lambda_without_filters', help='Compute lambda GC without any filters', action='store_true')
    parser.add_argument('--get_p_values', help='Generate P value table', action='store_true')
    parser.add_argument('--get_sig_cnts', help='Count the number of significant associations', action='store_true')
    parser.add_argument('--sig_without_filters', help='Get significant associations without any filters', action='store_true')
    parser.add_argument('--add_variant_info', help='Add extra annotation information on variants', action='store_true')
    parser.add_argument('--get_sig_betas', help='Get effect sizes for significant associations', action='store_true')
    parser.add_argument('--gene_results', help='Use gene-level results', action='store_true')
    parser.add_argument('--variant_results', help='Use single-variant results', action='store_true')
    parser.add_argument('--get_icd_pvalue', help='Get the minimum p-value for each ICD group', action='store_true')
    parser.add_argument('--compare_var_gene', help='Compare significant associations between genes and variants', action='store_true')
    parser.add_argument('--result_type', help='Test results to apply lambda filters on: skato OR burden', nargs='?', const='skato', type=str)
    parser.add_argument('--extension', help='Output file format', nargs='?', const='ht', type=str)
    parser.add_argument('--coverage_min', help='Keep genes with higher coverage', nargs='?', const=20, type=int)
    parser.add_argument('--n_var_min', help='Keep genes with larger number of variants', nargs='?', const=2, type=int)
    parser.add_argument('--caf_lower', help='Keep genes with higher cumulative allele frequency', nargs='?', const=1e-4, type=float)
    parser.add_argument('--af_lower', help='Keep variants with higher allele frequency', nargs='?', const=2e-5, type=float)
    parser.add_argument('--lambda_lower', help='Remove genes/phenotypes with lower lambda value', nargs='?', const=0.75, type=float)
    parser.add_argument('--lambda_upper', help='Remove genes/phenotypes with higher lambda value', nargs='?', const=1.5, type=float)
    parser.add_argument('--r2_cut', help='Remove related phenotypes with higher correlation', nargs='?', const=0.5, type=float)
    args = parser.parse_args()

    main(args)