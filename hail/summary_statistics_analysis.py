import sys
print(sys.path)
import hail as hl
from gnomad.resources.grch38.reference_data import clinvar
from ukb_exomes import *

sum_bucket = 'gs://ukbb-exome-public'
sumstats_folder = f'{sum_bucket}/summary_statistics_analysis'
lambda_folder = f'{sumstats_folder}/lambda_gc'
pvalue_folder = f'{sumstats_folder}/pvalue'
sig_folder = f'{sumstats_folder}/sig_cnt'
beta_folder = f'{sumstats_folder}/beta'


def get_ukb_exomes_sumstat_path(subdir:str, dataset:str, result_type: str = 'gene', extension: str = 'ht', tranche: str = CURRENT_TRANCHE):
    return f'{subdir}/{dataset}_{result_type}_{tranche}.{extension}'


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

    if args.get_p_values:
        get_pvalue_by_freq_interval_ht(result_type='gene', test_type='skato', random_phenos=args.random_phenos
                                       ).write(get_ukb_exomes_sumstat_path(subdir=pvalue_folder, dataset='rp_p_value_skato', 
                                                                           result_type='gene', extension=args.extension), overwrite=args.overwrite)
        get_pvalue_by_freq_interval_ht(result_type='gene', test_type='burden', random_phenos=args.random_phenos
                                       ).write(get_ukb_exomes_sumstat_path(subdir=pvalue_folder, dataset='rp_p_value_burden', 
                                                                           result_type='gene', extension=args.extension), overwrite=args.overwrite)
        get_pvalue_by_freq_interval_ht(result_type='variant', random_phenos=args.random_phenos
                                       ).write(get_ukb_exomes_sumstat_path(subdir=pvalue_folder, dataset='rp_p_value', 
                                                                           result_type='var', extension=args.extension), overwrite=args.overwrite)

    if args.get_sig_cnts:
        filter_flag = '_filtered'
        if args.sig_without_filters:
            phenos_to_keep = None
            genes_to_keep = None
            filter_flag = ''

        else:
            lambda_gene = hl.read_table(get_ukb_exomes_sumstat_path(subdir=lambda_folder, dataset='lambda_full_filtered'))
            lambda_by_gene = hl.read_table(get_ukb_exomes_sumstat_path(subdir=lambda_folder, dataset='lambda_by_gene_filtered', result_type=''))

            phenos_to_keep = get_lambda_filtered_ht(lambda_gene, lambda_name=f'lambda_gc_{args.result_type.lower()}', lower=args.lambda_lower, upper=args.lambda_upper)
            genes_to_keep = get_lambda_filtered_ht(lambda_by_gene, lambda_name=f'all_lambda_gc_{args.result_type.lower()}', lower=args.lambda_lower, upper=args.lambda_upper)
            phenos_to_remove = get_corr_phenos_ht(r_2=args.r2_cut, tie_breaker=more_cases_tie_breaker)
            phenos_to_keep = phenos_to_keep.filter(hl.is_defined( phenos_to_remove.key_by(trait_type=phenos_to_remove.node.trait_type, phenocode=phenos_to_remove.node.phenocode,
                                                                                          pheno_sex=phenos_to_remove.node.pheno_sex, coding=phenos_to_remove.node.coding,
                                                                                          modifier=phenos_to_remove.node.modifier, )[phenos_to_keep.key]), keep=False)
            genes_to_keep.write(get_ukb_exomes_sumstat_path(subdir=lambda_folder, dataset='genes_to_keep', result_type=''), overwrite=args.overwrite)
            phenos_to_keep.write(get_ukb_exomes_sumstat_path(subdir=lambda_folder, dataset='phenos_to_keep', result_type=''), overwrite=args.overwrite)

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


        if args.gene_results:
            gene = get_sig_cnt_mt(phenos_to_keep=phenos_to_keep, genes_to_keep=genes_to_keep)
            gene_sig = gene.rows()
            gene_sig.write(get_ukb_exomes_sumstat_path(subdir=sig_folder, dataset=f'sig_cnt{filter_flag}', result_type='gene', extension=args.extension), overwrite=args.overwrite)

            pheno_sig = gene.cols()
            pheno_sig.write(get_ukb_exomes_sumstat_path(subdir=sig_folder, dataset=f'pheno_sig_cnt{filter_flag}', result_type='gene', extension=args.extension), overwrite=args.overwrite)

        if args.variant_results:
            var = get_sig_cnt_mt('variant', phenos_to_keep=phenos_to_keep, var_min_freq=args.af_lower)
            if args.add_variant_info:
                clinvar_ht = clinvar.ht()
                clinvar_ht = clinvar_ht.explode(clinvar_ht.info.CLNSIG)
                clinvar_ht = clinvar_ht.annotate(pathogenicity=hl.case() \
                                                 .when(hl.literal({'Pathogenic', 'Likely_pathogenic', 'Pathogenic/Likely_pathogenic'}).contains(clinvar_ht.info.CLNSIG), 'P/LP') \
                                                 .when(hl.literal({'Uncertain_significance'}).contains(clinvar_ht.info.CLNSIG), 'VUS') \
                                                 .when(hl.literal({'Benign', 'Likely_benign', 'Benign/Likely_benign'}).contains(clinvar_ht.info.CLNSIG), 'B/LB') \
                                                 .or_missing())

                rg37 = hl.get_reference('GRCh37')
                rg38 = hl.get_reference('GRCh38')
                rg37.add_liftover('gs://hail-common/references/grch37_to_grch38.over.chain.gz', rg38)
                ht = hl.read_matrix_table('gs://gcp-public-data--gnomad/papers/2019-tx-annotation/pre_computed/all.possible.snvs.tx_annotated.021520.ht').rows()
                ht = ht.explode(ht.tx_annotation)
                ht = ht.annotate(new_locus=hl.liftover(ht.locus, 'GRCh38'))
                ht = ht.filter(hl.is_defined(ht.new_locus))
                ht = ht.key_by(locus=ht.new_locus, alleles=ht.alleles)

                vep = hl.read_table(var_annotations_ht_path('vep', *TRANCHE_DATA['300k']))
                vep = process_consequences(vep)
                vep = vep.explode(vep.vep.worst_csq_by_gene_canonical)
                vep = vep.annotate(tx_annotation_csq=ht[vep.key].tx_annotation.csq,
                                   mean_proportion=ht[vep.key].tx_annotation.mean_proportion)

                var = var.annotate_rows(pathogenicity=clinvar_ht[var.locus, var.alleles]['pathogenicity'], 
                                        annotation=hl.if_else(hl.literal({'missense', 'LC'}).contains(var.annotation), 'missense|LC', var.annotation), 
                                        polyphen2=vep[var.row_key].vep.worst_csq_by_gene_canonical.polyphen_prediction,
                                        amino_acids=vep[var.row_key].vep.worst_csq_by_gene_canonical.amino_acids,
                                        lof=vep[var.row_key].vep.worst_csq_by_gene_canonical.lof,
                                        most_severe_consequence=vep[var.row_key].vep.worst_csq_by_gene_canonical.most_severe_consequence,
                                        tx_annotation_csq=vep[var.row_key].tx_annotation_csq,
                                        mean_proportion_expressed=vep[var.row_key].mean_proportion)
            var_sig = var.rows()
            var_sig.write(get_ukb_exomes_sumstat_path(subdir=sig_folder, dataset=f'sig_cnt{filter_flag}', result_type='var', extension=args.extension), overwrite=args.overwrite)

            pheno_var_sig = var.cols()
            pheno_var_sig.write(get_ukb_exomes_sumstat_path(subdir=sig_folder, dataset=f'pheno_sig_cnt{filter_flag}', result_type='var', extension=args.extension), overwrite=args.overwrite)

    if args.get_sig_betas:
        if args.variant_results:
            var = hl.read_matrix_table(get_results_mt_path('variant'))
            var = var.annotate_entries(BETA_sig=hl.or_missing(var.Pvalue < 1e-6, var.BETA))
            var = var.select_cols('n_cases')
            var = var.annotate_cols(rare_beta=hl.agg.filter(var.AF < 0.001, hl.agg.mean(var.BETA_sig)), 
                                    common_beta=hl.agg.filter(var.AF >= 0.001, hl.agg.mean(var.BETA_sig)), )
            var.cols().write(get_ukb_exomes_sumstat_path(subdir=beta_folder, dataset='var_beta', result_type='var', extension=args.extension), overwrite=args.overwrite)

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
    parser.add_argument('--get_sig_betas', help='Add extra annotation information on variants', action='store_true')
    parser.add_argument('--gene_results', help='Use gene-level results', action='store_true')
    parser.add_argument('--variant_results', help='Use single-variant results', action='store_true')
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






