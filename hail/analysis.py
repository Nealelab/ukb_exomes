import hail as hl
import argparse
import subprocess
import pickle
import logging
from gnomad.utils.vep import process_consequences
from gnomad.resources.grch38.reference_data import clinvar
from gnomad.resources.resource_utils import DataException
from gnomad.utils.vep import vep_or_lookup_vep
from ukbb_common import *
from ukb_exomes import *

logging.basicConfig(format="%(levelname)s (%(name)s %(lineno)s): %(message)s")
logger = logging.getLogger(f"UKB {CURRENT_TRANCHE} exomes")
logger.setLevel(logging.INFO)


def main(args):
    hl._set_flags(no_whole_stage_codegen="1")
    filter_flag = "_filtered" if args.filters else ""
    if args.update_main_tables:
        gene_mt = modify_phenos_mt("gene")
        logger.info(f"Updated main gene result table: {gene_mt.count()}")
        gene_mt.write(get_results_mt_path("gene"), overwrite=args.overwrite)

        gene_mt.cols().write(
            get_results_mt_path("phenotype", extension="ht"), overwrite=args.overwrite
        )

        var_mt = modify_phenos_mt("variant")
        logger.info(f"Updated main variant result table: {var_mt.count()}")
        var_mt.write(get_results_mt_path("variant"), overwrite=args.overwrite)

    if args.update_util_tables:
        coverage = compute_mean_coverage_ht()
        logger.info(f"Updated gene coverage table: {coverage.count()}")
        coverage.describe()
        coverage.write(get_util_info_path("coverage"), overwrite=args.overwrite)
        caf = get_caf_info_ht()
        logger.info(f"Updated gene CAF table: {coverage.count()}")
        coverage.describe()
        caf.write(get_util_info_path("caf"), overwrite=args.overwrite)

    if args.update_lambda_tables:
        write_lambda_hts(
            result_type="gene",
            freq_min=args.freq_min,
            expected_AC_min=args.expected_AC_min,
            n_var_min=args.n_var_min,
            coverage_min=args.coverage_min,
            var_filter=args.filters,
            extension="ht",
            overwrite=args.overwrite,
        )
        compute_lambda_gc_ht(result_type="gene").write(
            get_ukb_exomes_sumstat_path(
                subdir="qc/lambda_gc", dataset="lambda_by_pheno_full"
            ),
            overwrite=args.overwrite,
        )
        compute_lambda_gc_ht(result_type="gene", by_gene=True).write(
            get_ukb_exomes_sumstat_path(
                subdir="qc/lambda_gc",
                dataset=f"lambda_by_gene",
                result_type="",
            ),
            overwrite=overwrite,
        )

    if args.update_pheno_corr_tables:
        phenos_to_remove = get_corr_phenos_ht(
            r_2=0.5, tie_breaker=more_cases_tie_breaker
        )
        phenos_to_remove.write(
            get_ukb_exomes_sumstat_path(
                subdir="qc", dataset="correlated", result_type="phenos"
            ),
            overwrite=args.overwrite,
        )

        pheno_mt = get_ukb_pheno_mt()
        pheno_mt = drop_pheno_fields_mt(pheno_mt)
        ht = hl.read_matrix_table(get_results_mt_path("variant")).cols()
        pheno_mt = pheno_mt.filter_cols(hl.is_defined(ht[pheno_mt.col_key]))
        corr = make_pairwise_ht(
            pheno_mt, pheno_field=pheno_mt.both_sexes, correlation=True
        )
        corr2 = corr.annotate(
            i_pheno=corr.i_data.phenocode,
            j_pheno=corr.j_data.phenocode,
            i_coding=corr.i_data.coding,
            j_coding=corr.j_data.coding,
        )
        corr2 = corr2.drop("i_data", "j_data")
        corr2.write(
            get_ukb_exomes_sumstat_path(
                subdir="qc", dataset="correlation_table", result_type="phenos"
            ),
            overwrite=args.overwrite,
        )
        corr.entry.export(
            get_ukb_exomes_sumstat_path(
                subdir="qc",
                dataset="pheno_correlation_before_filter",
                result_type="",
                extension="txt.bgz",
            )
        )

        if args.get_related_pheno_cnts:
            phenos_ht = hl.read_table(get_results_mt_path("pheno", extension="ht"))
            lst = get_related_pheno_cnt_list(phenos_ht)
            local = subprocess.check_output("pwd", shell=True)
            local = local.decode("utf-8").split("\n")[0]
            pickle.dump(lst, open("related_pheno_cnt_list", "wb"))
            hl.hadoop_copy(
                f"file://{local}/related_pheno_cnt_list",
                f"gs://ukbb-exome-public/{CURRENT_TRANCHE}/qc/pheno_correlated_cnt_list_{CURRENT_TRANCHE}",
            )

    if args.update_qc_tables:
        gene = annotate_qc_metric_mt(
            result_type="gene",
            lambda_min=args.lambda_min,
            expected_AC_min=args.expected_AC_min,
            coverage_min=args.coverage_min,
            n_var_min=args.n_var_min,
        )
        var = annotate_qc_metric_mt(
            "variant", expected_AC_min=args.expected_AC_min, lambda_min=args.lambda_min
        )

        gene.write(
            get_ukb_exomes_sumstat_path(
                subdir="qc",
                dataset="gene_qc_metrics_ukb_exomes",
                result_type="",
                extension="mt",
            ),
            overwrite=args.overwrite,
        )
        gene.rows().write(
            get_ukb_exomes_sumstat_path(
                subdir="qc",
                dataset="gene_qc_metrics_ukb_exomes",
                result_type="",
                extension="ht",
            ),
            overwrite=args.overwrite,
        )
        gene.cols().write(
            get_ukb_exomes_sumstat_path(
                subdir="qc",
                dataset="pheno_qc_metrics_ukb_exomes",
                result_type="",
                extension="ht",
            ),
            overwrite=args.overwrite,
        )
        var.write(
            get_ukb_exomes_sumstat_path(
                subdir="qc",
                dataset="variant_qc_metrics_ukb_exomes",
                result_type="",
                extension="mt",
            ),
            overwrite=args.overwrite,
        )
        var_ht = var.rows()
        var_ht = var_ht.annotate(**var_ht.call_stats)
        var_ht.write(
            get_ukb_exomes_sumstat_path(
                subdir="qc",
                dataset="variant_qc_metrics_ukb_exomes",
                result_type="",
                extension="ht",
            ),
            overwrite=args.overwrite,
        )

    if args.update_sig_cnt_tables:
        gene = get_sig_cnt_mt(
            result_type="gene", test_type=args.test_type, filters=args.filters
        )
        gene_sig = gene.rows()
        gene_sig.write(
            get_ukb_exomes_sumstat_path(
                subdir="analysis",
                dataset=f"gene_sig_cnt{filter_flag}_{args.test_type}",
                result_type="",
            ),
            overwrite=args.overwrite,
        )
        pheno_sig = gene.cols()
        pheno_sig.write(
            get_ukb_exomes_sumstat_path(
                subdir="analysis",
                dataset=f"pheno_sig_cnt{filter_flag}_{args.test_type}",
                result_type="gene",
            ),
            overwrite=args.overwrite,
        )

        var = get_sig_cnt_mt(
            result_type="variant", test_type=args.test_type, filters=args.filters
        )
        var = annotate_additional_info_mt(var, "variant")
        var_sig = var.rows()
        var_sig = var_sig.annotate(**var_sig.call_stats)
        var_sig.write(
            get_ukb_exomes_sumstat_path(
                subdir="analysis",
                dataset=f"var_sig_cnt{filter_flag}_{args.test_type}",
                result_type="",
            ),
            overwrite=args.overwrite,
        )
        pheno_var_sig = var.cols()
        pheno_var_sig.write(
            get_ukb_exomes_sumstat_path(
                subdir="analysis",
                dataset=f"pheno_sig_cnt{filter_flag}_{args.test_type}",
                result_type="var",
            ),
            overwrite=args.overwrite,
        )

        var_gene = compare_gene_var_sig_cnt_mt(
            test_type=args.test_type, filters=args.filters
        )
        var_gene.cols().write(
            get_ukb_exomes_sumstat_path(
                subdir="analysis",
                dataset=f"var_gene_comparison_by_pheno{filter_flag}_{args.test_type}_3annt",
                result_type="",
            ),
            overwrite=args.overwrite,
        )

        var_gene.rows().write(
            get_ukb_exomes_sumstat_path(
                subdir="analysis",
                dataset=f"var_gene_comparison_by_gene{filter_flag}_{args.test_type}_3annt",
                result_type="",
            ),
            overwrite=args.overwrite,
        )

    if args.update_icd_tables:
        icd_var = get_icd_min_p_ht(
            result_type="variant", test_type=args.test_type, filters=args.filters
        )
        icd_gene = get_icd_min_p_ht(
            result_type="gene", test_type=args.test_type, filters=args.filters
        )
        icd_gene.describe()
        icd_var.write(
            get_ukb_exomes_sumstat_path(
                subdir="analysis",
                dataset=f"icd_min_p_var{filter_flag}_{args.test_type}",
                result_type="",
            ),
            overwrite=args.overwrite,
        )
        icd_gene.write(
            get_ukb_exomes_sumstat_path(
                subdir="analysis",
                dataset=f"icd_min_p_gene{filter_flag}_{args.test_type}",
                result_type="",
            ),
            overwrite=args.overwrite,
        )

    if args.update_random_pheno_tables:
        rp_lambda_caf = compute_lambdas_by_freq_interval_ht(
            result_type="gene", random_phenos=True
        )
        rp_lambda_caf.write(
            get_ukb_exomes_sumstat_path(
                subdir="analysis", dataset=f"rp_lambda_gene_caf", result_type=""
            ),
            overwrite=args.overwrite,
        )
        rp_gene_p_skato = get_pvalue_by_freq_interval_ht(
            result_type="gene", test_type="skato", random_phenos=True
        )
        rp_gene_p_skato.write(
            get_ukb_exomes_sumstat_path(
                subdir="analysis", dataset=f"rp_gene_p_value_skato", result_type=""
            ),
            overwrite=args.overwrite,
        )
        rp_gene_p_burden = get_pvalue_by_freq_interval_ht(
            result_type="gene", test_type="burden", random_phenos=True
        )
        rp_gene_p_burden.write(
            get_ukb_exomes_sumstat_path(
                subdir="analysis", dataset=f"rp_gene_p_value_burden", result_type=""
            ),
            overwrite=args.overwrite,
        )
        rp_var_p = get_pvalue_by_freq_interval_ht(
            result_type="variant", random_phenos=True
        )
        rp_var_p.write(
            get_ukb_exomes_sumstat_path(
                subdir="analysis", dataset=f"rp_var_p_value", result_type=""
            ),
            overwrite=args.overwrite,
        )

    if args.convert_to_txt_bgz:
        export_all_ht_to_txt_bgz("qc/lambda_gc")
        export_all_ht_to_txt_bgz("qc")
        export_all_ht_to_txt_bgz("analysis")


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--overwrite", help="Overwrite everything", action="store_true")
    parser.add_argument(
        "--update_main_tables",
        help="Update the original result tables",
        action="store_true",
    )
    parser.add_argument(
        "--update_util_tables", help="Update the utility tables", action="store_true"
    )
    parser.add_argument(
        "--update_lambda_tables",
        help="Update the lambda gc tables",
        action="store_true",
    )
    parser.add_argument(
        "--update_pheno_corr_tables",
        help="Update the phenotype correlation tables",
        action="store_true",
    )
    parser.add_argument(
        "--update_qc_tables",
        help="Update the tables with QC metrics",
        action="store_true",
    )
    parser.add_argument(
        "--update_sig_cnt_tables",
        help="Update significant association tables",
        action="store_true",
    )
    parser.add_argument(
        "--update_icd_tables",
        help="Update minimum pvalue tables for icd phenotypes",
        action="store_true",
    )
    parser.add_argument(
        "--update_random_pheno_tables",
        help="Update tables for random phenotypes",
        action="store_true",
    )
    parser.add_argument(
        "--convert_to_txt_bgz",
        help="Convert all hail Tables to .txt.bgz file",
        action="store_true",
    )
    parser.add_argument(
        "--get_related_pheno_cnts",
        help="Count the number of correlated phenotypes to remove",
        action="store_true",
    )
    parser.add_argument(
        "--coverage_min", type=int, help="Keep genes with higher coverage"
    )
    parser.add_argument(
        "--n_var_min", type=int, help="Keep genes with larger number of variants"
    )
    parser.add_argument(
        "--freq_min",
        type=float,
        help="Keep genes/variants with higher cumulative allele frequency",
    )
    parser.add_argument(
        "--expected_AC_min",
        type=int,
        help="Keep genes/variants with higher expected AC",
    )
    parser.add_argument(
        "--lambda_min",
        type=float,
        help="Remove genes/phenotypes with lower lambda value",
    )
    parser.add_argument(
        "--test_type",
        type=str,
        help="Test results to apply lambda filters on: skato OR burden",
    )
    parser.add_argument("--filters", help="Apply filters", action="store_true")
    args = parser.parse_args()
    print(args)

    main(args)
