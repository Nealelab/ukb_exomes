from ukbb_qc.resources.sample_qc import interval_qc_path
from ukbb_qc.resources.variant_qc import var_annotations_ht_path
from ukbb_qc.resources.basics import release_ht_path
from ukbb_common import *
from gnomad.resources.grch38.reference_data import clinvar
from gnomad.utils.vep import process_consequences
from ..resources import *


ANNOTATIONS = ("pLoF", "missense|LC", "synonymous", "pLoF|missense|LC")
TESTS = ("skato", "skat", "burden")
TRAIT_TYPES = ("continuous", "categorical", "icd10")
P_VALUE_FIELDS = {"skato": "Pvalue", "skat": "Pvalue_SKAT", "burden": "Pvalue_Burden"}
EMPIRICAL_P_THRESHOLDS = {"skato": 2.5e-7, "burden": 6.7e-7, "variant": 8e-9}


def get_ukb_exomes_sumstat_path(
    subdir: str,
    dataset: str,
    result_type: str = "gene",
    extension: str = "ht",
    tranche: str = CURRENT_TRANCHE,
):
    """
    Get UKB downstream sumstats Table path
    :param str subdir: Path to sub directory
    :param str dataset: Type of association sumstats
    :param str result_type: Type of association test result, `gene` or `variant`
    :param str extension: Extension of output file
    :param str tranche: Tranche of data to use
    :return: Path to sumstats Table
    :rtype: str
    """
    result_type = result_type if result_type == "" else f"_{result_type}"
    return f"{public_bucket}/{tranche}/{subdir}/{dataset}{result_type}_{tranche}.{extension}"


def get_util_info_path(type: str, tranche: str = CURRENT_TRANCHE):
    """
    Get cumulative allele frequency or mean coverage Table path
    :param str type: Type of table to use, `coverage` or `other`
    :param str tranche: Tranche of data to use
    :return: Path to information Table
    :rtype: str
    """
    type = "_mean_coverage" if type == "coverage" else "_caf"
    return f"{public_bucket}/{tranche}/qc/gene{type}_{tranche}.ht"


def compute_lambda_gc_ht(
    result_type: str = "gene",
    by_annotation: bool = False,
    by_gene: bool = False,
    expected_AC_min: int = None,
    freq_min: float = None,
    freq_max: float = None,
    n_var_min: int = None,
    coverage_min: int = None,
    var_filter: bool = False,
    random_phenos: bool = False,
    tranche: str = CURRENT_TRANCHE,
):
    """
    Compute Lambda GC for variant-level test or gene-level tests (SKAT-O, SKAT, and Burden Test) under specified condition
    :param str result_type: Type of association test result, `gene` or `variant`
    :param bool by_annotation: Whether to group by annotation categories
    :param bool by_gene: Whether to compute lambda GC for each gene, applicable only for result_type == `gene`
    :param int expected_AC_min: the lower bound used to remove genes/variants with 0 phenotypes having expected AC greater than this.
    :param float freq_min: Lower bound of allele frequency/cumulative allele frequency used for variant/gene filtering
    :param float freq_max: Upper bound of allele frequency/cumulative allele frequency used for variant/gene filtering
    :param int n_var_min: Minimum number of variants included in a gene-based test used for gene filtering
    :param int coverage_min: Minimum coverage used for gene filtering
    :param bool var_filter: Whether to filter out test statistics with SE=0
    :param bool random_pheno: Whether to use randomly-generated phenotypes
    :param str tranche: Tranche of data to use
    :return: Hail table with lambda GC
    :rtype: Table
    """
    mt = hl.read_matrix_table(
        get_results_mt_path(result_type, tranche=tranche, random_phenos=random_phenos)
    )
    if result_type == "gene":  # Gene Info
        mt = annotate_additional_info_mt(
            mt=mt, result_type=result_type, tranche=tranche
        )
        if expected_AC_min is not None:
            logger.info(
                f"Removing {result_type}s with 0 phenotypes having expected AC >= %d...",
                expected_AC_min,
            )
            mt = mt.annotate_rows(
                expected_ac_row_filter=hl.agg.count_where(
                    mt.expected_AC >= expected_AC_min
                )
            )
            mt = mt.filter_entries(mt.expected_AC >= 50)
            # mt = mt.filter_rows(mt.expected_ac_row_filter > 0)
        if freq_min is not None:
            logger.info(f"Removing {result_type}s with CAF <= {freq_min}...")
            mt = mt.filter_rows(mt.CAF > freq_min)
        if freq_max is not None:
            logger.info(f"Removing {result_type}s with CAF > {freq_max}...")
            mt = mt.filter_rows(mt.CAF <= freq_max)
        if n_var_min is not None:
            logger.info(
                f"Removing {result_type}s with number of variants < {n_var_min}..."
            )
            mt = mt.filter_rows(mt.total_variants >= n_var_min)
        if coverage_min is not None:
            logger.info(
                f"Removing {result_type}s with mean coverage < {coverage_min}..."
            )
            mt = mt.filter_rows(mt.mean_coverage >= coverage_min)
        if by_annotation:
            logger.info(
                f"Computing lambda GC for {result_type}-based test results by annotation groups..."
            )
            mt = mt.select_cols(
                "n_cases",
                lambda_gc_skato=hl.agg.group_by(
                    mt.annotation, hl.methods.statgen._lambda_gc_agg(mt.Pvalue)
                ),
                lambda_gc_skat=hl.agg.group_by(
                    mt.annotation, hl.methods.statgen._lambda_gc_agg(mt.Pvalue_SKAT)
                ),
                lambda_gc_burden=hl.agg.group_by(
                    mt.annotation, hl.methods.statgen._lambda_gc_agg(mt.Pvalue_Burden)
                ),
                CAF_range=f"CAF:({freq_min}, {freq_max}]",
            )
            mt = mt.annotate_cols(
                **{
                    f"{annotation}_lambda_gc_{test}": mt[f"lambda_gc_{test}"].get(
                        annotation
                    )
                    for annotation in ANNOTATIONS
                    for test in TESTS
                }
            )
            ht = mt.cols()
            ht = ht.annotate(
                trait_type2=hl.if_else(
                    ht.trait_type == "icd_first_occurrence", "icd10", ht.trait_type
                ),
            )
        elif by_gene:
            logger.info(f"Computing gene-level lambda GC")
            mt = mt.annotate_cols(
                trait_type2=hl.if_else(
                    mt.trait_type == "continuous", mt.trait_type, "categorical"
                )
            )
            mt = mt.select_rows(
                "CAF",
                "mean_coverage",
                lambda_gc_skato=hl.agg.group_by(
                    mt.trait_type2, hl.methods.statgen._lambda_gc_agg(mt.Pvalue)
                ),
                lambda_gc_skat=hl.agg.group_by(
                    mt.trait_type2, hl.methods.statgen._lambda_gc_agg(mt.Pvalue_SKAT)
                ),
                lambda_gc_burden=hl.agg.group_by(
                    mt.trait_type2, hl.methods.statgen._lambda_gc_agg(mt.Pvalue_Burden)
                ),
            )
            mt = mt.annotate_rows(
                **{
                    f"{trait_type}_lambda_gc_{test}": mt[f"lambda_gc_{test}"].get(
                        trait_type
                    )
                    for trait_type in TRAIT_TYPES
                    for test in TESTS
                },
                all_lambda_gc_skato=hl.agg.filter(
                    hl.is_defined(mt.Pvalue),
                    hl.methods.statgen._lambda_gc_agg(mt.Pvalue),
                ),
                all_lambda_gc_skat=hl.agg.filter(
                    hl.is_defined(mt.Pvalue_SKAT),
                    hl.methods.statgen._lambda_gc_agg(mt.Pvalue_SKAT),
                ),
                all_lambda_gc_burden=hl.agg.filter(
                    hl.is_defined(mt.Pvalue_Burden),
                    hl.methods.statgen._lambda_gc_agg(mt.Pvalue_Burden),
                ),
            )
            ht = mt.rows()
        else:
            logger.info(
                f"Computing lambda GC for {result_type}-based test results without stratification..."
            )
            mt = mt.select_cols(
                "n_cases",
                lambda_gc_skato=hl.agg.filter(
                    hl.is_defined(mt.Pvalue),
                    hl.methods.statgen._lambda_gc_agg(mt.Pvalue),
                ),
                lambda_gc_skat=hl.agg.filter(
                    hl.is_defined(mt.Pvalue_SKAT),
                    hl.methods.statgen._lambda_gc_agg(mt.Pvalue_SKAT),
                ),
                lambda_gc_burden=hl.agg.filter(
                    hl.is_defined(mt.Pvalue_Burden),
                    hl.methods.statgen._lambda_gc_agg(mt.Pvalue_Burden),
                ),
                CAF_range=f"CAF:({freq_min}, {freq_max}]",
            )
            ht = mt.cols()
            ht = ht.annotate(
                trait_type2=hl.if_else(
                    ht.trait_type == "icd_first_occurrence", "icd10", ht.trait_type
                ),
            )
    else:  # Variant Info
        if freq_min is not None:
            mt = mt.filter_rows(mt.call_stats.AF > freq_min)
        if freq_max is not None:
            mt = mt.filter_rows(mt.call_stats.AF <= freq_max)
        if expected_AC_min is not None:
            mt = annotate_additional_info_mt(mt, "variant", True)
            mt = mt.annotate_rows(
                expected_ac_row_filter=hl.agg.count_where(
                    mt.expected_AC >= expected_AC_min
                )
            )
            mt = mt.filter_rows(mt.expected_ac_row_filter > 0)
        mt = mt.annotate_rows(
            annotation=hl.if_else(
                hl.literal({"missense", "LC"}).contains(mt.annotation),
                "missense|LC",
                mt.annotation,
            )
        )
        if var_filter:
            mt = mt.filter_rows(hl.is_defined(mt.annotation))
            if by_annotation:
                lambda_gc = hl.agg.group_by(
                    mt.annotation,
                    hl.agg.filter(
                        mt.SE != 0, hl.methods.statgen._lambda_gc_agg(mt.Pvalue)
                    ),
                )
                mt = mt.select_cols(
                    "n_cases",
                    lambda_gc=lambda_gc,
                    AF_range=f"AF:({freq_min}, {freq_max}]",
                )
                mt = mt.annotate_cols(
                    **{
                        f"{annotation}_lambda_gc_variant": mt.lambda_gc.get(annotation)
                        for annotation in ANNOTATIONS
                    },
                )
            else:
                lambda_gc = hl.agg.filter(
                    mt.SE != 0, hl.methods.statgen._lambda_gc_agg(mt.Pvalue)
                )
                mt = mt.select_cols(
                    "n_cases",
                    lambda_gc=lambda_gc,
                    AF_range=f"AF:({freq_min}, {freq_max}]",
                )
        else:
            if by_annotation:
                lambda_gc = hl.agg.group_by(
                    mt.annotation, hl.methods.statgen._lambda_gc_agg(mt.Pvalue)
                )
                mt = mt.select_cols(
                    "n_cases",
                    lambda_gc=lambda_gc,
                    AF_range=f"AF:({freq_min}, {freq_max}]",
                )
                mt = mt.annotate_cols(
                    **{
                        f"{annotation}_lambda_gc_variant": mt.lambda_gc[annotation]
                        for annotation in ANNOTATIONS
                    },
                )
            else:
                lambda_gc = hl.agg.filter(
                    hl.is_defined(mt.Pvalue),
                    hl.methods.statgen._lambda_gc_agg(mt.Pvalue),
                )
                mt = mt.select_cols(
                    "n_cases",
                    lambda_gc=lambda_gc,
                    AF_range=f"AF:({freq_min}, {freq_max}]",
                )
        ht = mt.cols()
        ht = ht.annotate(
            trait_type2=hl.if_else(
                ht.trait_type == "icd_first_occurrence", "icd10", ht.trait_type
            ),
        )
    return ht


def compute_lambdas_by_freq_interval_ht(
    result_type: str = "gene",
    by_annotation: bool = False,
    freq_breaks: list = [0.0001, 0.001, 0.01, 0.1],
    n_var_min: int = None,
    coverage_min: int = None,
    expected_AC_min: int = None,
    var_filter: bool = False,
    random_phenos: bool = False,
    tranche: str = CURRENT_TRANCHE,
):
    """
    Compute Lambda GC for variant-level test or gene-level tests (SKAT-O, SKAT, and Burden Test) by frequency intervals
    :param str result_type: Type of association test result, `gene` or `variant`
    :param bool by_annotation: Whether to group by annotation categories
    :param list freq_breaks: Breaks to divide frequency into bins
    :param int n_var_min: Minimum number of variants included in a gene-based test used for gene filtering
    :param int coverage_min: Minimum coverage used for gene filtering
    :param int expected_AC_min: the lower bound used to remove genes/variants with 0 phenotypes having expected AC greater than this.
    :param bool var_filter: Whether to filter out test statistics with SE=0
    :param bool random_pheno: Whether to use randomly-generated phenotypes
    :param str tranche: Tranche of data to use
    :return: Hail table with lambda GC by frequency intervals
    :rtype: Table
    """
    logger.info(
        f"Computing lambda GC for {result_type}-based test results across frequency intervals"
    )
    ht = compute_lambda_gc_ht(
        result_type=result_type,
        by_annotation=by_annotation,
        freq_max=freq_breaks[0],
        n_var_min=n_var_min,
        coverage_min=coverage_min,
        expected_AC_min=expected_AC_min,
        random_phenos=random_phenos,
        tranche=tranche,
    )
    for i in list(range(0, len(freq_breaks) - 1)):
        sub_ht = compute_lambda_gc_ht(
            result_type=result_type,
            by_annotation=by_annotation,
            freq_min=freq_breaks[i],
            freq_max=freq_breaks[i + 1],
            n_var_min=n_var_min,
            coverage_min=coverage_min,
            expected_AC_min=expected_AC_min,
            var_filter=var_filter,
            random_phenos=random_phenos,
            tranche=tranche,
        )
        ht = ht.union(sub_ht)
    ht = ht.union(
        compute_lambda_gc_ht(
            result_type=result_type,
            by_annotation=by_annotation,
            freq_min=freq_breaks[-1],
            n_var_min=n_var_min,
            coverage_min=coverage_min,
            expected_AC_min=expected_AC_min,
            var_filter=var_filter,
            random_phenos=random_phenos,
            tranche=tranche,
        )
    )
    return ht


def compute_lambdas_by_expected_ac_ht(
    result_type: str = "gene",
    ac_breaks: list = [0, 5, 50, 500, 5000, 50000],
    n_var_min: int = None,
    coverage_min: int = None,
    var_filter: bool = False,
    random_phenos: bool = False,
    tranche: str = CURRENT_TRANCHE,
):
    """
    Compute Lambda GC for variant-level test by expected allele count interval
    :param str result_type: Type of association test result, `gene` or `variant`
    :param list ac_breaks: Breaks to divide expected allele count into bins
    :param bool random_pheno: Whether to use randomly-generated phenotypes
    :param str tranche: Tranche of data to use
    :return: Hail table with lambda GC by expected allele count intervals
    :rtype: Table
    """
    logger.info(
        f"Computing lambda GC for {result_type}-based test results across expected AC bins"
    )
    mt = hl.read_matrix_table(
        get_results_mt_path(result_type, tranche=tranche, random_phenos=random_phenos)
    )
    mt = annotate_additional_info_mt(
        mt=mt, result_type=result_type, expected_AC_only=True, tranche=tranche
    )
    if result_type == "gene":
        if n_var_min is not None:
            logger.info(
                f"Removing {result_type}s with number of variants < {n_var_min}..."
            )
            mt = mt.filter_rows(mt.total_variants >= n_var_min)
        if coverage_min is not None:
            logger.info(
                f"Removing {result_type}s with mean coverage < {coverage_min}..."
            )
            mt = mt.filter_rows(mt.mean_coverage >= coverage_min)
        ht = mt.annotate_cols(
            expected_AC_range=f"expected_AC:({ac_breaks[-1]}, )",
            **{
                f"lambda_gc_{test_type}": hl.agg.filter(
                    hl.is_defined(mt[P_VALUE_FIELDS[test_type]])
                    & (mt.expected_AC > ac_breaks[-1]),
                    hl.methods.statgen._lambda_gc_agg(mt[P_VALUE_FIELDS[test_type]]),
                )
                for test_type in TESTS
            },
            **{
                f"n_{result_type}s_{test_type}": hl.agg.count_where(
                    hl.is_defined(mt[P_VALUE_FIELDS[test_type]])
                    & (mt.expected_AC > ac_breaks[-1])
                )
                for test_type in TESTS
            },
        ).cols()
        for i in list(range(0, len(ac_breaks) - 1)):
            ht_tmp = mt.annotate_cols(
                expected_AC_range=f"expected_AC:({ac_breaks[i]}, {ac_breaks[i + 1]}]",
                **{
                    f"lambda_gc_{test_type}": hl.agg.filter(
                        hl.is_defined(mt[P_VALUE_FIELDS[test_type]])
                        & (mt.expected_AC > ac_breaks[i])
                        & (mt.expected_AC <= ac_breaks[i + 1]),
                        hl.methods.statgen._lambda_gc_agg(
                            mt[P_VALUE_FIELDS[test_type]]
                        ),
                    )
                    for test_type in TESTS
                },
                **{
                    f"n_{result_type}s_{test_type}": hl.agg.count_where(
                        hl.is_defined(mt[P_VALUE_FIELDS[test_type]])
                        & (mt.expected_AC > ac_breaks[i])
                        & (mt.expected_AC <= ac_breaks[i + 1])
                    )
                    for test_type in TESTS
                },
            ).cols()
            ht = ht.union(ht_tmp)
        ht = ht.select(
            *{f"lambda_gc_{test_type}" for test_type in TESTS},
            *{f"n_genes_{test_type}" for test_type in TESTS},
            "expected_AC_range",
            "n_cases_defined",
        )
    else:  # variant
        if random_phenos:
            mt = (
                mt.annotate_rows(AF2=hl.agg.take(mt.AF, 1)[0])
                .drop("AF")
                .rename({"AF2": "AF"})
            )
        mt = mt.filter_rows(hl.is_defined(mt.annotation))
        ht = mt.annotate_cols(
            expected_AC_range=f"expected_AC:({ac_breaks[-1]}, )",
            lambda_gc=hl.agg.filter(
                hl.is_defined(mt.Pvalue)
                & (mt.expected_AC > ac_breaks[-1])
                & (mt.SE != 0),
                hl.methods.statgen._lambda_gc_agg(mt.Pvalue),
            ),
            n_vars=hl.agg.count_where(
                hl.is_defined(mt.Pvalue) & (mt.expected_AC > 100000) & (mt.SE != 0)
            ),
        ).cols()
        for i in list(range(0, len(ac_breaks) - 1)):
            ht_tmp = mt.annotate_cols(
                expected_AC_range=f"expected_AC:({ac_breaks[i]}, {ac_breaks[i + 1]}]",
                lambda_gc=hl.agg.filter(
                    hl.is_defined(mt.Pvalue)
                    & (mt.expected_AC > ac_breaks[i])
                    & (mt.expected_AC <= ac_breaks[i + 1])
                    & (mt.SE != 0),
                    hl.methods.statgen._lambda_gc_agg(mt.Pvalue),
                ),
                n_vars=hl.agg.count_where(
                    hl.is_defined(mt.Pvalue)
                    & (mt.expected_AC > ac_breaks[i])
                    & (mt.expected_AC <= ac_breaks[i + 1])
                    & (mt.SE != 0)
                ),
            ).cols()
            ht = ht.union(ht_tmp)
        ht = ht.select("lambda_gc", "n_vars", "expected_AC_range", "n_cases_defined")
    ht = ht.annotate(
        trait_type2=hl.if_else(
            ht.trait_type == "icd_first_occurrence", "icd10", ht.trait_type
        ),
    )
    return ht


def write_lambda_hts(
    result_type="gene",
    freq_min: float = None,
    expected_AC_min: int = None,
    n_var_min: int = None,
    coverage_min: int = None,
    var_filter: bool = False,
    random_phenos: bool = False,
    extension: str = "ht",
    overwrite: bool = False,
    tranche: str = CURRENT_TRANCHE,
):
    """
    Write lambda GC tables for variant-level test or gene-level tests (SKAT-O, SKAT, and Burden Test) under different condition
    :param str result_type: Type of association test result, `gene` or `variant`
    :param float freq_min: Lower bound of allele frequency/cumulative allele frequency used for variant/gene filtering
    :param int expected_AC_min: the lower bound used to remove genes/variants with 0 phenotypes having expected AC greater than this.
    :param int n_var_min: Minimum number of variants included in a gene-based test used for gene filtering
    :param int coverage_min: Minimum coverage used for gene filtering
    :param bool var_filter: Whether to filter out test statistics with SE=0
    :param bool random_pheno: Whether to use randomly-generated phenotypes
    :param str extension: Extension of output file
    :param bool overwrite: Whether to overwrite exsiting file
    :param str tranche: Tranche of data to use
    """
    kwargs = {
        "result_type": result_type,
        "n_var_min": n_var_min,
        "coverage_min": coverage_min,
        "var_filter": var_filter,
        "random_phenos": random_phenos,
    }
    pathargs = {
        "subdir": "qc/lambda_gc",
        "result_type": result_type,
        "tranche": tranche,
        "extension": extension,
    }
    rp = "rp_" if random_phenos else ""
    filter = (
        "_filtered"
        if any(
            v is not None for v in [freq_min, n_var_min, coverage_min, expected_AC_min]
        )
        else ""
    )
    freq_breaks = (
        [freq_min, 0.001, 0.01, 0.1]
        if freq_min is not None
        else [0.0001, 0.001, 0.01, 0.1]
    )

    compute_lambda_gc_ht(
        freq_min=freq_min, expected_AC_min=expected_AC_min, **kwargs
    ).write(
        get_ukb_exomes_sumstat_path(
            dataset=f"{rp}lambda_by_pheno_full{filter}", **pathargs
        ),
        overwrite=overwrite,
    )
    compute_lambda_gc_ht(
        by_annotation=True, freq_min=freq_min, expected_AC_min=expected_AC_min, **kwargs
    ).write(
        get_ukb_exomes_sumstat_path(
            dataset=f"{rp}lambda_by_pheno_annt{filter}", **pathargs
        ),
        overwrite=overwrite,
    )
    compute_lambdas_by_freq_interval_ht(
        freq_breaks=freq_breaks, expected_AC_min=expected_AC_min, **kwargs
    ).write(
        get_ukb_exomes_sumstat_path(
            dataset=f"{rp}lambda_by_pheno_freq{filter}", **pathargs
        ),
        overwrite=overwrite,
    )

    compute_lambdas_by_expected_ac_ht(
        ac_breaks=[0, 5, 50, 500, 5000, 50000], **kwargs
    ).write(
        get_ukb_exomes_sumstat_path(
            dataset=f"{rp}lambda_by_pheno_expected_ac{filter}", **pathargs
        ),
        overwrite=overwrite,
    )

    if result_type == "gene":
        compute_lambda_gc_ht(
            by_gene=True, expected_AC_min=expected_AC_min, **kwargs
        ).write(
            get_ukb_exomes_sumstat_path(
                subdir="qc/lambda_gc",
                dataset=f"{rp}lambda_by_gene{filter}",
                result_type="",
                extension=extension,
                tranche=tranche,
            ),
            overwrite=overwrite,
        )
        compute_lambda_gc_ht(
            result_type=result_type,
            by_gene=True,
            freq_min=freq_min,
            expected_AC_min=expected_AC_min,
            n_var_min=n_var_min,
        ).write(
            get_ukb_exomes_sumstat_path(
                subdir="qc/lambda_gc",
                dataset=f"{rp}lambda_by_gene_before_coverage{filter}",
                result_type="",
                extension=extension,
                tranche=tranche,
            ),
            overwrite=overwrite,
        )


def compute_ukb_pheno_moments_ht(pheno_sex="both_sexes", phenocode: list = None):
    """
    Generate first 4 moments Table for selected phenotypes
    :param str pheno_sex: Which sex group to use
    :param list phenocode: Subset of phenotypes to use
    :return: Hail Table of moments
    :rtype: Table
    """
    pheno = get_ukb_pheno_mt()
    if phenocode is not None:
        pheno = pheno.filter_cols(hl.literal(phenocode).contains(pheno.phenocode))
    pheno = pheno.annotate_cols(
        mean=hl.agg.mean(pheno[pheno_sex]), std=hl.agg.stats(pheno[pheno_sex]).stdev
    )
    pheno = pheno.annotate_cols(
        variance=pheno.std ** 2,
        skewness=hl.agg.mean((pheno[pheno_sex] - pheno.mean) ** 3) / (pheno.std ** 3),
        kurtosis=hl.agg.mean((pheno[pheno_sex] - pheno.mean) ** 4) / (pheno.std ** 4),
    )
    pheno = pheno.select_cols("mean", "variance", "skewness", "kurtosis").cols()
    return pheno


def get_caf_info_ht(tranche: str = CURRENT_TRANCHE):
    """
    Generate cumulative allele frequency Table
    :param str tranche: Tranche of data to use
    :return: Cumulative allele frequency Hail Table
    :rtype: Table
    """
    vep = hl.read_table(var_annotations_ht_path("vep", *TRANCHE_DATA[tranche]))
    vep = process_consequences(vep)
    vep = vep.explode(vep.vep.worst_csq_by_gene_canonical)
    annotation = annotation_case_builder(vep.vep.worst_csq_by_gene_canonical)
    vep = vep.annotate(
        annotation=hl.if_else(
            hl.literal({"missense", "LC"}).contains(annotation),
            "missense|LC",
            annotation,
        )
    )
    vep = vep.filter(hl.is_defined(vep.annotation) & (vep.annotation != "non-coding"))
    vep = vep.select(vep.vep.worst_csq_by_gene_canonical.gene_id, vep.annotation)

    freq = hl.read_table(release_ht_path(*TRANCHE_DATA[tranche]))
    var_af = vep.annotate(AF=freq[vep.key].freq[0].AF)
    sum_af = var_af.group_by(var_af.annotation, var_af.gene_id).aggregate(
        CAF=hl.agg.sum(var_af.AF)
    )
    if tranche == "500k":
        sub_af = sum_af.filter(
            hl.literal({"missense|LC", "pLoF"}).contains(sum_af.annotation)
        )
        sub_af = sub_af.key_by()
        sub_af = sub_af.annotate(annotation="pLoF|missense|LC")
        sub_af = sub_af.group_by("annotation", "gene_id").aggregate(
            CAF=hl.agg.sum(sub_af.CAF)
        )
        sum_af = sub_af.union(sum_af)
    return sum_af


def get_sig_cnt_mt(
    result_type: str = "gene",
    test_type: str = "skato",
    filters: bool = True,
    tranche: str = CURRENT_TRANCHE,
):
    """
    Generate the number of significant association for each phenotype and each gene/variant
    :param str result_type: Type of association test result, `gene` or `variant`
    :param str test_type: Type of gene-level association test result, `skato` or `burden`
    :param bool filters: Whether to use QCed set
    :param str tranche: Tranche of data to use
    :return: Hail MatrixTable with phenotype-level and gene/variant-level association count
    :rtype: MatrixTable
    """

    if filters:
        mt = get_qc_result_mt(
            result_type=result_type, test_type=test_type, tranche=tranche
        )
    else:
        mt = load_final_sumstats_table(
            result_type=result_type, extension="mt", tranche=tranche
        )
    mt = mt.annotate_cols(
        trait_type2=hl.if_else(
            mt.trait_type == "icd_first_occurrence", "icd10", mt.trait_type
        )
    )
    pvalue = P_VALUE_FIELDS[test_type.lower()]
    if result_type == "gene":
        mt = mt.annotate_rows(
            sig_pheno_cnt=hl.agg.group_by(
                mt.trait_type2,
                hl.agg.count_where(
                    mt[pvalue] < EMPIRICAL_P_THRESHOLDS[test_type.lower()]
                ),
            ),
            all_sig_pheno_cnt=hl.agg.filter(
                hl.is_defined(mt[pvalue]),
                hl.agg.count_where(
                    mt[pvalue] < EMPIRICAL_P_THRESHOLDS[test_type.lower()]
                ),
            ),
        )
        mt = mt.annotate_rows(
            **{
                f"{trait_type2}_sig_pheno_cnt_{test_type.lower()}": mt.sig_pheno_cnt.get(
                    trait_type2
                )
                for trait_type2 in TRAIT_TYPES
            }
        )

        mt = mt.annotate_cols(
            sig_gene_cnt=hl.agg.group_by(
                mt.annotation,
                hl.agg.count_where(
                    mt[pvalue] < EMPIRICAL_P_THRESHOLDS[test_type.lower()]
                ),
            ),
            all_sig_gene_cnt=hl.agg.filter(
                hl.is_defined(mt[pvalue]),
                hl.agg.count_where(
                    mt[pvalue] < EMPIRICAL_P_THRESHOLDS[test_type.lower()]
                ),
            ),
        )
        mt = mt.annotate_cols(
            **{
                f"{annotation}_sig_gene_cnt_{test_type.lower()}": mt.sig_gene_cnt.get(
                    annotation
                )
                for annotation in ANNOTATIONS
            }
        )
    else:
        mt = mt.annotate_rows(
            annotation=hl.if_else(
                hl.literal({"missense", "LC"}).contains(mt.annotation),
                "missense|LC",
                mt.annotation,
            )
        )
        if filters:
            mt = mt.annotate_rows(
                all_sig_pheno_cnt=hl.agg.filter(
                    mt.SE != 0,
                    hl.agg.count_where(mt.Pvalue < EMPIRICAL_P_THRESHOLDS["variant"]),
                ),
                sig_pheno_cnt=hl.agg.group_by(
                    mt.trait_type2,
                    hl.agg.filter(
                        mt.SE != 0,
                        hl.agg.count_where(
                            mt.Pvalue < EMPIRICAL_P_THRESHOLDS["variant"]
                        ),
                    ),
                ),
            )
            mt = mt.annotate_rows(
                **{
                    f"{trait_type}_sig_pheno_cnt": mt.sig_pheno_cnt.get(trait_type)
                    for trait_type in TRAIT_TYPES
                    for test in TESTS
                },
            )
            mt = mt.annotate_cols(
                all_sig_var_cnt=hl.agg.filter(
                    mt.SE != 0,
                    hl.agg.count_where(mt.Pvalue < EMPIRICAL_P_THRESHOLDS["variant"]),
                ),
                sig_var_cnt=hl.agg.group_by(
                    mt.annotation,
                    hl.agg.filter(
                        mt.SE != 0,
                        hl.agg.count_where(
                            mt.Pvalue < EMPIRICAL_P_THRESHOLDS["variant"]
                        ),
                    ),
                ),
            )
            mt = mt.annotate_cols(
                **{
                    f"{annotation}_sig_var_cnt": mt.sig_var_cnt.get(annotation)
                    for annotation in ANNOTATIONS
                },
            )
        else:
            mt = mt.annotate_rows(
                all_sig_pheno_cnt=hl.agg.sum(
                    mt.Pvalue < EMPIRICAL_P_THRESHOLDS["variant"]
                ),
                sig_pheno_cnt=hl.agg.group_by(
                    mt.trait_type2,
                    hl.agg.count_where(mt.Pvalue < EMPIRICAL_P_THRESHOLDS["variant"]),
                ),
            )
            mt = mt.annotate_rows(
                **{
                    f"{trait_type}_sig_pheno_cnt": mt.sig_pheno_cnt.get(trait_type)
                    for trait_type in TRAIT_TYPES
                    for test in TESTS
                },
            )
            mt = mt.annotate_cols(
                all_sig_var_cnt=hl.agg.sum(
                    mt.Pvalue < EMPIRICAL_P_THRESHOLDS["variant"]
                ),
                sig_var_cnt=hl.agg.group_by(
                    mt.annotation,
                    hl.agg.count_where(mt.Pvalue < EMPIRICAL_P_THRESHOLDS["variant"]),
                ),
            )
            mt = mt.annotate_cols(
                **{
                    f"{annotation}_sig_var_cnt": mt.sig_var_cnt.get(annotation)
                    for annotation in ANNOTATIONS
                },
            )
    return mt


def compare_gene_var_sig_cnt_mt(
    test_type: str = "skato", filters: bool = True, tranche: str = CURRENT_TRANCHE
):
    """
    Compare the number of significant association for genes with variants within corresponding genes
    :param str test_type: Type of gene-level association test result, `skato` or `burden`
    :param bool filters: Whether to use QCed sets
    :param str tranche: Tranche of data to use
    :return: Hail MatrixTable with variant-level association count grouped by genes at phenotype-level and gene-level
    :rtype: MatrixTable
    """
    if filters:
        mt = get_qc_result_mt(result_type="gene", test_type=test_type, tranche=tranche)
        var = get_qc_result_mt(
            result_type="variant", test_type=test_type, tranche=tranche
        )
    else:
        mt = load_final_sumstats_table(
            result_type="gene", extension="mt", tranche=tranche
        )
        var = load_final_sumstats_table(
            result_type="variant", extension="mt", tranche=tranche
        )
    if tranche == "500k":
        mt = mt.filter_rows(mt.annotation != "pLoF|missense|LC")
    vep = hl.read_table(var_annotations_ht_path("vep", *TRANCHE_DATA[tranche]))
    vep = process_consequences(vep)
    var = var.annotate_rows(
        gene_id=vep[var.row_key].vep.worst_csq_by_gene_canonical.gene_id
    )
    var = var.explode_rows(var.gene_id)
    var = var.annotate_rows(
        annotation=hl.if_else(
            hl.literal({"missense", "LC"}).contains(var.annotation),
            "missense|LC",
            var.annotation,
        ),
    )
    if filters:
        var = var.group_rows_by("gene_id", "gene", "annotation").aggregate(
            var_cnt=hl.agg.filter(
                var.SE != 0,
                hl.agg.count_where(var.Pvalue < EMPIRICAL_P_THRESHOLDS["variant"]),
            )
        )
    else:
        var = var.group_rows_by("gene_id", "gene", "annotation").aggregate(
            var_cnt=hl.agg.count_where(var.Pvalue < EMPIRICAL_P_THRESHOLDS["variant"])
        )
    pvalue = P_VALUE_FIELDS[test_type.lower()]
    mt = mt.annotate_entries(var_cnt=var[mt.row_key, mt.col_key]["var_cnt"])
    mt = mt.annotate_cols(
        pheno_gene_var_sig_cnt=hl.agg.count_where(
            hl.is_defined(mt.var_cnt)
            & (mt.var_cnt > 0)
            & hl.is_defined(mt[pvalue])
            & (mt[pvalue] < EMPIRICAL_P_THRESHOLDS[test_type.lower()])
        ),
        pheno_var_sig_cnt=hl.agg.count_where(
            (hl.is_defined(mt.var_cnt) & (mt.var_cnt > 0))
            & (
                ~(
                    hl.is_defined(mt[pvalue])
                    & (mt[pvalue] < EMPIRICAL_P_THRESHOLDS[test_type.lower()])
                )
            )
        ),
        pheno_gene_sig_cnt=hl.agg.count_where(
            (
                hl.is_defined(mt[pvalue])
                & (mt[pvalue] < EMPIRICAL_P_THRESHOLDS[test_type.lower()])
            )
            & (~(hl.is_defined(mt.var_cnt) & (mt.var_cnt > 0)))
        ),
        pheno_none_sig_cnt=hl.agg.count_where(
            ~(
                (
                    hl.is_defined(mt[pvalue])
                    & (mt[pvalue] < EMPIRICAL_P_THRESHOLDS[test_type.lower()])
                )
                | (hl.is_defined(mt.var_cnt) & (mt.var_cnt > 0))
            )
        ),
        trait_type2=hl.if_else(
            mt.trait_type == "icd_first_occurrence", "icd10", mt.trait_type
        ),
    )

    mt = mt.annotate_rows(
        gene_var_sig_cnt=hl.agg.count_where(
            hl.is_defined(mt.var_cnt)
            & (mt.var_cnt > 0)
            & hl.is_defined(mt[pvalue])
            & (mt[pvalue] < EMPIRICAL_P_THRESHOLDS[test_type.lower()])
        ),
        var_sig_cnt=hl.agg.count_where(
            (hl.is_defined(mt.var_cnt) & (mt.var_cnt > 0))
            & (
                ~(
                    hl.is_defined(mt[pvalue])
                    & (mt[pvalue] < EMPIRICAL_P_THRESHOLDS[test_type.lower()])
                )
            )
        ),
        gene_sig_cnt=hl.agg.count_where(
            (
                hl.is_defined(mt[pvalue])
                & (mt[pvalue] < EMPIRICAL_P_THRESHOLDS[test_type.lower()])
            )
            & (~(hl.is_defined(mt.var_cnt) & (mt.var_cnt > 0)))
        ),
        none_sig_cnt=hl.agg.count_where(
            (
                (
                    hl.is_missing(mt[pvalue])
                    | (mt[pvalue] >= EMPIRICAL_P_THRESHOLDS[test_type.lower()])
                )
                & (hl.is_missing(mt.var_cnt) | (mt.var_cnt == 0))
            )
        ),
    )
    return mt


def compute_mean_coverage_ht(tranche: str = CURRENT_TRANCHE):
    """
    Generate mean coverage Table
    :param str tranche: Tranche of data to use
    :return: Mean coverage Hail Table
    :rtype: Table
    """
    int_auto = hl.read_table(
        interval_qc_path(data_source="broad", freeze=7, chrom="autosomes")
    )
    int_sex = hl.read_table(
        interval_qc_path(data_source="broad", freeze=7, chrom="sex_chr")
    )
    int_sex = int_sex.select(
        target_mean_dp=int_sex.target_mean_dp["XX"],
        target_pct_gt_10x=int_sex.target_pct_gt_10x["XX"],
        target_pct_gt_20x=int_sex.target_pct_gt_20x["XX"],
        pct_samples_10x=int_sex.pct_samples_10x["XX"],
        pct_samples_20x=int_sex.pct_samples_20x["XX"],
    )
    int_full = int_auto.union(int_sex)

    var = hl.read_matrix_table(get_results_mt_path("variant", tranche=tranche)).rows()
    vep = hl.read_table(var_annotations_ht_path("vep", *TRANCHE_DATA[tranche]))
    vep = process_consequences(vep)
    var = var.annotate(csq=vep[var.key].vep.worst_csq_by_gene_canonical)
    var = var.explode(var.csq)
    var = var.annotate(
        gene_id=var.csq.gene_id,
        gene_symbol=var.csq.gene_symbol,
        annotation=annotation_case_builder(var.csq),
        coverage=int_full[var.locus].target_mean_dp,
    )
    var = var.annotate(
        annotation=hl.if_else(
            hl.literal({"missense", "LC"}).contains(var.annotation),
            "missense|LC",
            var.annotation,
        ),
    )
    if tranche == "500k":
        sub = var.filter(hl.literal({"missense|LC", "pLoF"}).contains(var.annotation))
        sub = sub.annotate(annotation="pLoF|missense|LC")
        var = sub.union(var)
    mean_coverage = var.group_by("gene_id", "gene_symbol", "annotation").aggregate(
        mean_coverage=hl.agg.mean(var.coverage)
    )
    return mean_coverage


def more_cases_tie_breaker(l, r):
    """
    Tie breaker function prefering phenotypes with more cases used in hl.maximal_independent_set()
    :param l: one phenotype
    :param r: the other phenotype
    :return: integer indicating which phenotype is preferred
    :rtype: int
    """
    return (
        hl.case()
        .when(l.n_cases_defined > r.n_cases_defined, -1)
        .when(l.n_cases_defined == r.n_cases_defined, 0)
        .when(l.n_cases_defined < r.n_cases_defined, 1)
        .or_missing()
    )


def get_corr_phenos_ht(
    r_2: float = None, tie_breaker=None, tranche: str = CURRENT_TRANCHE
):
    """
    Generate a subset of related phenotypes to remove
    :param float r_2: Correlation-square cutoff point, phenotype with correlation above which are included to be removed
    :param tie_breaker: Tie-breaker function to select phenotype based on certain preference
    :param str tranche: Tranche of data to use
    :return: Hail Table of correlated phenotypes to remove
    :rtype: Table
    """
    pheno_mt = get_ukb_pheno_mt()
    pheno_mt = drop_pheno_fields_mt(pheno_mt)
    ht = hl.read_table(get_results_mt_path("pheno", tranche=tranche, extension="ht"))
    pheno_mt = pheno_mt.filter_cols(hl.is_defined(ht[pheno_mt.col_key]))
    corr = make_pairwise_ht(pheno_mt, pheno_field=pheno_mt.both_sexes, correlation=True)
    related = corr.filter((corr.entry ** 2 >= r_2) & (corr.i != corr.j))
    pheno_to_remove = hl.maximal_independent_set(
        related.i_data, related.j_data, keep=False, tie_breaker=tie_breaker
    )
    return pheno_to_remove


def export_ht_to_txt_bgz(
    path_to_ht: str, sub_dir: str, output_filename: str, tranche: str = CURRENT_TRANCHE
):
    """
    Convert Hail Table to txt.bgz
    :param str path_to_ht: Path to input hail table
    :param str sub_dir: sub directory where the Hail Tables are stored
    :param str output_filename: Output file name
    :param str tranche: Tranche of data to use
    """
    ht = hl.read_table(path_to_ht)
    ht.export(f"{public_bucket}/{tranche}/{sub_dir}/{output_filename}.txt.bgz")


def export_all_ht_to_txt_bgz(sub_dir: str, tranche: str = CURRENT_TRANCHE):
    """
    Convert all Hail Tables to txt.bgz
    :param str sub_dir: sub directory where the Hail Tables are stored
    :param str tranche: Tranche of data to use
    """
    files = subprocess.check_output(
        f"gsutil ls {public_bucket}/{tranche}/{sub_dir}/", shell=True
    )
    files = files.decode("utf-8").split("\n")
    import re

    ht_list = [file for file in files if file.endswith(".ht/")]
    for ht in ht_list:
        name = re.split("/|\\.", ht)[-3]
        export_ht_to_txt_bgz(ht, sub_dir, name)


def get_related_pheno_cnt_list(pheno_ht: hl.Table):
    """
    Generate the number of related phenotypes to remove at different correlation cutoff points
    :param Table pheno_ht: Original phenotypes from association test result MatrixTable (hl.read_matrix_table(get_results_mt_path()).cols())
    :return: list of number of related phenotypes to remove
    :rtype: list
    """
    l = []
    pheno_mt = get_ukb_pheno_mt()
    pheno_mt = pheno_mt.filter_cols(hl.is_defined(pheno_ht[pheno_mt.col_key]))
    corr = make_pairwise_ht(pheno_mt, pheno_field=pheno_mt.both_sexes, correlation=True)
    import numpy as np

    for k in np.arange(0.1, 1.0, 0.1):
        print(k)
        related = corr.filter((corr.entry ** 2 >= k) & (corr.i != corr.j))
        pheno_to_remove = hl.maximal_independent_set(
            related.i_data,
            related.j_data,
            keep=False,
            tie_breaker=more_cases_tie_breaker,
        )
        l.append(pheno_to_remove.count())
    return l


def get_icd_min_p_ht(
    result_type: str = "gene",
    test_type: str = "skato",
    filters: bool = True,
    tranche: str = CURRENT_TRANCHE,
):
    """
    Generate min pvalue for each icd10 group from variant-level or gene-level test result
    :param str result_type: Type of association test result, `gene` or `variant`
    :param str test_type: Type of gene-level association test result, `skato` or `burden`
    :param bool filters: Whether to use QCed sets
    :param str tranche: Tranche of data to use
    :return: Hail Table of min p-values for each icd10 group
    :rtype: Table
    """
    pvalue = P_VALUE_FIELDS[test_type.lower()]
    if filters:
        mt = get_qc_result_mt(
            result_type=result_type, test_type=test_type, tranche=tranche
        )
    else:
        mt = load_final_sumstats_table(
            result_type="gene", extension="mt", tranche=tranche
        )

    icd = mt.filter_cols(
        hl.literal({"icd10", "icd_first_occurrence"}).contains(mt.trait_type)
    )
    icd = icd.annotate_cols(
        icd10_group=hl.if_else(
            icd.trait_type == "icd10",
            icd.phenocode[0],
            icd.description.split("first", 2)[0][5],
        ),
        icd10_index=hl.int32(
            hl.if_else(
                icd.trait_type == "icd10",
                icd.phenocode[1:3],
                icd.description.split("first", 2)[0][6:8],
            )
        ),
    )
    icd = icd.annotate_cols(
        icd10_group=hl.case()
        .when((icd.icd10_group == "B"), "A")
        .when((icd.icd10_group == "D") & (icd.icd10_index < 50), "C")
        .when((icd.icd10_group == "H") & (icd.icd10_index < 60), "H1")
        .when((icd.icd10_group == "H") & (icd.icd10_index >= 60), "H2")
        .default(icd.icd10_group)
    )
    ICD10_GROUPS = icd.aggregate_cols(hl.agg.collect_as_set(icd.icd10_group))
    if result_type == "gene":
        icd = icd.annotate_rows(
            min_p=hl.agg.group_by(icd.icd10_group, hl.agg.min(icd[pvalue]))
        )
        icd = icd.select_rows(
            "interval",
            **{
                f"{icd_group}_min_p": icd.min_p.get(icd_group)
                for icd_group in ICD10_GROUPS
            },
        )
    else:
        icd = icd.annotate_rows(
            min_p=hl.agg.group_by(
                icd.icd10_group, hl.agg.filter(icd.SE != 0, hl.agg.min(icd.Pvalue))
            )
        )
        icd = icd.select_rows(
            **{
                f"{icd_group}_min_p": icd.min_p.get(icd_group)
                for icd_group in ICD10_GROUPS
            }
        )
    return icd.rows()


def get_pheno_pvalue_ht(
    result_type: str = "gene",
    test_type: str = "skato",
    coding: list = None,
    phenos_to_keep: hl.Table = None,
    freq_min: float = None,
    freq_max: float = None,
    n_var_min: int = 2,
    coverage_min: int = 20,
    random_phenos: bool = False,
    tranche: str = CURRENT_TRANCHE,
):
    """
    Generate Pvalue table from variant-level or gene-level test result
    :param str result_type: Type of association test result, `gene` or `variant`
    :param str test_type: Type of gene-level association test result, `skato` or `burden`
    :param list coding: List of coding used for phenotype filtering
    :param Table phenos_to_keep: Subset of phenotypes to keep, keyed by ['trait_type', 'phenocode', 'pheno_sex', 'coding', 'modifier']
    :param float freq_min: Lower bound of allele frequency/cumulative allele frequency used for variant/gene filtering
    :param float freq_max: Upper bound of allele frequency/cumulative allele frequency used for variant/gene filtering
    :param int n_var_min: Minimum number of variants included in a gene-based test used for gene filtering
    :param int coverage_min: Minimum coverage used for gene filtering
    :param bool random_pheno: Whether to use randomly-generated phenotypes
    :param str tranche: Tranche of data to use
    :return: Hail Table of p-value
    :rtype: Table
    """
    mt = hl.read_matrix_table(
        get_results_mt_path(result_type, tranche=tranche, random_phenos=random_phenos)
    )
    if coding is not None:
        mt = mt.filter_cols(hl.literal(coding).contains(mt.coding))
    if phenos_to_keep is not None:
        mt = mt.filter_cols(hl.is_defined(phenos_to_keep.index(mt.col_key)))
    if result_type == "gene":  # Gene Info
        mt = annotate_additional_info_mt(
            mt=mt, result_type=result_type, tranche=tranche
        )
        if freq_min is not None:
            mt = mt.filter_rows(mt.CAF > freq_min)
        if freq_max is not None:
            mt = mt.filter_rows(mt.CAF <= freq_max)
        if n_var_min is not None:
            mt = mt.filter_rows(mt.total_variants >= n_var_min)
        if coverage_min is not None:
            mt = mt.filter_rows(mt.mean_coverage >= coverage_min)
        mt = mt.annotate_rows(CAF_range=f"CAF:({freq_min}, {freq_max}]")
        ht = mt.entries()
        pvalue = P_VALUE_FIELDS[test_type.lower()]
        ht = ht.select(pvalue, "CAF_range")
    else:  # Variant Info
        if random_phenos:
            mt = mt.annotate_rows(
                af=hl.agg.collect_as_set(mt.AF),
                annt=hl.agg.collect_as_set(mt.annotation),
            )
            mt = mt.explode_rows(mt.af)
            mt = mt.explode_rows(mt.annt)
            mt = mt.drop("annotation", "AF")
            mt = mt.annotate_rows(annotation=mt.annt, AF=mt.af)
        if freq_min is not None:
            mt = mt.filter_rows(mt.call_stats.AF > freq_min)
        if freq_max is not None:
            mt = mt.filter_rows(mt.call_stats.AF <= freq_max)
        mt = mt.annotate_rows(
            annotation=hl.if_else(
                hl.literal({"missense", "LC"}).contains(mt.annotation),
                "missense|LC",
                mt.annotation,
            ),
            AF_range=f"AF:({freq_min}, {freq_max}]",
        )
        ht = mt.entries()
        ht = ht.select("Pvalue", "AF_range")
    return ht


def get_pvalue_by_freq_interval_ht(
    result_type: str = "gene",
    test_type: str = "skato",
    freq_breaks: list = [0.0001, 0.001, 0.01, 0.1],
    phenos_to_keep: hl.Table = None,
    n_var_min: int = None,
    coverage_min: int = None,
    random_phenos: bool = False,
    tranche: str = CURRENT_TRANCHE,
):
    """
    Generate Pvalue Table from variant-level or gene-level test result by frequency interval
    :param str result_type: Type of association test result, `gene` or `variant`
    :param str test_type: Type of gene-level association test result, `skato` or `burden`
    :param list freq_breaks: Breaks to divide frequency into bins
    :param Table phenos_to_keep: Subset of phenotypes to keep, keyed by ['trait_type', 'phenocode', 'pheno_sex', 'coding', 'modifier']
    :param int n_var_min: Minimum number of variants included in a gene-based test used for gene filtering
    :param int coverage_min: Minimum coverage used for gene filtering
    :param bool random_pheno: Whether to use randomly-generated phenotypes
    :param str tranche: Tranche of data to use
    :return: Hail Table of p-value by frequency interval
    :rtype: Table
    """
    ht = get_pheno_pvalue_ht(
        result_type=result_type,
        test_type=test_type,
        phenos_to_keep=phenos_to_keep,
        freq_min=freq_breaks[-1],
        n_var_min=n_var_min,
        coverage_min=coverage_min,
        random_phenos=random_phenos,
        tranche=tranche,
    )
    for i in list(range(0, len(freq_breaks) - 1)):
        sub_ht = get_pheno_pvalue_ht(
            result_type=result_type,
            test_type=test_type,
            phenos_to_keep=phenos_to_keep,
            freq_min=freq_breaks[i],
            freq_max=freq_breaks[i + 1],
            n_var_min=n_var_min,
            coverage_min=coverage_min,
            random_phenos=random_phenos,
            tranche=tranche,
        )
        ht = ht.union(sub_ht)
    ht = ht.union(
        get_pheno_pvalue_ht(
            result_type=result_type,
            test_type=test_type,
            phenos_to_keep=phenos_to_keep,
            freq_max=freq_breaks[0],
            n_var_min=n_var_min,
            coverage_min=coverage_min,
            random_phenos=random_phenos,
            tranche=tranche,
        )
    )
    return ht


def annotate_clinvar_pathogenicity_ht():
    """
    Re-define pathogenicity group for clinvar variants
    :return: Hail Table with new definitions
    :rtype: Table
    """
    clinvar_ht = clinvar.ht()
    clinvar_ht = clinvar_ht.filter(
        ~clinvar_ht.info.CLNREVSTAT[0].contains("no_assertion")
    )  # remove clinvar zero star variants
    clinvar_ht = clinvar_ht.explode(clinvar_ht.info.CLNSIG)
    clinvar_ht = clinvar_ht.annotate(
        pathogenicity=hl.case()
        .when(
            hl.literal(
                {"Pathogenic", "Likely_pathogenic", "Pathogenic/Likely_pathogenic"}
            ).contains(clinvar_ht.info.CLNSIG),
            "P/LP",
        )
        .when(
            hl.literal({"Uncertain_significance"}).contains(clinvar_ht.info.CLNSIG),
            "VUS",
        )
        .when(
            hl.literal({"Benign", "Likely_benign", "Benign/Likely_benign"}).contains(
                clinvar_ht.info.CLNSIG
            ),
            "B/LB",
        )
        .or_missing()
    )
    return clinvar_ht


def add_liftover_rg37_to_rg38_ht(ht: hl.Table):
    """
    Add liftover to Hail Table from rg37 to rg38
    :param Table ht: Hail Table to add liftover on
    :return: Hail Table
    :rtype: hl.Table
    """
    rg37 = hl.get_reference("GRCh37")
    rg38 = hl.get_reference("GRCh38")
    if not rg37.has_liftover("GRCh38"):
        rg37.add_liftover(
            "gs://hail-common/references/grch37_to_grch38.over.chain.gz", rg38
        )
    ht = ht.annotate(new_locus=hl.liftover(ht.locus, "GRCh38"))
    ht = ht.filter(hl.is_defined(ht.new_locus))
    ht = ht.key_by(locus=ht.new_locus, alleles=ht.alleles)
    return ht


def annotate_additional_info_mt(
    mt: hl.MatrixTable,
    result_type: str = "gene",
    expected_AC_only: bool = False,
    expected_AC_min: int = 50,
    tranche: str = CURRENT_TRANCHE,
):
    """
    Attach additional information to the original gene/variant-level sumstats MatrixTable
    :param MatrixTable mt: sumstats MatrixTable to annotate informtaion to
    :param str result_type: Type of association test result, `gene` or `variant`
    :param bool expected_AC_only: only annotate expected AC for variant results.
    :param int expected_AC_min: the lower bound used to remove genes/variants with 0 phenotypes having expected AC greater than this.
    :param str tranche: Tranche of data to use
    :return: Hail MatrixTable with updated information
    :rtype: hl.MatrixTable
    """
    case_field = "n_cases_defined" if tranche == "500k" else "n_cases"
    if result_type == "gene":
        caf_ht = hl.read_table(get_util_info_path("caf", tranche=tranche))
        coverage_ht = hl.read_table(get_util_info_path("coverage", tranche=tranche))
        mt = mt.annotate_rows(
            CAF=caf_ht[mt.annotation, mt.gene_id]["CAF"],
            mean_coverage=coverage_ht[mt.row_key].mean_coverage,
        )
        mt = mt.annotate_entries(expected_AC=mt.CAF * mt[case_field])
    else:
        if tranche == "500k":
            mt = mt.annotate_entries(expected_AC=mt.call_stats.AF * mt[case_field])
        else:
            mt = mt.annotate_entries(expected_AC=mt.AF * mt[case_field])
        if not expected_AC_only:
            clinvar_ht = annotate_clinvar_pathogenicity_ht()

            ht = hl.read_table(
                "gs://gcp-public-data--gnomad/papers/2019-tx-annotation/pre_computed/all.possible.snvs.tx_annotated.021520.ht"
            )
            ht = ht.explode(ht.tx_annotation)
            ht = add_liftover_rg37_to_rg38_ht(ht)

            vep = hl.read_table(var_annotations_ht_path("vep", *TRANCHE_DATA[tranche]))
            vep = process_consequences(vep)
            vep = vep.explode(vep.vep.worst_csq_by_gene_canonical)
            vep = vep.annotate(
                tx_annotation_csq=ht[vep.key].tx_annotation.csq,
                mean_proportion=ht[vep.key].tx_annotation.mean_proportion,
            )
            mt = mt.annotate_rows(
                pathogenicity=clinvar_ht[mt.locus, mt.alleles]["pathogenicity"],
                annotation=hl.if_else(
                    hl.literal({"missense", "LC"}).contains(mt.annotation),
                    "missense|LC",
                    mt.annotation,
                ),
                polyphen2=vep[
                    mt.row_key
                ].vep.worst_csq_by_gene_canonical.polyphen_prediction,
                amino_acids=vep[mt.row_key].vep.worst_csq_by_gene_canonical.amino_acids,
                lof=vep[mt.row_key].vep.worst_csq_by_gene_canonical.lof,
                most_severe_consequence=vep[
                    mt.row_key
                ].vep.worst_csq_by_gene_canonical.most_severe_consequence,
                tx_annotation_csq=vep[mt.row_key].tx_annotation_csq,
                mean_proportion_expressed=vep[mt.row_key].mean_proportion,
            )
    mt = mt.annotate_rows(
        expected_ac_row_filter=hl.agg.count_where(mt.expected_AC >= expected_AC_min)
    )
    mt = mt.annotate_cols(
        expected_ac_col_filter=hl.agg.count_where(mt.expected_AC >= expected_AC_min)
    )
    return mt


def annotate_synonymous_lambda_ht(lambda_ht, test_type):
    """
    Annotate lambda gc of synonymous gene to the corresponding pLoF and missense gene
    :param Table lambda_ht: original lambda gc table
    :param str test_type: Type of gene-level association test result, `skato`, `skat` or `burden`
    :return: Hail Table with synonymous lambda gc
    :rtype: hl.Table
    """
    assert test_type.lower() in TESTS, "Invalid test type"
    lambda_ht_syn = lambda_ht.filter(lambda_ht.annotation == "synonymous")
    lambda_ht_syn = lambda_ht_syn.key_by("gene_id", "gene_symbol")
    lambda_ht = lambda_ht.annotate(
        **{
            f"synonymous_lambda_gc_{test_type.lower()}": lambda_ht_syn.index(
                lambda_ht.gene_id, lambda_ht.gene_symbol
            )[f"annotation_lambda_gc_{test_type.lower()}"]
        }
    )
    return lambda_ht


def load_final_sumstats_table(
    result_type: str, extension: str, tranche: str = CURRENT_TRANCHE
):
    """
    Load summary statistics table with all qc metrics
    :param str result_type: Type of association test result, `gene` or `variant`
    :param str extension: Extension of output file
    :param str tranche: Tranche of data to use
    :return: hl.Table or hl.MatrixTable with all qc metrics
    :rtype: hl.Table or hl.MatrixTable
    """
    path = f"{public_bucket}/{tranche}/qc/{result_type}_qc_metrics_ukb_exomes_{tranche}.{extension}"
    if extension == "mt":
        table = hl.read_matrix_table(path)
    else:
        table = hl.read_table(path)
    return table


def get_qc_result_mt(
    result_type: str = "gene", test_type: str = "skato", tranche: str = CURRENT_TRANCHE
):
    """
    Apply all filters and generate Completely QCed data
    :param str result_type: Type of association test result, `gene` or `variant`
    :param str test_type: Type of gene-level association test result, `skato`, `skat` or `burden`
    :param str tranche: Tranche of data to use
    :return: Completely QCed Hail MatrixTable
    :rtype: hl.MatrixTable
    """
    mt = load_final_sumstats_table(
        result_type=result_type, extension="mt", tranche=tranche
    )
    mt = mt.filter_cols(mt[f"keep_pheno_{test_type.lower()}"] & mt.keep_pheno_unrelated)
    if result_type == "gene":
        mt = mt.filter_rows(
            mt[f"keep_gene_{test_type.lower()}"]
            & mt.keep_gene_coverage
            & mt.keep_gene_expected_ac
            & mt.keep_gene_n_var
        )
    else:
        mt = mt.filter_rows(mt.keep_var_expected_ac & mt.keep_var_annt)
    mt = mt.filter_entries(mt.keep_entry_expected_ac)
    return mt


def modify_phenos_mt(result_type):
    """
    Modify and remove phenotypes from the original result Matrixtable
    :param str result_type: Type of association test result, `gene` or `variant`
    :return: Hail MatrixTable with phenotypes removed
    :rtype: hl.MatrixTable
    """
    mt = hl.read_matrix_table(
        f"{bucket}/{CURRENT_TRANCHE}/results/{'' if result_type=='gene' else 'variant_'}results.mt"
    )
    mt = drop_pheno_fields_mt(mt)
    mt = mt.filter_cols(
        (mt.n_cases_defined >= 100)
        & ~((mt.phenocode == "20004") & ((mt.coding == "1490") | (mt.coding == "1540")))
        & ~((mt.phenocode == "Allergy_pfe") | (mt.phenocode == "AnyAutoimmune_pfe"))
    )
    mt = mt.annotate_cols(
        description=hl.case()
        .when(mt.description.matches("AbbVie"), mt.description.replace("AbbVie ", ""))
        .when(mt.description.matches("pfe"), mt.description.replace(" \(pfe\)", ""))
        .default(mt.description)
    )
    mt = mt.annotate_cols(
        description=hl.if_else(
            mt.phenocode.startswith("WBFMadjBMI_")
            | mt.phenocode.startswith("WBfatmass_"),
            mt.description.replace("fat", "fat free"),
            mt.description,
        ),
        description_more=hl.if_else(
            mt.phenocode.startswith("WBFMadjBMI_")
            | mt.phenocode.startswith("WBfatmass_"),
            mt.description_more.replace("mass", "free mass"),
            mt.description_more,
        ),
    )
    mt = mt.annotate_cols(
        description=hl.if_else(
            mt.description.matches("WBFMadjBMI"),
            mt.description.replace("WBFMadjBMI", "WBFFMadjBMI"),
            mt.description,
        )
    )

    import random

    random.seed(2022)
    alz_order = [1, 2]
    random.shuffle(alz_order)
    ibd_order = [1, 2]
    random.shuffle(ibd_order)

    mt = mt.key_cols_by(
        phenocode=hl.case()
        .when(mt.phenocode == "AbbVie_Alzheimers", f"Alzheimers_custom{alz_order[0]}")
        .when(mt.phenocode == "Alzheimers_BI", f"Alzheimers_custom{alz_order[1]}")
        .when(mt.phenocode == "AbbVie_IBD", f"IBD_custom{ibd_order[0]}")
        .when(mt.phenocode == "IBD_pfe", f"IBD_custom{ibd_order[1]}")
        .when(
            mt.phenocode.startswith("AbbVie_"),
            mt.phenocode.replace("AbbVie_", "") + "_custom",
        )
        .when(mt.phenocode.endswith("_pfe"), mt.phenocode.replace("_pfe", "_custom"))
        .when(mt.phenocode.endswith("_BI"), mt.phenocode.replace("_BI", "_custom"))
        .default(mt.phenocode),
    )
    mt = mt.key_cols_by(
        trait_type=mt.trait_type,
        phenocode=hl.case()
        .when(
            mt.phenocode.startswith("WBFMadjBMI_"),
            mt.phenocode.replace("WBFMadjBMI_", "WBFFMadjBMI_"),
        )
        .when(
            mt.phenocode.startswith("WBfatmass_"),
            mt.phenocode.replace("WBfatmass_", "WBfatfreemass_"),
        )
        .default(mt.phenocode),
        pheno_sex=mt.pheno_sex,
        coding=mt.coding,
        modifier=hl.if_else(
            hl.set({"biogen", "abbvie", "pfizer"}).contains(mt.modifier),
            "custom",
            mt.modifier,
        ),
    )
    return mt


def drop_pheno_fields_mt(mt: hl.MatrixTable):
    """
    Drop a set of fieldss from the original result Matrixtable
    :param MatrixTable mt: Input original result Matrixtable
    :return: Hail MatrixTable with fields removed
    :rtype: hl.MatrixTable
    """
    return mt.drop("Abbvie_Priority", "Biogen_Priority", "Pfizer_Priority", "score")

def modify_n_case_control(mt):
    """
    Reannotate n_cases and n_controls
    :param MatrixTable mt: Input original result Matrixtable
    :return: hl.MatrixTable with fields defined
    :rtype: hl.MatrixTable
    """
    mt = mt.annotate_cols(n_cases = hl.if_else(hl.is_defined(mt.n_cases), mt.n_cases, mt['N.Cases']),
                          n_controls = hl.if_else(hl.is_defined(mt.n_controls), mt.n_controls, mt['N.Controls']))
    mt = mt.drop('N.Cases', 'N.Controls')
    return mt


def annotate_pheno_qc_metric_mt(
    mt: hl.MatrixTable, lambda_min: float = 0.75, tranche: str = CURRENT_TRANCHE
):
    """
    Annotate QC metrics to phenotype table
    :param MatrixTable mt: Input original result Matrixtable
    :param float lambda_min: Lower bound of lambda GC
    :param str tranche: Tranche of data to use
    :return: Hail MatrixTable with phenotype QC metrics
    :rtype: MatrixTable
    """
    pheno_corr = hl.read_table(
        get_ukb_exomes_sumstat_path(
            subdir="qc", dataset="correlated", result_type="phenos", tranche=tranche
        )
    )
    mt = mt.annotate_cols(
        **{
            f"keep_pheno_{test}": (mt[f"lambda_gc_{test}"] > lambda_min)
            for test in TESTS
        },
        keep_pheno_unrelated=hl.is_missing(
            pheno_corr.key_by(
                **{field: pheno_corr.node[field] for field in mt.col_key}
            )[mt.col_key]
        ),
    )
    return mt


def annotate_qc_metric_mt(
    result_type: str,
    lambda_min: float = 0.75,
    coverage_min: int = 20,
    expected_AC_min: int = 50,
    n_var_min: int = 2,
    tranche: str = CURRENT_TRANCHE,
):
    """
    Annotate QC metrics to gene table
    :param str result_type: Type of association test result, `gene` or `variant`
    :param float lambda_min: Lower bound of synonymous lambda GC
    :param int coverage_min: Minimum coverage
    :param int expected_AC_min: the lower bound used to remove genes/variants with 0 phenotypes having expected AC greater than this.
    :param int n_var_min: Minimum number of variants
    :param str tranche: Tranche of data to use
    :return: Hail MatrixTable with gene QC metrics
    :rtype: MatrixTable
    """
    mt = hl.read_matrix_table(
        get_results_mt_path(result_type=result_type, tranche=tranche)
    )
    mt = annotate_additional_info_mt(
        mt,
        result_type=result_type,
        expected_AC_only=True,
        expected_AC_min=expected_AC_min,
        tranche=tranche,
    )
    pheno_lambda = hl.read_table(
        get_ukb_exomes_sumstat_path(
            subdir="qc/lambda_gc",
            dataset="lambda_by_pheno_full_filtered",
            tranche=tranche,
        )
    )
    pheno_lambda = pheno_lambda.select(*{f"lambda_gc_{test}" for test in TESTS})
    mt = mt.annotate_cols(**pheno_lambda[mt.col_key])
    if result_type == "gene":
        gene_lambda = hl.read_table(
            get_ukb_exomes_sumstat_path(
                subdir="qc/lambda_gc",
                dataset="lambda_by_gene",
                result_type="",
                tranche=tranche,
            )
        )
        gene_lambda = gene_lambda.rename(
            {f"all_lambda_gc_{test}": f"annotation_lambda_gc_{test}" for test in TESTS}
        )
        gene_lambda = gene_lambda.drop(*{f"lambda_gc_{test}" for test in TESTS})
        for test in TESTS:
            gene_lambda = annotate_synonymous_lambda_ht(gene_lambda, test)

        mt = mt.annotate_rows(**gene_lambda[mt.row_key])
        mt = mt.annotate_rows(
            **{
                f"keep_gene_{test}": (mt[f"synonymous_lambda_gc_{test}"] > lambda_min)
                for test in TESTS
            },
            keep_gene_coverage=mt.mean_coverage >= coverage_min,
            keep_gene_expected_ac=mt.expected_ac_row_filter > 0,
            keep_gene_n_var=mt.total_variants >= n_var_min,
        )
        mt = mt.annotate_globals(
            coverage_min=coverage_min,
            expected_AC_min=expected_AC_min,
            n_var_min=n_var_min,
            gene_syn_lambda_min=lambda_min,
        )
    else:
        mt = mt.annotate_rows(
            annotation=hl.if_else(
                hl.literal({"missense", "LC"}).contains(mt.annotation),
                "missense|LC",
                mt.annotation,
            ),
            keep_var_expected_ac=mt.expected_ac_row_filter > 0,
            keep_var_annt=hl.is_defined(mt.annotation),
        )
        mt = mt.annotate_globals(
            expected_AC_min=expected_AC_min,
        )
    mt = annotate_pheno_qc_metric_mt(mt, lambda_min=lambda_min)
    mt = mt.annotate_globals(pheno_lambda_min=lambda_min)
    mt = mt.annotate_entries(keep_entry_expected_ac = mt.expected_AC >= expected_AC_min)

    return mt

def get_constrained_con_pheno_group_ht(test_type:str='skato', tranche:str = CURRENT_TRANCHE):
    ht = hl.import_table(
        'gs://gcp-public-data--gnomad/release/2.1.1/constraint/gnomad.v2.1.1.lof_metrics.by_gene.txt.bgz',
        key=['gene', 'gene_id'], types={'oe_lof_upper_bin': hl.tint64})
    mt = get_qc_result_mt('gene', test_type, tranche)
    mt = mt.filter_rows(mt.annotation != 'pLoF|missense|LC')
    mt = mt.filter_cols(mt.trait_type == 'continuous')
    mt = mt.annotate_cols(pheno_group=hl.case()
                          .when(mt.category.matches('Touchscreen'), 'Touchscreen')
                          .when(mt.category.matches('Physical measures'), 'Physical measures')
                          .when(mt.category.matches('Imaging'), 'Imaging')
                          .when((mt.phenocode.startswith('30')) & (hl.len(mt.phenocode) == 5), 'Biomarkers')
                          .or_missing())
    mt = mt.filter_cols(hl.is_defined(mt.pheno_group))
    mt = mt.annotate_entries(pheno_sig=hl.if_else(mt[P_VALUE_FIELDS[test_type]] <
                                                  EMPIRICAL_P_THRESHOLDS[test_type], 1, 0))
    mt = mt.annotate_entries(pheno_sig=hl.if_else(hl.is_missing(mt.pheno_sig), 0, mt.pheno_sig))
    sub = mt.group_cols_by(mt.pheno_group).aggregate(pheno_group_sig_cnt=hl.agg.sum(mt.pheno_sig))
    sub = sub.select_rows('CAF')
    return sub

