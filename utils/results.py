from ukbb_qc.resources.sample_qc import interval_qc_path
from ukbb_qc.resources.variant_qc import var_annotations_ht_path
from ukbb_qc.resources.basics import release_ht_path
from ukb_common.utils.annotations import annotation_case_builder
from gnomad.resources.grch38.reference_data import clinvar
from gnomad.utils.vep import process_consequences
from ..resources import *

ANNOTATIONS = ('pLoF', 'missense|LC', 'synonymous')
TESTS = ('skato', 'skat', 'burden')
TRAIT_TYPES = ('continuous', 'categorical', 'icd10')
P_VALUE_FIELDS = {'skato': 'Pvalue', 'skat': 'Pvalue_SKAT', 'burden': 'Pvalue_Burden'}
EMPIRICAL_P_THRESHOLDS = {'skato': 2.5e-8, 'burden': 6.7e-7, 'variant': 8e-9}
sumstats_folder = f'{public_bucket}summary_statistics_analysis'
lambda_folder = f'{sumstats_folder}/lambda_gc'
pvalue_folder = f'{sumstats_folder}/pvalue'
sig_folder = f'{sumstats_folder}/sig_cnt'
beta_folder = f'{sumstats_folder}/beta'

def get_ukb_exomes_sumstat_path(subdir:str, dataset:str, result_type: str = 'gene', extension: str = 'ht', tranche: str = CURRENT_TRANCHE):
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
    return f'{subdir}/{dataset}_{result_type}_{tranche}.{extension}'

def get_util_info_path(type:str, tranche: str = CURRENT_TRANCHE):
    """
    Get cumulative allele frequency or mean coverage Table path
    :param str type: Type of table to use, `coverage` or `other`
    :param str tranche: Tranche of data to use
    :return: Path to information Table
    :rtype: str
    """
    type = '_mean_coverage' if type=='coverage' else '_caf'
    return f'{sumstats_folder}/gene{type}_{tranche}.ht'

def compute_lambda_gc_ht(result_type: str = 'gene', by_annotation: bool = False, by_gene: bool = False, freq_lower: float = None, freq_upper: float = None, 
                         n_var_min: int = None, coverage_min: int=None, var_filter: bool = False, random_phenos: bool = False, tranche: str = CURRENT_TRANCHE):
    """
    Compute Lambda GC for variant-level test or gene-level tests (SKAT-O, SKAT, and Burden Test) under specified condition
    :param str result_type: Type of association test result, `gene` or `variant`
    :param bool by_annotation: Whether to group by annotation categories
    :param bool by_gene: Whether to compute lambda GC for each gene, applicable only for result_type == `gene`
    :param float freq_lower: Lower bound of allele frequency/cumulative allele frequency used for variant/gene filtering
    :param float freq_upper: Upper bound of allele frequency/cumulative allele frequency used for variant/gene filtering
    :param int n_var_min: Minimum number of variants included in a gene-based test used for gene filtering
    :param int coverage_min: Minimum coverage used for gene filtering
    :param bool var_filter: Whether to filter out test statistics with SE=0
    :param bool random_pheno: Whether to use randomly-generated phenotypes
    :param str tranche: Tranche of data to use
    :return: Hail table with lambda GC
    :rtype: Table
    """
    mt = hl.read_matrix_table(get_results_mt_path(result_type, tranche, random_phenos=random_phenos))
    if result_type == 'gene':  # Gene Info
        mt = annotate_additional_info_mt(result_mt=mt, result_type=result_type, tranche=tranche)
        if freq_lower is not None:
            mt = mt.filter_rows(mt.CAF > freq_lower)
        if freq_upper is not None:
            mt = mt.filter_rows(mt.CAF <= freq_upper)
        if n_var_min is not None:
            mt = mt.filter_rows(mt.total_variants >= n_var_min)
        if coverage_min is not None:
            mt = mt.filter_rows(mt.mean_coverage > coverage_min)
        if by_annotation:
            mt = mt.select_cols('n_cases', 
                                lambda_gc_skato=hl.agg.group_by(mt.annotation, hl.methods.statgen._lambda_gc_agg(mt.Pvalue)), 
                                lambda_gc_skat=hl.agg.group_by(mt.annotation, hl.methods.statgen._lambda_gc_agg(mt.Pvalue_SKAT)), 
                                lambda_gc_burden=hl.agg.group_by(mt.annotation, hl.methods.statgen._lambda_gc_agg(mt.Pvalue_Burden)), 
                                CAF_range=f'CAF:({freq_lower}, {freq_upper}]')
            mt = mt.annotate_cols(**{f'{annotation}_lambda_gc_{test}': mt[f'lambda_gc_{test}'].get(annotation) for annotation in ANNOTATIONS for test in TESTS})
            ht = mt.cols()
            ht = ht.annotate(trait_type2=hl.if_else(ht.trait_type == 'icd_first_occurrence', 'icd10', ht.trait_type), )
        elif by_gene:
            mt = mt.annotate_cols(trait_type2=hl.if_else(mt.trait_type == 'continuous', mt.trait_type, 'categorical'))
            mt = mt.select_rows('CAF', 'mean_coverage',
                                lambda_gc_skato=hl.agg.group_by(mt.trait_type2, hl.methods.statgen._lambda_gc_agg(mt.Pvalue)), 
                                lambda_gc_skat=hl.agg.group_by(mt.trait_type2, hl.methods.statgen._lambda_gc_agg(mt.Pvalue_SKAT)), 
                                lambda_gc_burden=hl.agg.group_by(mt.trait_type2, hl.methods.statgen._lambda_gc_agg(mt.Pvalue_Burden)), )
            mt = mt.annotate_rows(**{f'{trait_type}_lambda_gc_{test}': mt[f'lambda_gc_{test}'].get(trait_type) for trait_type in TRAIT_TYPES for test in TESTS},
                                  all_lambda_gc_skato=hl.agg.filter(hl.is_defined(mt.Pvalue), hl.methods.statgen._lambda_gc_agg(mt.Pvalue)), 
                                  all_lambda_gc_skat=hl.agg.filter(hl.is_defined(mt.Pvalue_SKAT), hl.methods.statgen._lambda_gc_agg(mt.Pvalue_SKAT)), 
                                  all_lambda_gc_burden=hl.agg.filter(hl.is_defined(mt.Pvalue_Burden), hl.methods.statgen._lambda_gc_agg(mt.Pvalue_Burden)), )
            ht = mt.rows()
        else:
            mt = mt.select_cols('n_cases', 
                                lambda_gc_skato=hl.agg.filter(hl.is_defined(mt.Pvalue), hl.methods.statgen._lambda_gc_agg(mt.Pvalue)), 
                                lambda_gc_skat=hl.agg.filter(hl.is_defined(mt.Pvalue_SKAT), hl.methods.statgen._lambda_gc_agg(mt.Pvalue_SKAT)), 
                                lambda_gc_burden=hl.agg.filter(hl.is_defined(mt.Pvalue_Burden), hl.methods.statgen._lambda_gc_agg(mt.Pvalue_Burden)), 
                                CAF_range=f'CAF:({freq_lower}, {freq_upper}]')
            ht = mt.cols()
            ht = ht.annotate(trait_type2=hl.if_else(ht.trait_type == 'icd_first_occurrence', 'icd10', ht.trait_type), )
    else:  # Variant Info
        if freq_lower is not None:
            mt = mt.filter_rows(mt.AF > freq_lower)
        if freq_upper is not None:
            mt = mt.filter_rows(mt.AF <= freq_upper)
        mt = mt.annotate_rows(annotation=hl.if_else(hl.literal({'missense', 'LC'}).contains(mt.annotation), 'missense|LC', mt.annotation))
        if var_filter:
            mt = mt.filter_rows(hl.is_defined(mt.annotation))
            if by_annotation:
                lambda_gc = hl.agg.group_by(mt.annotation, hl.agg.filter(mt.SE!=0, hl.methods.statgen._lambda_gc_agg(mt.Pvalue)))
                mt = mt.select_cols('n_cases', lambda_gc=lambda_gc, AF_range=f'AF:({freq_lower}, {freq_upper}]')
                mt = mt.annotate_cols(**{f'{annotation}_lambda_gc_variant': mt.lambda_gc.get(annotation) for annotation in ANNOTATIONS}, )
            else:
                lambda_gc = hl.agg.filter(mt.SE!=0, hl.methods.statgen._lambda_gc_agg(mt.Pvalue))
                mt = mt.select_cols('n_cases', lambda_gc=lambda_gc, AF_range=f'AF:({freq_lower}, {freq_upper}]')
        else:
            if by_annotation:
                lambda_gc = hl.agg.group_by(mt.annotation, hl.methods.statgen._lambda_gc_agg(mt.Pvalue))
                mt = mt.select_cols('n_cases', lambda_gc=lambda_gc, AF_range=f'AF:({freq_lower}, {freq_upper}]')
                mt = mt.annotate_cols(**{f'{annotation}_lambda_gc_variant': mt.lambda_gc[annotation] for annotation in ANNOTATIONS}, )
            else:
                lambda_gc = hl.agg.filter(hl.is_defined(mt.Pvalue), hl.methods.statgen._lambda_gc_agg(mt.Pvalue))
                mt = mt.select_cols('n_cases', lambda_gc=lambda_gc, AF_range=f'AF:({freq_lower}, {freq_upper}]')
        ht = mt.cols()
        ht = ht.annotate(trait_type2=hl.if_else(ht.trait_type == 'icd_first_occurrence', 'icd10', ht.trait_type), )
    return ht

def compute_lambdas_by_freq_interval_ht(result_type='gene', by_annotation: bool = False, freq_breaks: list=[0.0001, 0.001, 0.01, 0.1], n_var_min: int = None, coverage_min: int=None,
                                        var_filter: bool = False, random_phenos: bool = False, tranche: str = CURRENT_TRANCHE):
    """
    Compute Lambda GC for variant-level test or gene-level tests (SKAT-O, SKAT, and Burden Test) by frequency intervals
    :param str result_type: Type of association test result, `gene` or `variant`
    :param bool by_annotation: Whether to group by annotation categories
    :param list freq_breaks: Breaks to divide frequency into bins
    :param int n_var_min: Minimum number of variants included in a gene-based test used for gene filtering
    :param int coverage_min: Minimum coverage used for gene filtering
    :param bool var_filter: Whether to filter out test statistics with SE=0
    :param bool random_pheno: Whether to use randomly-generated phenotypes
    :param str tranche: Tranche of data to use
    :return: Hail table with lambda GC by frequency intervals
    :rtype: Table
    """
    ht = compute_lambda_gc_ht(result_type=result_type, by_annotation=by_annotation, freq_upper=freq_breaks[0], n_var_min=n_var_min, coverage_min=coverage_min, random_phenos=random_phenos, tranche=tranche)
    for i in list(range(0, len(freq_breaks)-1)):
        sub_ht = compute_lambda_gc_ht(result_type=result_type, by_annotation=by_annotation, freq_lower=freq_breaks[i], freq_upper=freq_breaks[i+1], n_var_min=n_var_min, coverage_min=coverage_min, var_filter=var_filter, random_phenos=random_phenos, tranche=tranche)
        ht = ht.union(sub_ht)
    ht = ht.union(compute_lambda_gc_ht(result_type=result_type, by_annotation=by_annotation, freq_lower=freq_breaks[-1], n_var_min=n_var_min, coverage_min=coverage_min, var_filter=var_filter, random_phenos=random_phenos, tranche=tranche))
    return ht

def compute_lambdas_by_expected_ac_ht(freq_lower: float=None, ac_breaks: list=[1, 10, 100, 1000, 10000, 100000], var_filter: bool = False, random_phenos: bool = False, tranche: str = CURRENT_TRANCHE):
    """
    Compute Lambda GC for variant-level test by expected allele count interval
    :param float freq_lower: Lower bound of allele frequency used for variant filtering
    :param list ac_breaks: Breaks to divide expected allele count into bins
    :param bool var_filter: Whether to filter out test statistics with SE=0
    :param bool random_pheno: Whether to use randomly-generated phenotypes
    :param str tranche: Tranche of data to use
    :return: Hail table with lambda GC by expected allele count intervals
    :rtype: Table
    """
    mt = hl.read_matrix_table(get_results_mt_path('variant', random_phenos=random_phenos))
    if freq_lower is not None:
        mt = mt.filter_rows(mt.AF > freq_lower)
    if random_phenos:
        mt = mt.annotate_rows(AF2=hl.agg.take(mt.AF, 1)[0]).drop('AF').rename({'AF2': 'AF'})
    if var_filter:
        mt = mt.filter_rows(hl.is_defined(mt.annotation))
        mt = mt.annotate_entries(expected_AC=mt.AF * mt.n_cases)
        mt = mt.annotate_cols(lambda_gc0 = hl.agg.filter((mt.expected_AC <= ac_breaks[0]) & (mt.SE != 0), hl.methods.statgen._lambda_gc_agg(mt.Pvalue)))
        for i in list(range(0, len(ac_breaks)-1)):
            mt = mt.annotate_cols(**{f'lambda_gc{ac_breaks[i]}': hl.agg.filter((mt.expected_AC > ac_breaks[i]) & (mt.expected_AC <= ac_breaks[i+1]) & (mt.SE != 0), hl.methods.statgen._lambda_gc_agg(mt.Pvalue))})
        mt = mt.annotate_cols(**{f'lambda_gc{ac_breaks[-1]}' : hl.agg.filter((mt.expected_AC > ac_breaks[-1])& (mt.SE != 0), hl.methods.statgen._lambda_gc_agg(mt.Pvalue))},
                              trait_type2=hl.if_else(mt.trait_type == 'icd_first_occurrence', 'icd10', mt.trait_type))
    else:
        mt = mt.annotate_entries(expected_AC=mt.AF * mt.n_cases)
        mt = mt.annotate_cols(lambda_gc0 = hl.agg.filter(mt.expected_AC <= ac_breaks[0], hl.methods.statgen._lambda_gc_agg(mt.Pvalue)))
        for i in list(range(0, len(ac_breaks)-1)):
            mt = mt.annotate_cols(**{f'lambda_gc{ac_breaks[i]}': hl.agg.filter((mt.expected_AC > ac_breaks[i]) & (mt.expected_AC <= ac_breaks[i+1]), hl.methods.statgen._lambda_gc_agg(mt.Pvalue))})
        mt = mt.annotate_cols(**{f'lambda_gc{ac_breaks[-1]}' : hl.agg.filter(mt.expected_AC > ac_breaks[-1], hl.methods.statgen._lambda_gc_agg(mt.Pvalue))},
                              trait_type2=hl.if_else(mt.trait_type == 'icd_first_occurrence', 'icd10', mt.trait_type))
    return mt.cols()

def write_lambda_hts(result_type='gene', freq_lower: float = None, n_var_min: int = None, coverage_min: int = None, var_filter: bool = False,
                     random_phenos: bool = False, extension: str = 'ht', overwrite: bool = False, tranche: str = CURRENT_TRANCHE):
    """
    Write lambda GC tables for variant-level test or gene-level tests (SKAT-O, SKAT, and Burden Test) under different condition
    :param str result_type: Type of association test result, `gene` or `variant`
    :param float freq_lower: Lower bound of allele frequency/cumulative allele frequency used for variant/gene filtering
    :param int n_var_min: Minimum number of variants included in a gene-based test used for gene filtering
    :param int coverage_min: Minimum coverage used for gene filtering
    :param bool var_filter: Whether to filter out test statistics with SE=0
    :param bool random_pheno: Whether to use randomly-generated phenotypes
    :param str extension: Extension of output file
    :param bool overwrite: Whether to overwrite exsiting file
    :param str tranche: Tranche of data to use
    """
    kwargs = {
        'result_type': result_type,
        'n_var_min': n_var_min,
        'coverage_min': coverage_min,
        'var_filter': var_filter,
        'random_phenos': random_phenos,
    }
    pathargs = {
        'subdir': lambda_folder,
        'result_type': result_type,
        'tranche': tranche,
        'extension': extension,
    }
    rp = 'rp_' if random_phenos else ''
    filter = '_filtered' if any(v is not None for v in [freq_lower, n_var_min, coverage_min]) else ''
    freq_breaks = [freq_lower, 0.001, 0.01, 0.1] if freq_lower is not None else [0.0001, 0.001, 0.01, 0.1]

    compute_lambda_gc_ht(freq_lower=freq_lower, **kwargs).write(get_ukb_exomes_sumstat_path(dataset=f'{rp}lambda_full{filter}', **pathargs), overwrite=overwrite)
    compute_lambda_gc_ht(by_annotation=True, freq_lower=freq_lower, **kwargs).write(get_ukb_exomes_sumstat_path(dataset=f'{rp}lambda_annt{filter}', **pathargs), overwrite=overwrite)
    compute_lambdas_by_freq_interval_ht(freq_breaks=freq_breaks, **kwargs).write(get_ukb_exomes_sumstat_path(dataset=f'{rp}lambda_freq{filter}', **pathargs), overwrite=overwrite)
    compute_lambdas_by_freq_interval_ht(by_annotation=True, freq_breaks=freq_breaks, **kwargs).write(get_ukb_exomes_sumstat_path(dataset=f'{rp}lambda_freq_annt{filter}', **pathargs), overwrite=overwrite)

    if result_type=='gene':
        compute_lambda_gc_ht(by_gene=True, **kwargs).write(get_ukb_exomes_sumstat_path(subdir=lambda_folder, dataset=f'{rp}lambda_by_gene{filter}', result_type='', extension=extension), overwrite=overwrite)
    else:
        compute_lambdas_by_expected_ac_ht(freq_lower=freq_lower, var_filter=var_filter, random_phenos=random_phenos).write(get_ukb_exomes_sumstat_path(subdir=lambda_folder, dataset=f'{rp}lambda_expectedAC{filter}', result_type='var', extension=extension), overwrite=overwrite)


def compute_ukb_pheno_moments_ht(pheno_sex='both_sexes', phenocode: list = None):
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
    pheno = pheno.annotate_cols(mean=hl.agg.mean(pheno[pheno_sex]), 
                                std=hl.agg.stats(pheno[pheno_sex]).stdev)
    pheno = pheno.annotate_cols(variance=pheno.std ** 2, 
                                skewness=hl.agg.mean((pheno[pheno_sex] - pheno.mean) ** 3) / (pheno.std ** 3), 
                                kurtosis=hl.agg.mean((pheno[pheno_sex] - pheno.mean) ** 4) / (pheno.std ** 4))
    pheno = pheno.select_cols('mean', 'variance', 'skewness', 'kurtosis').cols()
    return pheno


def get_caf_info_ht(tranche: str = CURRENT_TRANCHE):
    """
    Generate cumulative allele frequency Table
    :param str tranche: Tranche of data to use
    :return: Cumulative allele frequency Hail Table
    :rtype: Table
    """
    vep = hl.read_table(var_annotations_ht_path('vep', *TRANCHE_DATA[tranche]))
    vep = process_consequences(vep)
    vep = vep.explode(vep.vep.worst_csq_by_gene_canonical)
    annotation = annotation_case_builder(vep.vep.worst_csq_by_gene_canonical)
    vep = vep.annotate(annotation=hl.if_else(hl.literal({'missense', 'LC'}).contains(annotation), 'missense|LC', annotation))
    vep = vep.filter(hl.is_defined(vep.annotation) & (vep.annotation != 'non-coding'))
    vep = vep.select(vep.vep.worst_csq_by_gene_canonical.gene_id, vep.annotation)

    freq = hl.read_table(release_ht_path(*TRANCHE_DATA[tranche]))
    var_af = vep.annotate(AF=freq[vep.key].freq[0].AF)
    sum_af = var_af.group_by(var_af.annotation, var_af.gene_id).aggregate(CAF=hl.agg.sum(var_af.AF))
    return sum_af


def get_sig_cnt_mt(result_type: str = 'gene', test_type: str = 'skato', phenos_to_keep: hl.Table = None, genes_to_keep: hl.Table = None, var_min_freq: float=None, tranche: str = CURRENT_TRANCHE):
    """
    Generate the number of significant association for each phenotype and each gene/variant
    :param str result_type: Type of association test result, `gene` or `variant`
    :param str test_type: Type of gene-level association test result, `skato` or `burden`
    :param Table phenos_to_keep: Subset of phenotypes to keep, keyed by ['trait_type', 'phenocode', 'pheno_sex', 'coding', 'modifier']
    :param Table genes_to_keep: Subset of genes to keep, keyed by ['gene_id', 'gene_symbol', 'annotation']
    :param float var_min_freq: Minimum allele frequency used for variant filtering
    :param str tranche: Tranche of data to use
    :return: Hail MatrixTable with phenotype-level and gene/variant-level association count
    :rtype: MatrixTable
    """
    mt = hl.read_matrix_table(get_results_mt_path(result_type, tranche=tranche))
    mt = mt.annotate_cols(trait_type2=hl.if_else(mt.trait_type == 'icd_first_occurrence', 'icd10', mt.trait_type))
    pvalue = P_VALUE_FIELDS[test_type.lower()]
    if phenos_to_keep is not None:
        mt = mt.filter_cols(hl.is_defined(phenos_to_keep.index(mt.col_key)))
    if result_type == 'gene':
        af = get_caf_info_ht(tranche=tranche)
        mt = mt.annotate_rows(CAF=af[mt.annotation, mt.gene_id]['CAF'])
        if genes_to_keep is not None:
            mt = mt.filter_rows(hl.is_defined(genes_to_keep.index(mt.row_key)))
        mt = mt.annotate_rows(sig_pheno_cnt=hl.agg.group_by(mt.trait_type2, hl.agg.count_where(mt[pvalue] < EMPIRICAL_P_THRESHOLDS[test_type.lower()])),
                              all_sig_pheno_cnt=hl.agg.filter(hl.is_defined(mt[pvalue]), hl.agg.count_where(mt[pvalue] < EMPIRICAL_P_THRESHOLDS[test_type.lower()])))
        mt = mt.annotate_rows(**{f'{trait_type2}_sig_pheno_cnt_{test_type.lower()}': mt.sig_pheno_cnt.get(trait_type2) for trait_type2 in TRAIT_TYPES})

        mt = mt.annotate_cols(sig_gene_cnt=hl.agg.group_by(mt.annotation, hl.agg.count_where(mt[pvalue] < EMPIRICAL_P_THRESHOLDS[test_type.lower()])),
                              all_sig_gene_cnt=hl.agg.filter(hl.is_defined(mt[pvalue]), hl.agg.count_where(mt[pvalue]< EMPIRICAL_P_THRESHOLDS[test_type.lower()])))
        mt = mt.annotate_cols(**{f'{annotation}_sig_gene_cnt_{test_type.lower()}': mt.sig_gene_cnt.get(annotation) for annotation in ANNOTATIONS})
    else:
        mt = mt.annotate_rows(annotation=hl.if_else(hl.literal({'missense', 'LC'}).contains(mt.annotation), 'missense|LC', mt.annotation))
        if var_min_freq is not None:
            mt = mt.filter_rows((mt.AF > var_min_freq) & hl.is_defined(mt.annotation))
            mt = mt.annotate_rows(all_sig_pheno_cnt=hl.agg.filter(mt.SE != 0, hl.agg.count_where(mt.Pvalue < EMPIRICAL_P_THRESHOLDS['variant'])),
                                  sig_pheno_cnt=hl.agg.group_by(mt.trait_type2, hl.agg.filter(mt.SE != 0, hl.agg.count_where(mt.Pvalue < EMPIRICAL_P_THRESHOLDS['variant']))),)
            mt = mt.annotate_rows(**{f'{trait_type}_sig_pheno_cnt': mt.sig_pheno_cnt.get(trait_type) for trait_type in TRAIT_TYPES for test in TESTS},)
            mt = mt.annotate_cols(all_sig_var_cnt=hl.agg.filter(mt.SE != 0, hl.agg.count_where(mt.Pvalue < EMPIRICAL_P_THRESHOLDS['variant'])),
                                  sig_var_cnt=hl.agg.group_by(mt.annotation, hl.agg.filter(mt.SE != 0, hl.agg.count_where(mt.Pvalue < EMPIRICAL_P_THRESHOLDS['variant']))))
            mt = mt.annotate_cols(**{f'{annotation}_sig_var_cnt': mt.sig_var_cnt.get(annotation) for annotation in ANNOTATIONS},)
        else:
            mt = mt.annotate_rows(all_sig_pheno_cnt=hl.agg.sum(mt.Pvalue < EMPIRICAL_P_THRESHOLDS['variant']),
                                  sig_pheno_cnt=hl.agg.group_by(mt.trait_type2, hl.agg.count_where(mt.Pvalue < EMPIRICAL_P_THRESHOLDS['variant'])),)
            mt = mt.annotate_rows(**{f'{trait_type}_sig_pheno_cnt': mt.sig_pheno_cnt.get(trait_type) for trait_type in TRAIT_TYPES for test in TESTS},)
            mt = mt.annotate_cols(all_sig_var_cnt=hl.agg.sum(mt.Pvalue < EMPIRICAL_P_THRESHOLDS['variant']),
                                  sig_var_cnt=hl.agg.group_by(mt.annotation, hl.agg.count_where(mt.Pvalue < EMPIRICAL_P_THRESHOLDS['variant'])))
            mt = mt.annotate_cols(**{f'{annotation}_sig_var_cnt': mt.sig_var_cnt.get(annotation) for annotation in ANNOTATIONS},)
    return mt

def compare_gene_var_sig_cnt_mt(test_type: str = 'skato', phenos_to_keep: hl.Table = None, genes_to_keep: hl.Table = None, var_min_freq: float = None, tranche: str = CURRENT_TRANCHE):
    """
    Compare the number of significant association for genes with variants within corresponding genes
    :param str test_type: Type of gene-level association test result, `skato` or `burden`
    :param Table phenos_to_keep: Subset of phenotypes to keep, keyed by ['trait_type', 'phenocode', 'pheno_sex', 'coding', 'modifier']
    :param Table genes_to_keep: Subset of genes to keep, keyed by ['gene_id', 'gene_symbol', 'annotation']
    :param float var_min_freq: Minimum allele frequency used for variant filtering
    :param str tranche: Tranche of data to use
    :return: Hail MatrixTable with variant-level association count grouped by genes at phenotype-level and gene-level
    :rtype: MatrixTable
    """
    var = hl.read_matrix_table(get_results_mt_path('variant', tranche=tranche))
    mt = hl.read_matrix_table(get_results_mt_path(tranche=tranche))
    if phenos_to_keep is not None:
        mt = mt.filter_cols(hl.is_defined(phenos_to_keep.index(mt.col_key)))
        var = var.filter_cols(hl.is_defined(phenos_to_keep.index(var.col_key)))
    if genes_to_keep is not None:
        mt = mt.filter_rows(hl.is_defined(genes_to_keep.index(mt.row_key)))
    vep = hl.read_table(var_annotations_ht_path('vep', *TRANCHE_DATA[tranche]))
    vep = process_consequences(vep)
    var = var.annotate_rows(gene_id=vep[var.row_key].vep.worst_csq_by_gene_canonical.gene_id)
    var = var.explode_rows(var.gene_id)
    var = var.annotate_rows(annotation=hl.if_else(hl.literal({'missense', 'LC'}).contains(var.annotation), 'missense|LC', var.annotation), )
    if var_min_freq is not None:
        var = var.filter_rows((var.AF>var_min_freq) & hl.is_defined(var.annotation))
        var = var.group_rows_by('gene_id', 'gene', 'annotation').aggregate(var_cnt=hl.agg.filter(var.SE != 0, hl.agg.count_where(var.Pvalue < EMPIRICAL_P_THRESHOLDS['variant'])))
    else:
        var = var.group_rows_by('gene_id', 'gene', 'annotation').aggregate(var_cnt=hl.agg.count_where(var.Pvalue < EMPIRICAL_P_THRESHOLDS['variant']))
    pvalue = P_VALUE_FIELDS[test_type.lower()]
    mt = mt.annotate_entries(var_cnt=var[mt.row_key, mt.col_key]['var_cnt'])
    mt = mt.annotate_cols(pheno_gene_var_sig_cnt=hl.agg.count_where(hl.is_defined(mt.var_cnt) & (mt.var_cnt > 0) &
                                                                  hl.is_defined(mt[pvalue]) & (mt[pvalue] < EMPIRICAL_P_THRESHOLDS[test_type.lower()])),
                          pheno_var_sig_cnt=hl.agg.count_where((hl.is_defined(mt.var_cnt) & (mt.var_cnt > 0)) &
                                                             (~(hl.is_defined(mt[pvalue]) & (mt[pvalue] < EMPIRICAL_P_THRESHOLDS[test_type.lower()])))),
                          pheno_gene_sig_cnt=hl.agg.count_where((hl.is_defined(mt[pvalue]) & (mt[pvalue] < EMPIRICAL_P_THRESHOLDS[test_type.lower()])) &
                                                              (~(hl.is_defined(mt.var_cnt) & (mt.var_cnt > 0)))),
                          pheno_none_sig_cnt=hl.agg.count_where(~((hl.is_defined(mt[pvalue]) & (mt[pvalue] < EMPIRICAL_P_THRESHOLDS[test_type.lower()])) |
                                                                (hl.is_defined(mt.var_cnt) & (mt.var_cnt > 0)))),
                          trait_type2=hl.if_else(mt.trait_type == 'icd_first_occurrence', 'icd10', mt.trait_type))

    mt = mt.annotate_rows(gene_var_sig_cnt=hl.agg.count_where(hl.is_defined(mt.var_cnt) & (mt.var_cnt > 0) &
                                                                  hl.is_defined(mt[pvalue]) & (mt[pvalue] < EMPIRICAL_P_THRESHOLDS[test_type.lower()])),
                          var_sig_cnt=hl.agg.count_where((hl.is_defined(mt.var_cnt) & (mt.var_cnt > 0)) &
                                                             (~(hl.is_defined(mt[pvalue]) & (mt[pvalue] < EMPIRICAL_P_THRESHOLDS[test_type.lower()])))),
                          gene_sig_cnt=hl.agg.count_where((hl.is_defined(mt[pvalue]) & (mt[pvalue] < EMPIRICAL_P_THRESHOLDS[test_type.lower()])) &
                                                              (~(hl.is_defined(mt.var_cnt) & (mt.var_cnt > 0)))),
                          none_sig_cnt=hl.agg.count_where(((hl.is_missing(mt[pvalue]) | (mt[pvalue] >= EMPIRICAL_P_THRESHOLDS[test_type.lower()])) &
                                                           (hl.is_missing(mt.var_cnt) | (mt.var_cnt == 0)))))
    return mt

def compute_mean_coverage_ht(tranche: str = CURRENT_TRANCHE):
    """
    Generate mean coverage Table
    :param str tranche: Tranche of data to use
    :return: Mean coverage Hail Table
    :rtype: Table
    """
    int_auto = hl.read_table(interval_qc_path(data_source='broad', freeze=7, chrom='autosomes'))
    int_sex = hl.read_table(interval_qc_path(data_source='broad', freeze=7, chrom='sex_chr'))
    int_sex = int_sex.select(target_mean_dp=int_sex.target_mean_dp['XX'], 
                            target_pct_gt_10x=int_sex.target_pct_gt_10x['XX'], 
                            target_pct_gt_20x=int_sex.target_pct_gt_20x['XX'], 
                            pct_samples_10x=int_sex.pct_samples_10x['XX'], 
                            pct_samples_20x=int_sex.pct_samples_20x['XX'])
    int_full = int_auto.union(int_sex)

    var = hl.read_matrix_table(get_results_mt_path('variant', tranche=tranche)).rows()
    vep = hl.read_table(var_annotations_ht_path('vep', *TRANCHE_DATA[tranche]))
    vep = process_consequences(vep)
    var = var.annotate(csq=vep[var.key].vep.worst_csq_by_gene_canonical)
    var = var.explode(var.csq)
    var = var.annotate(gene_id = var.csq.gene_id, 
                       gene_symbol = var.csq.gene_symbol, 
                       annotation = annotation_case_builder(var.csq), 
                       coverage=int_full[var.locus].target_mean_dp)
    var = var.annotate(annotation=hl.if_else(hl.literal({'missense', 'LC'}).contains(var.annotation), 'missense|LC', var.annotation), )
    mean_coverage = var.group_by('gene_id', 'gene_symbol', 'annotation').aggregate(mean_coverage=hl.agg.mean(var.coverage))
    return mean_coverage

def more_cases_tie_breaker(l, r):
    """
    Tie breaker function prefering phenotypes with more cases used in hl.maximal_independent_ser()
    :param l: one phenotype
    :param r: the other phenotype
    :return: integer indicating which phenotype is preferred
    :rtype: int
    """
    return (hl.case()\
        .when(l.n_cases_both_sexes > r.n_cases_both_sexes, -1)\
        .when(l.n_cases_both_sexes == r.n_cases_both_sexes, 0)\
        .when(l.n_cases_both_sexes < r.n_cases_both_sexes, 1)\
        .or_missing())

def get_lambda_filtered_ht(lambda_ht: hl.Table = None, lambda_name: str = 'all_lambda_gc_skato', lower:float = None, upper: float = None):
    """
    Filter lambda GC to specific range
    :param Table lambda_ht: Hail table with original lambda GC information
    :param str lambda_name: name of the lambda GC column to apply filters on
    :param float lower: Lower bound of lambda GC
    :param float upper: Upper bound of lambda GC
    :return: Hail Table of filtered lambda GC
    :rtype: Table
    """
    if lower is not None:
        lambda_ht = lambda_ht.filter(lambda_ht[lambda_name] > lower)
    if upper is not None:
        lambda_ht = lambda_ht.filter(lambda_ht[lambda_name] < upper)
    return lambda_ht

def get_corr_phenos_ht(r_2: float = None, tie_breaker = None, tranche: str = CURRENT_TRANCHE):
    """
    Generate a subset of related phenotypes to remove
    :param float r_2: Correlation-square cutoff point, phenotype with correlation above which are included to be removed
    :param tie_breaker: Tie-breaker function to select phenotype based on certain preference
    :param str tranche: Tranche of data to use
    :return: Hail Table of correlated phenotypes to remove
    :rtype: Table
    """
    pheno_mt = get_ukb_pheno_mt()
    ht = hl.read_matrix_table(get_results_mt_path(tranche)).cols()
    pheno_mt = pheno_mt.filter_cols(hl.is_defined(ht[pheno_mt.col_key]))
    corr = make_pairwise_ht(pheno_mt, pheno_field=pheno_mt.both_sexes, correlation=True)
    related = corr.filter((corr.entry ** 2 >= r_2) & (corr.i != corr.j))
    pheno_to_remove = hl.maximal_independent_set(related.i_data, related.j_data, keep=False, tie_breaker=tie_breaker)
    return pheno_to_remove

def export_ht_to_txt(path_to_ht: str, output_filename: str):
    """
    Convert Hail Table to txt.bgz
    :param str path_to_ht: Path to input hail table
    :param str output_filename: Output file name
    """
    ht = hl.read_table(path_to_ht)
    ht.export(f'{sumstats_folder}/{output_filename}.txt.bgz')

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
        pheno_to_remove = hl.maximal_independent_set(related.i_data, related.j_data, keep=False)
        l.append(pheno_to_remove.count())
    return l

def get_icd_min_p_ht(result_type: str = 'gene', test_type: str = 'skato', phenos_to_keep: hl.Table = None, genes_to_keep: hl.Table = None, var_min_freq: float=None, tranche: str = CURRENT_TRANCHE):
    """
    Generate min pvalue for each icd10 group from variant-level or gene-level test result
    :param str result_type: Type of association test result, `gene` or `variant`
    :param str test_type: Type of gene-level association test result, `skato` or `burden`
    :param Table phenos_to_keep: Subset of phenotypes to keep, keyed by ['trait_type', 'phenocode', 'pheno_sex', 'coding', 'modifier']
    :param Table genes_to_keep: Subset of genes to keep, keyed by ['gene_id', 'gene_symbol', 'annotation']
    :param float var_min_freq: Minimum allele frequency used for variant filtering
    :param str tranche: Tranche of data to use
    :return: Hail Table of min p-values for each icd10 group
    :rtype: Table
    """
    mt = hl.read_matrix_table(get_results_mt_path(result_type, tranche=tranche))
    pvalue = P_VALUE_FIELDS[test_type.lower()]
    if phenos_to_keep is not None:
        mt = mt.filter_cols(hl.is_defined(phenos_to_keep.index(mt.col_key)))
    icd = mt.filter_cols(mt.trait_type == 'icd_first_occurrence')
    icd = icd.annotate_cols(icd10_group = icd.description.split('first', 2)[0][5])
    ICD10_GROUPS = icd.aggregate_cols(hl.agg.collect_as_set(icd.icd10_group))
    if result_type == 'gene':
        if genes_to_keep is not None:
            icd = icd.filter_rows(hl.is_defined(genes_to_keep.index(icd.row_key)))
        icd = icd.annotate_rows(min_p = hl.agg.group_by(icd.icd10_group, hl.agg.min(icd[pvalue])))
    else:
        if var_min_freq is not None:
            icd = icd.filter_rows((icd.AF > var_min_freq) & hl.is_defined(icd.annotation))
        icd = icd.annotate_rows(min_p = hl.agg.group_by(icd.icd10_group, hl.agg.filter(icd.SE != 0, hl.agg.min(icd.Pvalue))))
    icd = icd.annotate_rows(**{f'{icd_group}_min_p': icd.min_p.get(icd_group) for icd_group in ICD10_GROUPS})
    return icd.rows()

def get_pheno_pvalue_ht(result_type: str = 'gene', test_type: str = 'skato', coding: list = None, phenos_to_keep: hl.Table = None, freq_lower: float = None, freq_upper: float = None,
                        n_var_min: int = 2, coverage_min: int = 20, random_phenos: bool = False, tranche: str = CURRENT_TRANCHE):
    """
    Generate Pvalue table from variant-level or gene-level test result
    :param str result_type: Type of association test result, `gene` or `variant`
    :param str test_type: Type of gene-level association test result, `skato` or `burden`
    :param list coding: List of coding used for phenotype filtering
    :param Table phenos_to_keep: Subset of phenotypes to keep, keyed by ['trait_type', 'phenocode', 'pheno_sex', 'coding', 'modifier']
    :param float freq_lower: Lower bound of allele frequency/cumulative allele frequency used for variant/gene filtering
    :param float freq_upper: Upper bound of allele frequency/cumulative allele frequency used for variant/gene filtering
    :param int n_var_min: Minimum number of variants included in a gene-based test used for gene filtering
    :param int coverage_min: Minimum coverage used for gene filtering
    :param bool random_pheno: Whether to use randomly-generated phenotypes
    :param str tranche: Tranche of data to use
    :return: Hail Table of p-value
    :rtype: Table
    """
    mt = hl.read_matrix_table(get_results_mt_path(result_type, tranche, random_phenos=random_phenos))
    if coding is not None:
        mt = mt.filter_cols(hl.literal(coding).contains(mt.coding))
    if phenos_to_keep is not None:
        mt = mt.filter_cols(hl.is_defined(phenos_to_keep.index(mt.col_key)))
    if result_type == 'gene':  # Gene Info
        mt = annotate_additional_info_mt(result_mt=mt, result_type=result_type, tranche=tranche)
        if freq_lower is not None:
            mt = mt.filter_rows(mt.CAF > freq_lower)
        if freq_upper is not None:
            mt = mt.filter_rows(mt.CAF <= freq_upper)
        if n_var_min is not None:
            mt = mt.filter_rows(mt.total_variants >= n_var_min)
        if coverage_min is not None:
            mt = mt.filter_rows(mt.mean_coverage > coverage_min)
        mt = mt.annotate_rows(CAF_range=f'CAF:({freq_lower}, {freq_upper}]')
        ht = mt.entries()
        pvalue = P_VALUE_FIELDS[test_type.lower()]
        ht = ht.select(pvalue, 'CAF_range')
    else:  # Variant Info
        if freq_lower is not None:
            mt = mt.filter_rows(mt.AF > af_lower)
        if freq_upper is not None:
            mt = mt.filter_rows(mt.AF <= af_upper)
        mt = mt.filter_rows(hl.is_defined(af.index(mt['locus'], mt['alleles'])))
        mt = mt.annotate_rows(annotation=hl.if_else(hl.literal({'missense', 'LC'}).contains(mt.annotation), 'missense|LC', mt.annotation),
                              AF_range=f'AF:({freq_lower}, {freq_upper}]')
        ht = mt.entries()
        ht = ht.select('Pvalue', 'AF_range')
    return ht

def get_pvalue_by_freq_interval_ht(result_type: str = 'gene', test_type: str = 'skato', freq_breaks: list = [0.0001, 0.001, 0.01, 0.1], phenos_to_keep: hl.Table = None,
                                   n_var_min: int = None, coverage_min: int = None, random_phenos: bool = False, tranche: str = CURRENT_TRANCHE):
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
    ht = get_pheno_pvalue_ht(result_type = result_type, test_type = test_type, phenos_to_keep = phenos_to_keep,freq_lower=freq_breaks[-1], n_var_min=n_var_min, coverage_min=coverage_min, random_phenos=random_phenos, tranche=tranche)
    for i in list(range(0, len(freq_breaks)-1)):
        sub_ht = get_pheno_pvalue_ht(result_type = result_type, test_type = test_type, phenos_to_keep = phenos_to_keep, freq_lower=freq_breaks[i], freq_upper=freq_breaks[i+1],
                                     n_var_min=n_var_min, coverage_min=coverage_min, random_phenos=random_phenos, tranche=tranche)
        ht = ht.union(sub_ht)
    return ht

def annotate_clinvar_pathogenicity_ht():
    """
    Re-define pathogenicity group for clinvar variants
    :return: Hail Table with new definitions
    :rtype: Table
    """
    clinvar_ht = clinvar.ht()
    clinvar_ht = clinvar_ht.filter(~clinvar_ht.info.CLNREVSTAT[0].contains('no_assertion'))  # remove clinvar zero star variants
    clinvar_ht = clinvar_ht.annotate(pathogenicity=hl.case() \
                                     .when(hl.literal({'Pathogenic', 'Likely_pathogenic', 'Pathogenic/Likely_pathogenic'}).contains(clinvar_ht.info.CLNSIG), 'P/LP') \
                                     .when(hl.literal({'Uncertain_significance'}).contains(clinvar_ht.info.CLNSIG), 'VUS') \
                                     .when(hl.literal({'Benign', 'Likely_benign', 'Benign/Likely_benign'}).contains(clinvar_ht.info.CLNSIG), 'B/LB') \
                                     .or_missing())
    return clinvar_ht

def add_liftover_rg37_to_rg38_ht(ht: hl.Table):
    """
    Add liftover to Hail Table from rg37 to rg38
    :param Table ht: Hail Table to add liftover on
    :return: Hail Table
    :rtype: hl.Table
    """
    rg37 = hl.get_reference('GRCh37')
    rg38 = hl.get_reference('GRCh38')
    if not rg37.has_liftover("GRCh38"):
        rg37.add_liftover('gs://hail-common/references/grch37_to_grch38.over.chain.gz', rg38)
    ht = ht.annotate(new_locus=hl.liftover(ht.locus, 'GRCh38'))
    ht = ht.filter(hl.is_defined(ht.new_locus))
    ht = ht.key_by(locus=ht.new_locus, alleles=ht.alleles)
    return ht

def annotate_additional_info_mt(result_mt: hl.MatrixTable, result_type: str = 'gene', tranche: str = CURRENT_TRANCHE):
    """
    Attach additional information to the original gene/variant-level sumstats MatrixTable
    :param MatrixTable result_mt: sumstats MatrixTable to annotate informtaion to
    :param str result_type: Type of association test result, `gene` or `variant`
    :param str tranche: Tranche of data to use
    :return: Hail MatrixTable with updated information
    :rtype: hl.MatrixTable
    """
    if result_type == 'gene':
        af = hl.read_table(get_util_info_path('caf', tranche))
        coverage_ht = hl.read_table(get_util_info_path('coverage', tranche))
        mt = result_mt.annotate_rows(CAF=af[result_mt.annotation, result_mt.gene_id]['CAF'],
                                     mean_coverage=coverage_ht[result_mt.row_key].mean_coverage)
    else:
        clinvar_ht = annotate_clinvar_pathogenicity_ht()

        ht = hl.read_table('gs://gcp-public-data--gnomad/papers/2019-tx-annotation/pre_computed/all.possible.snvs.tx_annotated.021520.ht')
        ht = ht.explode(ht.tx_annotation)
        ht = add_liftover_rg37_to_rg38_ht(ht)

        vep = hl.read_table(var_annotations_ht_path('vep', *TRANCHE_DATA[tranche]))
        vep = process_consequences(vep)
        vep = vep.explode(vep.vep.worst_csq_by_gene_canonical)
        vep = vep.annotate(tx_annotation_csq=ht[vep.key].tx_annotation.csq,
                           mean_proportion=ht[vep.key].tx_annotation.mean_proportion)
        mt = result_mt.annotate_rows(pathogenicity=clinvar_ht[result_mt.locus, result_mt.alleles]['pathogenicity'],
                                     annotation=hl.if_else(hl.literal({'missense', 'LC'}).contains(result_mt.annotation),'missense|LC', result_mt.annotation),
                                     polyphen2=vep[result_mt.row_key].vep.worst_csq_by_gene_canonical.polyphen_prediction,
                                     amino_acids=vep[result_mt.row_key].vep.worst_csq_by_gene_canonical.amino_acids,
                                     lof=vep[result_mt.row_key].vep.worst_csq_by_gene_canonical.lof,
                                     most_severe_consequence=vep[result_mt.row_key].vep.worst_csq_by_gene_canonical.most_severe_consequence,
                                     tx_annotation_csq=vep[result_mt.row_key].tx_annotation_csq,
                                     mean_proportion_expressed=vep[result_mt.row_key].mean_proportion)
    return mt