from ukbb_qc.resources.sample_qc import interval_qc_path
from ukbb_qc.resources.variant_qc import var_annotations_ht_path
from ukbb_qc.resources.basics import release_ht_path
from ukb_common.utils.annotations import annotation_case_builder
from gnomad.utils.vep import process_consequences
from ukb_exomes.resources import *

ANNOTATIONS = ('pLoF', 'missense|LC', 'synonymous')
TESTS = ('skato', 'skat', 'burden')
TRAIT_TYPES = ('continuous', 'categorical')

def compute_lambda_gc_ht(result_type: str = 'gene', by_annotation: bool = False, by_gene: bool = False, tranche: str = CURRENT_TRANCHE,
                         af_lower: float = None, af_upper: float = None, n_var_min: int = None, random_phenos: bool = False):
    mt = hl.read_matrix_table(get_results_mt_path(result_type, tranche, random_phenos=random_phenos))
    if result_type == 'gene':  # Gene Info
        af = get_caf_info_ht(tranche)
        mt = mt.annotate_rows(CAF=af[mt.annotation, mt.gene_id]['CAF'])
        if af_lower is not None:
            mt = mt.filter_rows(mt.CAF > af_lower)
        if af_upper is not None:
            mt = mt.filter_rows(mt.CAF <= af_upper)
        if n_var_min is not None:
            mt = mt.filter_rows(mt.total_variants >= n_var_min)
        if by_annotation:
            mt = mt.select_cols('n_cases', 
                                lambda_gc_skato=hl.agg.group_by(mt.annotation, hl.methods.statgen._lambda_gc_agg(mt.Pvalue)), 
                                lambda_gc_skat=hl.agg.group_by(mt.annotation, hl.methods.statgen._lambda_gc_agg(mt.Pvalue_SKAT)), 
                                lambda_gc_burden=hl.agg.group_by(mt.annotation, hl.methods.statgen._lambda_gc_agg(mt.Pvalue_Burden)), 
                                CAF_range=f'CAF:({af_lower}, {af_upper}]')
            mt = mt.annotate_cols(**{f'{annotation}_lambda_gc_{test}': mt[f'lambda_gc_{test}'][annotation] for annotation in ANNOTATIONS for test in TESTS})
            ht = mt.cols()
        elif by_gene:
            mt = mt.annotate_cols(trait_type2=hl.if_else(mt.trait_type == 'continuous', mt.trait_type, 'categorical'))
            mt = mt.select_rows('CAF',
                                lambda_gc_skato=hl.agg.group_by(mt.trait_type2, hl.methods.statgen._lambda_gc_agg(mt.Pvalue)),
                                lambda_gc_skat=hl.agg.group_by(mt.trait_type2, hl.methods.statgen._lambda_gc_agg(mt.Pvalue_SKAT)),
                                lambda_gc_burden=hl.agg.group_by(mt.trait_type2, hl.methods.statgen._lambda_gc_agg(mt.Pvalue_Burden)), )
            mt = mt.annotate_rows(**{f'{trait_type}_lambda_gc_{test}': mt[f'lambda_gc_{test}'][trait_type] for trait_type in TRAIT_TYPES for test in TESTS},
                                  all_lambda_gc_skato=hl.agg.filter(hl.is_defined(mt.Pvalue), hl.methods.statgen._lambda_gc_agg(mt.Pvalue)),
                                  all_lambda_gc_skat=hl.agg.filter(hl.is_defined(mt.Pvalue_SKAT), hl.methods.statgen._lambda_gc_agg(mt.Pvalue_SKAT)),
                                  all_lambda_gc_burden=hl.agg.filter(hl.is_defined(mt.Pvalue_Burden), hl.methods.statgen._lambda_gc_agg(mt.Pvalue_Burden)), )
            ht = mt.rows()
        else:
            mt = mt.select_cols('n_cases', 
                                lambda_gc_skato=hl.agg.filter(hl.is_defined(mt.Pvalue), hl.methods.statgen._lambda_gc_agg(mt.Pvalue)), 
                                lambda_gc_skat=hl.agg.filter(hl.is_defined(mt.Pvalue_SKAT), hl.methods.statgen._lambda_gc_agg(mt.Pvalue_SKAT)), 
                                lambda_gc_burden=hl.agg.filter(hl.is_defined(mt.Pvalue_Burden), hl.methods.statgen._lambda_gc_agg(mt.Pvalue_Burden)), 
                                CAF_range=f'CAF:({af_lower}, {af_upper}]')
            ht = mt.cols()
    else:  # Variant Info
        if af_lower is not None:
            mt = mt.filter_rows(mt.AF > af_lower)
        if af_upper is not None:
            mt = mt.filter_rows(mt.AF <= af_upper)
        if by_annotation:
            lambda_gc = hl.agg.group_by(mt.annotation, hl.methods.statgen._lambda_gc_agg(mt.Pvalue))
            mt = mt.select_cols('n_cases', lambda_gc=lambda_gc, AF_range=f'AF:({af_lower}, {af_upper}]')
            mt = mt.annotate_rows(annotation=hl.if_else(hl.literal({'missense', 'LC'}).contains(mt.annotation), 'missense|LC', mt.annotation))
            mt = mt.annotate_cols(**{f'{annotation}_lambda_gc_variant': mt.lambda_gc[annotation] for annotation in ('pLoF', 'missense|LC', 'synonymous')}, )
        else:
            lambda_gc = hl.agg.filter(hl.is_defined(mt.Pvalue), hl.methods.statgen._lambda_gc_agg(mt.Pvalue))
            mt = mt.select_cols('n_cases', lambda_gc=lambda_gc, AF_range=f'AF:({af_lower}, {af_upper}]')
        ht = mt.cols()
    return ht

def compute_ukb_pheno_moments_ht(pheno_sex='both_sexes', phenocode: list = None):
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
    vep = hl.read_table(var_annotations_ht_path('vep', *TRANCHE_DATA[tranche]))
    vep = process_consequences(vep)
    vep = vep.explode(vep.vep.worst_csq_by_gene_canonical)
    annotation = annotation_case_builder(vep.vep.worst_csq_by_gene_canonical)
    vep = vep.annotate(annotation=hl.if_else(hl.literal({'missense', 'LC'}).contains(annotation), 'missense|LC', annotation))
    vep = vep.filter((hl.is_defined(vep.annotation)) & (vep.annotation != 'non-coding'))
    vep = vep.select(vep.vep.worst_csq_by_gene_canonical.gene_id, vep.annotation)

    freq = hl.read_table(release_ht_path(*TRANCHE_DATA[tranche]))
    var_af = vep.annotate(AF=freq[vep.key].freq[0].AF)
    sum_af = var_af.group_by(var_af.annotation, var_af.gene_id).aggregate(CAF=hl.agg.sum(var_af.AF))
    return sum_af


def get_sig_cnt_mt(result_type: str = 'gene', phenos_to_keep: hl.Table = None, genes_to_keep: hl.Table = None, level: float = 1e-6, tranche: str = CURRENT_TRANCHE):
    mt = hl.read_matrix_table(get_results_mt_path(result_type, tranche=tranche))
    if phenos_to_keep is not None:
        mt = mt.filter_cols(hl.is_defined(phenos_to_keep.index(mt.col_key)))
    if result_type == 'gene':
        af = get_caf_info_ht(tranche=tranche)
        mt = mt.annotate_rows(CAF=af[mt.annotation, mt.gene_id]['CAF'])
        if genes_to_keep is not None:
            mt = mt.filter_rows(hl.is_defined(genes_to_keep.index(mt.row_key)))
        # Count Significant Hits per Gene for each Test
        mt = mt.annotate_rows(sig_pheno_cnt_skato=hl.agg.sum(mt.Pvalue < level), 
                              sig_pheno_cnt_skat=hl.agg.sum(mt.Pvalue_SKAT < level), 
                              sig_pheno_cnt_burden=hl.agg.sum(mt.Pvalue_Burden < level))
        # Count Significant Hits per Phenotype for each Test
        mt = mt.annotate_cols(sig_gene_cnt_skato=hl.agg.sum(mt.Pvalue < level), 
                              sig_gene_cnt_skat=hl.agg.sum(mt.Pvalue_SKAT < level), 
                              sig_gene_cnt_burden=hl.agg.sum(mt.Pvalue_Burden < level))
    else:
        mt = mt.annotate_rows(annotation=hl.if_else(hl.literal({'missense', 'LC'}).contains(mt.annotation), 'missense|LC', mt.annotation),
                              sig_pheno_cnt=hl.agg.sum(mt.Pvalue < level))
        # Count Significant Hits per Phenotype for each Test
        mt = mt.annotate_cols(sig_var_cnt=hl.agg.sum(mt.Pvalue < level))
    return mt

def get_sig_cnt_annt_ht(result_type: str = 'gene', phenos_to_keep: hl.Table = None, genes_to_keep: hl.Table = None, level: float = 1e-6, tranche: str = CURRENT_TRANCHE):
    mt = hl.read_matrix_table(get_results_mt_path(result_type, tranche))
    if phenos_to_keep is not None:
        mt = mt.filter_cols(hl.is_defined(phenos_to_keep.index(mt.col_key)))
    if result_type == 'gene':
        if genes_to_keep is not None:
            mt = mt.filter_rows(hl.is_defined(genes_to_keep.index(mt.row_key)))
        mt = mt.annotate_cols(sig_cnt_skato=hl.agg.group_by(mt.annotation, hl.agg.count_where(mt.Pvalue < level)), 
                              sig_cnt_skat=hl.agg.group_by(mt.annotation, hl.agg.count_where(mt.Pvalue_SKAT < level)), 
                              sig_cnt_burden=hl.agg.group_by(mt.annotation, hl.agg.count_where(mt.Pvalue_Burden < level)), )
        ht = mt.cols().select('n_cases', 'description', 'sig_cnt_skato', 'sig_cnt_skat', 'sig_cnt_burden')
        ht = ht.annotate(**{f'{annotation}_sig_cnt_{test}': ht[f'sig_cnt_{test}'][annotation] for annotation in ANNOTATIONS for test in TESTS})
    else:
        mt = mt.annotate_cols(sig_cnt=hl.agg.group_by(mt.annotation, hl.agg.count_where(mt.Pvalue < level)), )
        mt = mt.annotate_rows(annotation=hl.if_else(hl.literal({'missense', 'LC'}).contains(mt.annotation), 'missense|LC', mt.annotation))
        ht = mt.cols().select('n_cases', 'description', 'sig_cnt')
        ht = ht.annotate(**{f'{annotation}_sig_cnt': ht.sig_cnt[annotation] for annotation in ('pLoF', 'missense|LC', 'synonymous')}, )
    return ht

def compare_gene_var_sig_cnt_mt(test_type: str = 'skato', phenos_to_keep: hl.Table = None, genes_to_keep: hl.Table = None, level: float = 1e-6, tranche: str = CURRENT_TRANCHE):
    var = hl.read_matrix_table(get_results_mt_path('variant', tranche=tranche))
    mt = hl.read_matrix_table(get_results_mt_path(tranche=tranche))
    if phenos_to_keep is not None:
        mt = mt.filter_cols(hl.is_defined(phenos_to_keep.index(mt.col_key)))
        var = var.filter_cols(hl.is_defined(phenos_to_keep.index(var.col_key)))
    vep = hl.read_table(var_annotations_ht_path('vep', *TRANCHE_DATA[tranche]))
    vep = process_consequences(vep)
    var = var.annotate_rows(gene_id=vep[var.row_key].vep.worst_csq_by_gene_canonical.gene_id)
    var = var.explode_rows(var.gene_id)
    var = var.annotate_rows(annotation=hl.if_else(hl.literal({'missense', 'LC'}).contains(var.annotation), 'missense|LC', var.annotation), )
    var = var.group_rows_by('gene_id', 'gene', 'annotation').aggregate(var_cnt=hl.agg.count_where(var.Pvalue < level))
    test_type = test_type.lower()
    if test_type == 'skat':
        pvalue = 'Pvalue_SKAT'
    elif test_type == 'burden':
        pvalue = 'Pvalue_Burden'
    else:
        pvalue = 'Pvalue'
    mt = mt.annotate_entries(var_cnt=var[mt.row_key, mt.col_key]['var_cnt'])
    mt = mt.annotate_cols(pheno_gene_var_sig_cnt=hl.agg.count_where(hl.is_defined(mt.var_cnt) & (mt.var_cnt > 0) &
                                                                  hl.is_defined(mt[pvalue]) & (mt[pvalue] < level)),
                          pheno_var_sig_cnt=hl.agg.count_where((hl.is_defined(mt.var_cnt) & (mt.var_cnt > 0)) &
                                                             (~(hl.is_defined(mt[pvalue]) & (mt[pvalue] < level)))),
                          pheno_gene_sig_cnt=hl.agg.count_where((hl.is_defined(mt[pvalue]) & (mt[pvalue] < level)) &
                                                              (~(hl.is_defined(mt.var_cnt) & (mt.var_cnt > 0)))),
                          pheno_none_sig_cnt=hl.agg.count_where(~((hl.is_defined(mt[pvalue]) & (mt[pvalue] < level)) |
                                                                (hl.is_defined(mt.var_cnt) & (mt.var_cnt > 0)))))

    mt = mt.annotate_rows(gene_var_sig_cnt=hl.agg.count_where(hl.is_defined(mt.var_cnt) & (mt.var_cnt > 0) &
                                                                  hl.is_defined(mt[pvalue]) & (mt[pvalue] < level)),
                          var_sig_cnt=hl.agg.count_where((hl.is_defined(mt.var_cnt) & (mt.var_cnt > 0)) &
                                                             (~(hl.is_defined(mt[pvalue]) & (mt[pvalue] < level)))),
                          gene_sig_cnt=hl.agg.count_where((hl.is_defined(mt[pvalue]) & (mt[pvalue] < level)) &
                                                              (~(hl.is_defined(mt.var_cnt) & (mt.var_cnt > 0)))),
                          none_sig_cnt=hl.agg.count_where(((hl.is_missing(mt[pvalue]) | (mt[pvalue] >= level)) &
                                                           (hl.is_missing(mt.var_cnt) | (mt.var_cnt == 0)))))
    return mt

def compute_mean_coverage_ht(tranche: str = CURRENT_TRANCHE):
    int_auto = hl.read_table(interval_qc_path(data_source='broad', freeze=7, chrom='autosomes'))
    int_sex = hl.read_table(interval_qc_path(data_source='broad', freeze=7, chrom='sex_chr'))
    int_sex = int_sex.select(target_mean_dp=int_sex.target_mean_dp['XX'], 
                            target_pct_gt_10x=int_sex.target_pct_gt_10x['XX'], 
                            target_pct_gt_20x=int_sex.target_pct_gt_20x['XX'], 
                            pct_samples_10x=int_sex.pct_samples_10x['XX'], 
                            pct_samples_20x=int_sex.pct_samples_20x['XX'])
    int_full = int_auto.union(int_sex)

    var = hl.read_matrix_table(get_results_mt_path('variant', tranche=tranche))
    vep = hl.read_table(var_annotations_ht_path('vep', *TRANCHE_DATA[tranche]))
    vep = process_consequences(vep)
    var = var.annotate_rows(gene_id=vep[var.row_key].vep.worst_csq_by_gene_canonical.gene_id, 
                            annotation=hl.if_else(hl.is_definedhl.literal({'missense', 'LC'}).contains(var.annotation), 'missense|LC', var.annotation),
                            coverage=int_full[var.locus].target_mean_dp)
    var = var.explode_rows(var.gene_id)
    var = var.rows()
    mean_coverage = var.group_by('gene_id', 'gene', 'annotation').aggregate(mean_coverage=hl.agg.mean(var.coverage))
    return mean_coverage

def tie_breaker(l, r):
    return (hl.case()\
        .when(l.n_cases_both_sexes > r.n_cases_both_sexes, -1)\
        .when(l.n_cases_both_sexes == r.n_cases_both_sexes, 0)\
        .when(l.n_cases_both_sexes < r.n_cases_both_sexes, 1)\
        .or_missing())

def get_random_pheno_pvalue_ht(result_type: str = 'gene', test_type: str = 'skato', coding: list = None, tranche: str = CURRENT_TRANCHE, af_lower: float = None, af_upper: float = None, random_phenos: bool = False):
    af = get_caf_info_ht(result_type, tranche)
    mt = hl.read_matrix_table(get_results_mt_path(result_type, tranche, random_phenos=random_phenos))
    if coding is not None:
        mt = mt.filter_cols(hl.literal(coding).contains(mt.coding))
    if result_type == 'gene':  # Gene Info
        if af_lower is not None:
            af = af.filter(af.CAF > af_lower)
        if af_upper is not None:
            af = af.filter(af.CAF <= af_upper)
        mt = mt.filter_rows(hl.is_defined(af.index(mt['annotation'], mt['gene_id'])))
        mt = mt.annotate_rows(af=f'CAF:({af_lower}, {af_upper}]')
    else:  # Variant Info
        if af_lower is not None:
            af = af.filter(af.AF > af_lower)
        if af_upper is not None:
            af = af.filter(af.AF <= af_upper)
        mt = mt.filter_rows(hl.is_defined(af.index(mt['locus'], mt['alleles'])))
        mt = mt.annotate_rows(annotation=hl.if_else(hl.literal({'missense', 'LC'}).contains(mt.annotation), 'missense|LC', mt.annotation),
                              af=f'AF:({af_lower}, {af_upper}]')
    ht = mt.entries()
    test_type = test_type.lower()
    if test_type == 'skat':
        pvalue = 'Pvalue_SKAT'
    elif test_type == 'burden':
        pvalue = 'Pvalue_Burden'
    else:
        pvalue = 'Pvalue'
    ht = ht.select(pvalue, 'af')
    return ht

def get_lambda_filtered_ht(lambda_ht: hl.Table = None, lambda_name: str = 'lambda_gc_skato', lower:float = None, upper: float = None):
    if lower is not None:
        ht = lambda_ht.filter(lambda_ht[lambda_name] > lower)
    if upper is not None:
        ht = lambda_ht.filter(lambda_ht[lambda_name] < upper)
    return ht

def get_corr_phenos_ht(r_2: float = None, tie_breaker = None, tranche: str = CURRENT_TRANCHE):
    pheno_mt = get_ukb_pheno_mt()
    ht = hl.read_matrix_table(get_results_mt_path(tranche)).cols()
    pheno_mt = pheno_mt.filter_cols(hl.is_defined(ht[pheno_mt.col_key]))
    corr = make_pairwise_ht(pheno_mt, pheno_field=pheno_mt.both_sexes, correlation=True)
    related = corr.filter((corr.entry > hl.sqrt(r_2)) & (corr.i != corr.j))
    pheno_to_remove = hl.maximal_independent_set(related.i_data, related.j_data, keep=False, tie_breaker=tie_breaker)
    return pheno_to_remove
