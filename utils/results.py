from ukbb_qc.resources.sample_qc import interval_qc_path
from ukbb_qc.resources.variant_qc import var_annotations_ht_path
from ukbb_qc.resources.basics import release_ht_path
from ukb_common.utils.annotations import annotation_case_builder
from gnomad.utils.vep import process_consequences
from ukb_exomes.resources import *


def compute_lambda_gc_ht(result_type: str = 'gene', tranche: str = CURRENT_TRANCHE, cutoff: float = 0.0001):
    af = get_af_info_ht(result_type, tranche)
    mt = hl.read_matrix_table(get_results_mt_path(result_type, tranche))

    if result_type == 'gene':  # Gene Info
        af = af.filter(af.caf > cutoff)
        mt = mt.filter_rows(hl.is_defined(af.index(mt['annotation'], mt['gene_id'])))
        lambda_gc = hl.agg.filter(hl.is_defined(mt.Pvalue), hl.methods.statgen._lambda_gc_agg(mt.Pvalue))
        mt = mt.select_cols('n_cases', lambda_gc_skato=lambda_gc)
    else:  # Variant Info
        af = af.filter(af.AF > cutoff)
        mt = mt.filter_rows(hl.is_defined(af.index(mt['locus'], mt['alleles'])))
        lambda_gc = hl.agg.filter(hl.is_defined(mt.Pvalue), hl.methods.statgen._lambda_gc_agg(mt.Pvalue))
        mt = mt.select_cols('n_cases', lambda_gc=lambda_gc)

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


def get_af_info_ht(result_type: str = 'gene', tranche: str = CURRENT_TRANCHE):
    vep = hl.read_table(var_annotations_ht_path('vep', *TRANCHE_DATA[tranche]))
    vep = process_consequences(vep)
    vep = vep.explode(vep.vep.worst_csq_by_gene_canonical)
    annotation = annotation_case_builder(vep.vep.worst_csq_by_gene_canonical)
    vep = vep.annotate(annotation=hl.if_else(hl.literal({'missense', 'LC'}).contains(annotation), 'missense|LC', annotation))
    vep = vep.filter((hl.is_defined(vep.annotation)) & (vep.annotation != 'non-coding'))
    vep = vep.select(vep.vep.worst_csq_by_gene_canonical.gene_id, vep.annotation)

    freq = hl.read_table(release_ht_path(*TRANCHE_DATA[tranche]))
    var_af = vep.annotate(AF=freq[vep.key].freq[0].AF)

    if result_type == 'gene':
        sum_af = var_af.group_by(var_af.annotation, var_af.gene_id).aggregate(caf=hl.agg.sum(var_af.AF))
        return sum_af
    else:
        return var_af


def get_sig_cnt_mt(phenocode: list = None, gene_symbol: list = None, level: float = 1e-6, tranche: str = CURRENT_TRANCHE):
    # Load and Combine Data
    af = get_af_info_ht(result_type='gene', tranche=tranche)
    mt = hl.read_matrix_table(get_results_mt_path(tranche=tranche))
    mt = mt.annotate_rows(caf=af[mt.gene_id, mt.annotation]['caf'])
    # Count Significant Hits per Gene for each Test
    mt = mt.annotate_rows(sig_pheno_cnt_skato=hl.agg.sum(mt.Pvalue < level),
                          sig_pheno_cnt_skat=hl.agg.sum(mt.Pvalue_SKAT < level),
                          sig_pheno_cnt_burden=hl.agg.sum(mt.Pvalue_Burden < level))
    # Count Significant Hits per Phenotype for each Test
    mt = mt.annotate_cols(sig_gene_cnt_skato=hl.agg.sum(mt.Pvalue < level),
                          sig_gene_cnt_skat=hl.agg.sum(mt.Pvalue_SKAT < level),
                          sig_gene_cnt_burden=hl.agg.sum(mt.Pvalue_Burden < level))
    if phenocode is not None:
        mt = mt.filter_cols(hl.literal(phenocode).contains(mt.phenocode))
    if gene_symbol is not None:
        mt = mt.filter_rows(hl.literal(gene_symbol).contains(mt.gene_symbol))
    return mt

def compare_gene_var_sig_cnt_ht(test_type: str = 'skato', level: float = 1e-6, tranche: str = CURRENT_TRANCHE):
    var = hl.read_matrix_table(get_results_mt_path('variant', tranche=tranche))
    gene = hl.read_matrix_table(get_results_mt_path(tranche=tranche))

    vep = hl.read_table(var_annotations_ht_path('vep', *TRANCHE_DATA[tranche]))
    vep = process_consequences(vep)
    var = var.annotate_rows(gene_id=vep[var.row_key].vep.worst_csq_by_gene_canonical.gene_id)
    var = var.explode_rows(var.gene_id)
    var = var.annotate_rows(annotation=hl.if_else(hl.literal({'missense', 'LC'}).contains(var.annotation), 'missense|LC', var.annotation), )
    var = var.group_rows_by('gene_id', 'gene', 'annotation').aggregate(var_cnt=hl.agg.count_where(var.Pvalue<level))

    test_type = test_type.lower()
    if test_type == 'skat':
        pvalue = 'Pvalue_SKAT'
    elif test_type == 'burden':
        pvalue = 'Pvalue_Burden'
    else:
        pvalue = 'Pvalue'

    gene = gene.annotate_entries(var_cnt=var[gene.row_key, gene.col_key]['var_cnt'])
    gene = gene.annotate_cols(gene_var_sig_cnt=hl.agg.count_where(hl.is_defined(gene.var_cnt) & (gene.var_cnt > 0) &
                                                                  hl.is_defined(gene[pvalue]) & (gene[pvalue] < level)),
                              var_sig_cnt=hl.agg.count_where((hl.is_defined(gene.var_cnt) & (gene.var_cnt > 0)) &
                                                             (~(hl.is_defined(gene[pvalue]) & (gene[pvalue] < level)))),
                              gene_sig_cnt=hl.agg.count_where((hl.is_defined(gene[pvalue]) & (gene[pvalue] < level)) &
                                                              (~(hl.is_defined(gene.var_cnt) & (gene.var_cnt > 0)))),
                              none_sig_cnt=hl.agg.count_where(~((hl.is_defined(gene[pvalue]) & (gene[pvalue] < level)) |
                                                                (hl.is_defined(gene.var_cnt) & (gene.var_cnt > 0)))))
    return gene.cols()

def compute_mean_coverage_ht(tranche: str = CURRENT_TRANCHE):
    int_auto = hl.read_table(interval_qc_path(data_source='broad', freeze=7, chrom='autosomes')) # freeze = 7
    int_sex = hl.read_table(interval_qc_path(data_source='broad', freeze=7, chrom='sex_chr')) # freeze = 7
    int_xx = int_sex.select(**{'target_mean_dp': int_sex.target_mean_dp[chrom] for chrom in {'XX'}},
                            **{'target_pct_gt_10x': int_sex.target_pct_gt_10x[chrom] for chrom in {'XX'}},
                            **{'target_pct_gt_20x': int_sex.target_pct_gt_20x[chrom] for chrom in {'XX'}},
                            **{'pct_samples_10x': int_sex.pct_samples_10x[chrom] for chrom in {'XX'}},
                            **{'pct_samples_20x': int_sex.pct_samples_20x[chrom] for chrom in {'XX'}})
    int_xy = int_sex.select(**{'target_mean_dp': int_sex.target_mean_dp[chrom] for chrom in {'XY'}},
                            **{'target_pct_gt_10x': int_sex.target_pct_gt_10x[chrom] for chrom in {'XY'}},
                            **{'target_pct_gt_20x': int_sex.target_pct_gt_20x[chrom] for chrom in {'XY'}},
                            **{'pct_samples_10x': int_sex.pct_samples_10x[chrom] for chrom in {'XY'}},
                            **{'pct_samples_20x': int_sex.pct_samples_20x[chrom] for chrom in {'XY'}})
    int_sex = int_xx.union(int_xy)
    int_full = int_auto.union(int_sex)

    var = hl.read_matrix_table(get_results_mt_path('variant', tranche=tranche))
    vep = hl.read_table(var_annotations_ht_path('vep', *TRANCHE_DATA[tranche]))
    vep = process_consequences(vep)
    # Previous processing:
    # vep = vep.explode(vep.vep.worst_csq_by_gene_canonical)
    # vep = vep.select(vep.vep.worst_csq_by_gene_canonical.gene_id)
    # var = var.annotate_rows(gene_id=vep[var.row_key].gene_id)
    var = var.annotate_rows(gene_id=vep[var.row_key].vep.worst_csq_by_gene_canonical.gene_id)
    var = var.explode_rows(var.gene_id)
    var = var.annotate_rows(annotation=hl.if_else(hl.literal({'missense', 'LC'}).contains(var.annotation), 'missense|LC', var.annotation), )
    var = var.key_rows_by('locus')
    var = var.annotate_rows(coverage=int_full[var.row_key].target_mean_dp)
    var = var.rows()
    mean_coverage = var.group_by('gene_id', 'gene', 'annotation').aggregate(mean_coverage=hl.agg.mean(var.coverage))
    return mean_coverage
