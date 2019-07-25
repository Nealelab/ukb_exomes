#!/usr/bin/env python3

__author__ = 'konradk'

from gnomad_hail import *
# import ukbb_qc.resources as ukb
from ukb_exomes import *
import shlex
from pipeline import *
import pipeline.pipeline
import time

bucket = 'gs://ukbb-pharma-exome-analysis'
root = f'{bucket}/data'
ukb_for_grm_plink_path = f'{bucket}/mt/ukb.for_grm.pruned.plink'
pheno_data_index_path = 'gs://ukbb-pharma-exome-analysis/pheno_{data_type}/index.tsv'
saige_pheno_types = {
    'continuous': 'quantitative',
    'categorical': 'binary',
    'icd': 'binary'
}

HAIL_DOCKER_IMAGE = 'konradjk/hail_utils:latest'
SAIGE_DOCKER_IMAGE = 'wzhou88/saige:0.35.8.5'
# SAIGE_DOCKER_IMAGE = 'konradjk/saige:0.35.8.2.2'


def read_pheno_index_file(pheno_index_path: str):
    files = {}
    with hl.hadoop_open(pheno_index_path, 'r') as f:
        for line in f:
            fname, line_data = line.strip().split('\t')
            line_data = json.loads(line_data)
            if 'pheno' in line_data and 'coding' in line_data:
                line_data['id'] = f'{line_data["pheno"]}-{line_data["coding"]}'
            else:
                line_data['id'] = line_data['icd_code']
            files[fname] = line_data
    return files


def create_sparse_grm(p: Pipeline, output_path: str, plink_file_root: str,
                      relatedness_cutoff: str = '0.125', num_markers: int = 2000,
                      n_threads: int = 8, storage = '1500Mi'):
    in_bfile = p.read_input_group(
        **{ext: f'{plink_file_root}.{ext}' for ext in ('bed', 'bim', 'fam')})
    create_sparse_grm_task: pipeline.pipeline.Task = p.new_task(name='create_sparse_grm')
    create_sparse_grm_task.cpu(n_threads).storage(storage)
    create_sparse_grm_task.declare_resource_group(
        sparse_grm={ext: f'{{root}}{ext}' for ext in
                   (f'_relatednessCutoff_{relatedness_cutoff}_{num_markers}_randomMarkersUsed.sparseGRM.mtx',
                    f'_relatednessCutoff_{relatedness_cutoff}_{num_markers}_randomMarkersUsed.sparseGRM.mtx.sampleIDs.txt')})
    command = (f'Rscript /usr/local/bin/createSparseGRM.R '
               f'--plinkFile={in_bfile} '
               f'--nThreads={n_threads} '
               f'--outputPrefix={create_sparse_grm_task.sparse_grm}	'
               f'--numRandomMarkerforSparseKin={num_markers} '
               f'--relatednessCutoff={relatedness_cutoff}')
    create_sparse_grm_task.command(command)
    p.write_output(create_sparse_grm_task.sparse_grm, output_path)
    # Runtime: ~40 minutes on 8 cores
    return create_sparse_grm_task.sparse_grm


def extract_vcf_from_mt(p: Pipeline, output_root: str, gene_map_ht_path: str, mt_path: str, sample_mapping_file: str,
                        gene: str = None, interval: str = None, groups=None,
                        set_missing_to_hom_ref: bool = False, callrate_filter: bool = False, adj: bool = True,
                        export_bgen: bool = True,
                        n_threads: int = 2, storage: str = '500Mi', memory: str = ''):
    if groups is None:
        # groups = {'pLoF', 'missense|LC', 'pLoF|missense|LC', 'synonymous'}
        groups = {'pLoF', 'missense|LC', 'synonymous'}
    extract_task: pipeline.pipeline.Task = p.new_task(name='extract_vcf',
                                                      attributes={
                                                          'interval': interval
                                                      })
    extract_task.image(HAIL_DOCKER_IMAGE).cpu(n_threads).storage(storage)
    if export_bgen:
        extract_task.declare_resource_group(out={'bgen': '{root}.bgen',
                                                 'sample': '{root}.sample',
                                                 'bgen.bgi': '{root}.bgen.bgi'})
    else:
        extract_task.declare_resource_group(out={'vcf.gz': f'{{root}}.vcf.gz',
                                                 'vcf.gz.tbi': f'{{root}}.vcf.gz.tbi'})

    python_command = f"""import hail as hl
hl.init(master='local[{n_threads}]', default_reference='GRCh38')
gene_ht = hl.read_table('{gene_map_ht_path}')
mt = hl.read_matrix_table('{mt_path}')
meta_ht = hl.read_table('gs://broad-ukbb/regeneron.freeze_4/sample_qc/meta_w_pop_adj.ht')
meta_ht = meta_ht.filter(~meta_ht.is_filtered & (meta_ht.hybrid_pop == '12'))
sample_mapping_ht = hl.import_table('{sample_mapping_file}', key='eid_sample', delimiter=',')
qual_ht = hl.read_table('gs://broad-ukbb/broad.freeze_4/variant_qc/score_rankings/vqsr.ht')
qual_ht = qual_ht.filter(hl.len(qual_ht.filters) == 0)"""

    if gene is not None:
        python_command += f"""
gene_ht = gene_ht.filter(gene_ht.gene_symbol == '{gene}')
interval = gene_ht.aggregate(hl.agg.take(gene_ht.interval, 1), _localize=False)[0]"""
    elif interval is not None:
        python_command += f"""
interval = [hl.parse_locus_interval('{interval}')]
gene_ht = hl.filter_intervals(gene_ht, interval)"""

    python_command += (f"""
gene_ht = gene_ht.filter(hl.set({{'{"', '".join(groups)}'}}).contains(gene_ht.annotation))
gene_ht.select(group=gene_ht.gene_id + '_' + gene_ht.gene_symbol + '_' + gene_ht.annotation, variant=hl.delimit(gene_ht.variants, '"\\\\t"')).key_by().drop('start').export('{extract_task.group_file}', header=False)""")

    if adj:
        python_command += "\nmt = mt.filter_entries(mt.adj)"

    python_command += (f"""
mt = hl.filter_intervals(mt, interval).select_entries('GT')
mt = mt.annotate_cols(exome_id=sample_mapping_ht[mt.s.split('_')[1]].eid_26041, meta=meta_ht[mt.col_key]).key_cols_by('exome_id')
mt = mt.filter_cols(hl.is_defined(mt.meta) & hl.is_defined(mt.exome_id))
mt = mt.filter_rows(hl.is_defined(qual_ht[mt.row_key]) & (hl.agg.count_where(mt.GT.is_non_ref()) > 0))
mt = mt.annotate_rows(rsid=mt.locus.contig + ':' + hl.str(mt.locus.position) + '_' + mt.alleles[0] + '/' + mt.alleles[1])""")

    if callrate_filter:
        python_command += "\nmt = mt.filter_rows(hl.agg.fraction(hl.is_defined(mt.GT)) > 0.95)"

    if export_bgen:
        dosage_command = """mt = mt.annotate_entries(gp=hl.or_missing(
                                                         hl.is_defined(mt.GT),
                                                         hl.map(lambda i: hl.cond(mt.GT.unphased_diploid_gt_index() == i, 1.0, 0.0),
                                                                hl.range(0, hl.triangle(hl.len(mt.alleles))))))"""
        python_command += '\n' + dosage_command.replace("\n", " ")
        if set_missing_to_hom_ref:
            impute_string = "[1.0, 0.0, 0.0]"
        else:
            impute_string = "mt.mean_gp"
            python_command += "\nmt = mt.annotate_rows(mean_gp=hl.agg.array_agg(lambda x: hl.agg.mean(x), mt.gp))"
        python_command += f"\nmt = mt.annotate_entries(gp=hl.or_else(mt.gp, {impute_string}))"
        python_command += f"\nhl.export_bgen(mt, '{extract_task.out}', gp=mt.gp, varid=mt.rsid)"
    else:
        if set_missing_to_hom_ref:
            python_command += "\nmt = mt.annotate_entries(GT=hl.or_else(mt.GT, hl.call(0, 0)))"
        # Note: no mean-imputation for VCF
        python_command += f"\nhl.export_vcf(mt, '{extract_task.bgz}.bgz')"

    python_command = python_command.replace('\n', '; ').strip()
    command = f'python3 -c "{python_command}";'

    if export_bgen:
        command += f'\n/bgen_v1.1.4-Ubuntu16.04-x86_64/bgenix -g {extract_task.out.bgen} -index -clobber'
    else:
        command += f'\nmv {extract_task.bgz}.bgz {extract_task.out["vcf.gz"]}; tabix {extract_task.out["vcf.gz"]};'
    extract_task.command(command.replace('\n', ' '))

    p.write_output(extract_task.out, output_root)
    p.write_output(extract_task.group_file, f'{output_root}.gene.txt')
    return extract_task


def gt_to_gp(mt, location: str = 'gp'):
    return mt.annotate_entries(**{location: hl.or_missing(
        hl.is_defined(mt.GT),
        hl.map(lambda i: hl.cond(mt.GT.unphased_diploid_gt_index() == i, 1.0, 0.0),
               hl.range(0, hl.triangle(hl.len(mt.alleles)))))})


def impute_missing_gp(mt, location: str = 'gp', mean_impute: bool = True):
    mt = mt.annotate_entries(_gp = mt.location)
    if mean_impute:
        mt = mt.annotate_rows(_mean_gp=hl.agg.array_agg(lambda x: hl.agg.mean(x), mt._gp))
        gp_expr = mt._mean_gp
    else:
        gp_expr = [1.0, 0.0, 0.0]
    return mt.annotate_entries(**{location: hl.or_else(mt._gp, gp_expr)}).drop('_gp')


def fit_null_glmm(p: Pipeline, output_root: str, pheno_path: str, pheno_name: str, trait_type: str, covariates: str,
                  plink_file_root: str, sparse_grm: pipeline.pipeline.Resource = None,
                  sparse_grm_extension: str = None, skip_model_fitting: bool = False,
                  n_threads: int = 8, storage: str = '1500Mi'):
    analysis_type = "variant" if sparse_grm is None else "gene"
    pheno_col = 'any_codes' if trait_type == 'icd' else 'both_sexes'
    user_id_col = 'userId' if trait_type == 'icd' else 'userID'  # TODO: fix on next load
    pheno_file = p.read_input_group(**{'gz': pheno_path})
    in_bfile = p.read_input_group(**{ext: f'{plink_file_root}.{ext}' for ext in ('bed', 'bim', 'fam')})
    fit_null_task = p.new_task(name=f'fit_null_model',
                               attributes={
                                   'analysis_type': analysis_type,
                                   'pheno': pheno_name
                               }).cpu(n_threads).storage(storage)
    output_files = {ext: f'{{root}}{ext if ext.startswith("_") else "." + ext}' for ext in
                   ('rda', '_30markers.SAIGE.results.txt', f'{analysis_type}.varianceRatio.txt')}
    if analysis_type == 'gene':
        sparse_sigma_extension = sparse_grm_extension.replace("GRM", "Sigma")
        output_files[f'{analysis_type}.varianceRatio.txt{sparse_sigma_extension}'] = \
            f'{{root}}.{analysis_type}.varianceRatio.txt{sparse_sigma_extension}'
    fit_null_task.declare_resource_group(null_glmm=output_files)
    bim_fix_command = f'perl -pi -e s/^chr// {in_bfile.bim}'
    if trait_type == 'icd':
        bim_fix_command += (f"; zcat {pheno_file.gz} | perl -p -e 's/true/1/g' | perl -p -e 's/false/0/g' "
                            f"| gzip -c > {pheno_file.gz}.temp.gz; mv {pheno_file.gz}.temp.gz {pheno_file.gz}")

    command = (f'Rscript /usr/local/bin/step1_fitNULLGLMM.R '
               f'--plinkFile={in_bfile} '
               f'--phenoFile={pheno_file.gz} '
               f'--covarColList={covariates} '
               f'--phenoCol={pheno_col} '
               f'--sampleIDColinphenoFile={user_id_col} '
               f'--traitType={saige_pheno_types[trait_type]} '
               f'--outputPrefix={fit_null_task.null_glmm} '
               f'--outputPrefix_varRatio={fit_null_task.null_glmm}.{analysis_type} '
               f'--skipModelFitting={str(skip_model_fitting).upper()} ')
    if analysis_type == "gene":
        fit_null_task.declare_resource_group(sparse_sigma={sparse_sigma_extension: f'{{root}}.{sparse_sigma_extension}'})
        command += (f'--IsSparseKin=TRUE '
                    f'--sparseGRMFile={sparse_grm[sparse_grm_extension]} '
                    f'--sparseGRMSampleIDFile={sparse_grm[f"{sparse_grm_extension}.sampleIDs.txt"]} '
                    f'--isCateVarianceRatio=TRUE ')
    command += f'--nThreads={n_threads} --LOCO=FALSE 2>&1 | tee {fit_null_task.stdout}'
    command = '; '.join([bim_fix_command, command])
    fit_null_task.command(command)
    p.write_output(fit_null_task.null_glmm, output_root)
    p.write_output(fit_null_task.stdout, f'{output_root}.log')
    # Runtimes: 8 threads: ~5 minutes of 100% CPU (~3G RAM), followed by ~9 minutes of 800% (~6G RAM)
    return fit_null_task


def run_saige(p: Pipeline, output_root: str, model_file: str, variance_ratio_file: str,
              vcf_file: pipeline.pipeline.ResourceGroup, samples_file: pipeline.pipeline.ResourceGroup,
              group_file: str = None, sparse_sigma_file: str = None, use_bgen: bool = True,
              trait_type: str = 'continuous',
              chrom: str = 'chr1', min_mac: int = 1, min_maf: float = 0, max_maf: float = 0.5,
              storage: str = '500Mi'):

    analysis_type = "gene" if sparse_sigma_file is not None else "variant"
    run_saige_task: pipeline.pipeline.Task = p.new_task(name=f'run_saige',
                                                        attributes={
                                                            'analysis_type': analysis_type,
                                                            'output_path': output_root,
                                                            'chromosome': chrom
                                                        }).cpu(1).storage(storage)  # Step 2 is single-threaded only

    run_saige_task.declare_resource_group(result={'gene.txt': '{root}',
                                                  'single.txt': '{root}_single'})

    command = (f'head -2 {group_file} | tac > {group_file}.tmp; mv {group_file}.tmp {group_file}; '
               f'Rscript /usr/local/bin/step2_SPAtests.R '
               f'--minMAF={min_maf} '
               f'--minMAC={min_mac} '
               f'--maxMAFforGroupTest={max_maf} '
               f'--sampleFile={samples_file} '
               f'--GMMATmodelFile={model_file} '
               f'--varianceRatioFile={variance_ratio_file} '
               f'--SAIGEOutputFile={run_saige_task.result} ')
    if saige_pheno_types[trait_type] == 'binary':
        command += f'--IsOutputPvalueNAinGroupTestforBinary=TRUE '

    if use_bgen:
        command += (f'--bgenFile={vcf_file.bgen} '
                    f'--bgenFileIndex={vcf_file["bgen.bgi"]} ')
    else:
        command += (f'--vcfFile={vcf_file["vcf.gz"]} '
                    f'--vcfFileIndex={vcf_file["vcf.gz.tbi"]} '
                    f'--chrom={chrom} '
                    f'--vcfField=GT ')
    if analysis_type == "gene":
        command += (f'--groupFile={group_file} '
                    f'--sparseSigmaFile={sparse_sigma_file} '
                    f'--IsSingleVarinGroupTest=TRUE ')
    command += f'--IsOutputAFinCaseCtrl=TRUE 2>&1 | tee {run_saige_task.stdout}; '
    command += f"input_length=$(wc -l {group_file} | awk '{{print $1}}'); " \
        f"output_length=$(wc -l {run_saige_task.result['gene.txt']} | awk '{{print $1}}'); " \
        f"echo 'Got input:' $input_length 'output:' $output_length | tee -a {run_saige_task.stdout}; " \
        f"if [[ $input_length > 0 ]]; then echo 'got input' | tee -a {run_saige_task.stdout}; " \
        f"if [[ $output_length == 1 ]]; then echo 'but not enough output' | tee -a {run_saige_task.stdout}; exit 1; fi; fi"
    run_saige_task.command(command)
    p.write_output(run_saige_task.result, output_root)
    p.write_output(run_saige_task.stdout, f'{output_root}.log')
    # Runtimes: PCSK9 (4 groups, 31v, 289v, 320v, 116v): 22 minutes


def get_tasks_from_pipeline(p):
    return dict(Counter(map(lambda x: x.name, p.select_tasks(""))))


# TODO:
#  continuous: height: 12144/50, bmi: 23104
#  dichotomous: asthma: ICD J45.9, CF: E84.*, crohn's: K50.*, UC: K51.*, t1d: E10.*, t2d: E11.*
def main(args):
    trait_type = args.trait_type

    if trait_type == 'icd':
        phenos_to_run = {'K519'} #, 'K509', 'E10.9', 'E11.9', 'J45.9'}
    else:
        phenos_to_run = {'50-raw'}
        # phenos_to_run = {'50-raw', '699-raw', '23104-raw'}

    num_pcs = 20
    hl.init(log='/tmp/saige_temp_hail.log')
    start_time = time.time()
    pheno_data = read_pheno_index_file(pheno_data_index_path.format(data_type=trait_type))
    covariates = ','.join(['sex', 'age', 'age2', 'age_sex'] + [f'pc{x}' for x in range(1, num_pcs + 1)])
    relatedness_cutoff = '0.125'
    num_markers = 2000
    n_threads = 6
    chromosomes = range(7, 8)  # 23
    # chromosomes = range(1, 23)  # 23

    sparse_grm_root = f'{root}/tmp/sparse'
    sparse_grm_extension = f'_relatednessCutoff_{relatedness_cutoff}_{num_markers}_randomMarkersUsed.sparseGRM.mtx'

    if args.local_test:
        backend = LocalBackend(gsa_key_file='/Users/konradk/.hail/ukb_exomes.json')
    else:
        backend = BatchBackend(url='https://batch.hail.is')
    p = Pipeline(name='saige', backend=backend, default_image=SAIGE_DOCKER_IMAGE,
                 # default_memory='1Gi',
                 default_storage='500Mi', default_cpu=n_threads)

    # if args.variant_test:
    #     model_created = True  # hl.hadoop_exists(model_file_path) and not args.overwrite
    #     variant_var_ratio_created = True  # hl.hadoop_exists(variant_variance_ratio_file_path) and not args.overwrite
    #     vcf_created = False  # hl.hadoop_exists(f'{vcf_root}.vcf.gz') and not args.overwrite
    #     pheno = phenos[0]
    #     variant_variance_ratio_file_path = f'{null_glmm_root}.variant.varianceRatio.txt'
    #     variant_var_ratio_created = True  # hl.hadoop_exists(variant_variance_ratio_file_path) and not args.overwrite
    #     if not (model_created and variant_var_ratio_created):
    #         fit_null_task = fit_null_glmm(p, null_glmm_root, pheno_path, pheno, trait_type, covariates,
    #                                       ukb_for_grm_plink_path,
    #                                       skip_model_fitting=model_created, n_threads=n_threads)
    #
    #     model_file = p.read_input(model_file_path) if model_created else fit_null_task.null_glmm.rda
    #     variance_ratio_file = p.read_input(variant_variance_ratio_file_path) if variant_var_ratio_created else \
    #         fit_null_task.null_glmm['variant.varianceRatio.txt']
    #
    #     vcf_task = None
    #     if not vcf_created:
    #         vcf_task = extract_vcf_from_mt(p, vcf_root, get_ukb_gene_map_ht_path(), get_ukb_exomes_mt_path(), sample_mapping_file, gene)
    #
    #     vcf_file = p.read_input_group(vcf={'vcf.gz': f'{vcf_root}.vcf.gz',
    #                                        'vcf.gz.tbi': f'{vcf_root}.vcf.gz.tbi'})
    #     run_saige(p, results_path, model_file, variance_ratio_file, vcf_file.vcf, get_ukb_samples_file_path(), dependency=vcf_task)

    if args.create_sparse_grm:
        sparse_grm = create_sparse_grm(p, sparse_grm_root, ukb_for_grm_plink_path,
                                       relatedness_cutoff, num_markers, n_threads=n_threads)
    else:
        sparse_grm = p.read_input_group(**{ext: f'{sparse_grm_root}.{ext}' for ext in
                                           (sparse_grm_extension, f'{sparse_grm_extension}.sampleIDs.txt')})

    overwrite_null_models = args.create_null_models
    null_model_dir = f'{root}/null_glmm'
    if not overwrite_null_models:
        null_models_already_created = {x['path'] for x in hl.hadoop_ls(null_model_dir)}
    null_models = {}

    for pheno_path, pheno_meta in pheno_data.items():
        if pheno_meta['id'] not in phenos_to_run:
            continue
        pheno = pheno_meta['id']
        null_glmm_root = f'{null_model_dir}/{pheno}'
        model_file_path = f'{null_glmm_root}.rda'
        gene_variance_ratio_file_path = f'{null_glmm_root}.gene.varianceRatio.txt'
        sparse_sigma_file_path = gene_variance_ratio_file_path + sparse_grm_extension.replace("GRM", "Sigma")

        if not overwrite_null_models and model_file_path in null_models_already_created:
            model_file = p.read_input(model_file_path)
            variance_ratio_file = p.read_input(gene_variance_ratio_file_path)
            sparse_sigma_file = p.read_input(sparse_sigma_file_path)
        else:
            fit_null_task = fit_null_glmm(p, null_glmm_root, pheno_path, pheno, trait_type, covariates,
                                          ukb_for_grm_plink_path,
                                          sparse_grm=sparse_grm, sparse_grm_extension=sparse_grm_extension,
                                          n_threads=n_threads)
            model_file = fit_null_task.null_glmm.rda
            variance_ratio_file = fit_null_task.null_glmm['gene.varianceRatio.txt']
            sparse_sigma_file = fit_null_task.null_glmm[
                f'gene.varianceRatio.txt{sparse_grm_extension.replace("GRM", "Sigma")}']
        null_models[pheno] = (model_file, variance_ratio_file, sparse_sigma_file)

    chrom_lengths = hl.get_reference('GRCh38').lengths

    use_bgen = args.use_bgen
    vcf_dir = f'{root}/vcf'
    test_extension = 'bgen' if use_bgen else 'vcf.gz'
    overwrite_vcfs = args.create_vcfs
    if not overwrite_vcfs and hl.hadoop_exists(vcf_dir):
        vcfs_already_created = {x['path'] for x in hl.hadoop_ls(vcf_dir)}
    chunk_size = int(1e6)
    vcfs = {}
    for chrom in chromosomes:
        chromosome = f'chr{chrom}'
        chrom_length = chrom_lengths[chromosome]
        for start_pos in range(1, chrom_length, chunk_size):
            if start_pos != 90000001: continue
            end_pos = chrom_length if start_pos + chunk_size > chrom_length else (start_pos + chunk_size)
            interval = f'{chromosome}:{start_pos}-{end_pos}'
            vcf_root = f'{vcf_dir}/test_{chromosome}_{str(start_pos).zfill(9)}'
            if not overwrite_vcfs and f'{vcf_root}.{test_extension}' in vcfs_already_created:
                if use_bgen:
                    vcf_file = p.read_input_group(**{'bgen': f'{vcf_root}.bgen',
                                                     'bgen.bgi': f'{vcf_root}.bgen.bgi',
                                                     'sample': f'{vcf_root}.sample'})
                else:
                    vcf_file = p.read_input_group(**{'vcf.gz': f'{vcf_root}.vcf.gz',
                                                     'vcf.gz.tbi': f'{vcf_root}.vcf.gz.tbi'})
                group_file = p.read_input(f'{vcf_root}.gene.txt')
            else:
                vcf_task = extract_vcf_from_mt(p, vcf_root, get_ukb_gene_map_ht_path(), get_ukb_exomes_mt_path(),
                                               sample_mapping_file, interval=interval, export_bgen=use_bgen)
                vcf_file = vcf_task.out
                group_file = vcf_task.group_file
            vcfs[interval] = (vcf_file, group_file)
            if args.local_test:
                break
        if args.local_test:
            break

    result_dir = f'{root}/result'
    overwrite_results = args.overwrite_results
    for pheno_path, pheno_meta in pheno_data.items():
        pheno = pheno_meta['id']
        if pheno not in null_models:
            continue

        pheno_results_dir = f'{result_dir}/{pheno}'
        results_already_created = {}
        if not overwrite_results and hl.hadoop_exists(pheno_results_dir):
            results_already_created = {x['path'] for x in hl.hadoop_ls(pheno_results_dir)}

        model_file, variance_ratio_file, sparse_sigma_file = null_models[pheno]
        for chrom in chromosomes:
            chromosome = f'chr{chrom}'
            chrom_length = chrom_lengths[chromosome]
            for start_pos in range(1, chrom_length, chunk_size):
                if start_pos != 90000001: continue
                end_pos = chrom_length if start_pos + chunk_size > chrom_length else (start_pos + chunk_size)
                interval = f'{chromosome}:{start_pos}-{end_pos}'
                vcf_file, group_file = vcfs[interval]
                results_path = f'{pheno_results_dir}/result_{pheno}_{chromosome}_{str(start_pos).zfill(9)}'
                if overwrite_results or f'{results_path}.txt' not in results_already_created:
                    samples_file = p.read_input(get_ukb_samples_file_path())
                    run_saige(p, results_path, model_file, variance_ratio_file, vcf_file, samples_file,
                              group_file, sparse_sigma_file, trait_type=trait_type, use_bgen=use_bgen)
                if args.local_test:
                    break
            if args.local_test:
                break
        if args.local_test:
            break

    print(f'Setup took: {time.strftime("%H:%M:%S", time.gmtime(time.time() - start_time))}')
    pprint(get_tasks_from_pipeline(p))
    p.run(dry_run=args.dry_run, verbose=True, delete_scratch_on_exit=False)
    pprint(get_tasks_from_pipeline(p))


if __name__ == '__main__':
    parser = argparse.ArgumentParser()

    parser.add_argument('--create_sparse_grm', help='Run single variant SAIGE', action='store_true')
    parser.add_argument('--create_null_models', help='Run single variant SAIGE', action='store_true')
    parser.add_argument('--create_vcfs', help='Run single variant SAIGE', action='store_true')
    parser.add_argument('--overwrite_results', help='Run single variant SAIGE', action='store_true')
    parser.add_argument('--local_test', help='Run single variant SAIGE', action='store_true')
    parser.add_argument('--use_bgen', help='Run single variant SAIGE', action='store_true')
    parser.add_argument('--trait_type', help='Run single variant SAIGE')
    parser.add_argument('--dry_run', help='Dry run only', action='store_true')
    parser.add_argument('--send_slack', help='Dry run only', action='store_true')
    args = parser.parse_args()

    if args.local_test:
        try_slack('@konradjk', main, args)
    else:
        main(args)


