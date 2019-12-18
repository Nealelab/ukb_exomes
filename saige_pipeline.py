#!/usr/bin/env python3

__author__ = 'konradk'

import time

# from pipeline import *
from gnomad_hail import *
from hailtop import pipeline
from hailtop.pipeline.pipeline import *
from ukb_exomes import *

bucket = 'gs://ukbb-pharma-exome-analysis'
# old_root = f'{bucket}/data'
root = f'{bucket}/data'
ukb_for_grm_plink_path = f'{bucket}/mt/ukb.for_grm.pruned.plink'
pheno_data_index_path = 'gs://ukbb-pharma-exome-analysis/pheno_{data_type}/index.tsv'
saige_pheno_types = {
    'continuous': 'quantitative',
    'biomarkers': 'quantitative',
    'categorical': 'binary',
    'icd': 'binary'
}

MIN_CASES = 200

HAIL_DOCKER_IMAGE = 'gcr.io/ukbb-exome-pharma/hail_utils:1.0'
SAIGE_DOCKER_IMAGE = 'wzhou88/saige:0.36.2'
SCRIPT_DIR = '/ukb_exomes/hail'


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

    output_file = f'{extract_task.bgz}.bgz' if export_bgen else extract_task.out
    command = f"""python3 {SCRIPT_DIR}/extract_vcf_from_mt.py
    --sample_mapping_file {sample_mapping_file}
    --mt_path {mt_path}
    --gene_map_ht_path {gene_map_ht_path}
    --meta_ht_path gs://broad-ukbb/regeneron.freeze_4/sample_qc/meta_w_pop_adj.ht
    --qual_ht_path gs://broad-ukbb/broad.freeze_4/variant_qc/score_rankings/vqsr.ht
    {"--gene " + gene if gene else ""}
    {"--interval " + interval if interval else ""}
    --groups {','.join(groups)}
    {"--callrate_filter" if callrate_filter else ""} 
    {"--export_bgen" if export_bgen else ""} 
    {"" if set_missing_to_hom_ref else "--mean_impute_missing"}
    {"" if adj else "--no_adj"} 
    --group_output_file {extract_task.group_file}
    --output_file {output_file} | tee {extract_task.stdout}
    ;""".replace('\n', ' ')

    if export_bgen:
        command += f'\n/bgen_v1.1.4-Ubuntu16.04-x86_64/bgenix -g {extract_task.out.bgen} -index -clobber'
    else:
        command += f'\nmv {extract_task.bgz}.bgz {extract_task.out["vcf.gz"]}; tabix {extract_task.out["vcf.gz"]};'
    extract_task.command(command.replace('\n', ' '))

    p.write_output(extract_task.out, output_root)
    p.write_output(extract_task.group_file, f'{output_root}.gene.txt')
    p.write_output(extract_task.stdout, f'{output_root}.log')
    return extract_task


def export_pheno(p: Pipeline, output_path: str, pheno: str, input_mt_path: str, covariates_path: str,
                 data_type: str = 'icd', n_threads: int = 8, storage: str = '500Mi'):
    extract_task: pipeline.pipeline.Task = p.new_task(name='extract_pheno',
                                                      attributes={
                                                          'pheno': pheno
                                                      })
    extract_task.image(HAIL_DOCKER_IMAGE).cpu(n_threads).storage(storage)
    coding = None
    if '-' in pheno:
        pheno, coding = pheno.split('-')
    python_command = f"""python3 {SCRIPT_DIR}/export_pheno.py
    --input_file {input_mt_path} 
    --covariates_path {covariates_path}
    --data_type {data_type}
    --pheno {pheno}
    {"--coding " + coding if coding else ''}
    --output_file {extract_task.out}
    --n_threads {n_threads} | tee {extract_task.stdout}
    ; """.replace('\n', ' ')

    extract_task.command(python_command)

    p.write_output(extract_task.out, output_path)
    p.write_output(extract_task.stdout, f'{output_path}.log')
    return extract_task.out


def fit_null_glmm(p: Pipeline, output_root: str, pheno_file: pipeline.pipeline.Resource, pheno_name: str,
                  trait_type: str, covariates: str,
                  plink_file_root: str, sparse_grm: pipeline.pipeline.Resource = None,
                  sparse_grm_extension: str = None, skip_model_fitting: bool = False,
                  n_threads: int = 8, storage: str = '1500Mi'):
    analysis_type = "variant" if sparse_grm is None else "gene"
    pheno_col = 'any_codes' if trait_type == 'icd' else 'both_sexes'
    user_id_col = 'userId' if trait_type == 'icd' else 'userID'  # TODO: fix on next load
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
    # if trait_type == 'icd':
    #     bim_fix_command += (f"; zcat {pheno_file.gz} | perl -p -e 's/true/1/g' | perl -p -e 's/false/0/g' "
    #                         f"| gzip -c > {pheno_file.gz}.temp.gz; mv {pheno_file.gz}.temp.gz {pheno_file.gz}")

    command = (f'Rscript /usr/local/bin/step1_fitNULLGLMM.R '
               f'--plinkFile={in_bfile} '
               f'--phenoFile={pheno_file} '
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
    p.write_output(fit_null_task.stdout, f'{output_root}.{analysis_type}.log')
    # Runtimes: 8 threads: ~5 minutes of 100% CPU (~3G RAM), followed by ~9 minutes of 800% (~6G RAM)
    return fit_null_task


def run_saige(p: Pipeline, output_root: str, model_file: str, variance_ratio_file: str,
              vcf_file: pipeline.pipeline.ResourceGroup, samples_file: pipeline.pipeline.ResourceGroup,
              group_file: str = None, sparse_sigma_file: str = None, use_bgen: bool = True,
              trait_type: str = 'continuous',
              chrom: str = 'chr1', min_mac: int = 1, min_maf: float = 0, max_maf: float = 0.5,
              memory: str = '', storage: str = '500Mi'):

    analysis_type = "gene" if sparse_sigma_file is not None else "variant"
    run_saige_task: pipeline.pipeline.Task = p.new_task(name=f'run_saige',
                                                        attributes={
                                                            'analysis_type': analysis_type,
                                                            'output_path': output_root,
                                                            'chromosome': chrom
                                                        }).cpu(1).storage(storage)  # Step 2 is single-threaded only

    if analysis_type == 'gene':
        run_saige_task.declare_resource_group(result={'gene.txt': '{root}',
                                                      'single.txt': '{root}_single'})
    else:
        run_saige_task.declare_resource_group(result={'single_variant.txt': '{root}'})

    command = (f'Rscript /usr/local/bin/step2_SPAtests.R '
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
    if analysis_type == 'gene':
        command += f"input_length=$(wc -l {group_file} | awk '{{print $1}}'); " \
            f"output_length=$(wc -l {run_saige_task.result['gene.txt']} | awk '{{print $1}}'); " \
            f"echo 'Got input:' $input_length 'output:' $output_length | tee -a {run_saige_task.stdout}; " \
            f"if [[ $input_length > 0 ]]; then echo 'got input' | tee -a {run_saige_task.stdout}; " \
            f"if [[ $output_length == 1 ]]; then echo 'but not enough output' | tee -a {run_saige_task.stdout}; exit 1; fi; fi"
    run_saige_task.command(command)
    p.write_output(run_saige_task.result, output_root)
    p.write_output(run_saige_task.stdout, f'{output_root}.{analysis_type}.log')
    return run_saige_task


def load_results_into_hail(p: Pipeline, output_root: str, pheno: str, tasks_to_hold,
                           n_threads: int = 8, storage: str = '500Mi'):

    load_data_task: pipeline.pipeline.Task = p.new_task(name=f'load_data',
                                                        attributes={
                                                            'output_path': output_root
                                                        }).image(HAIL_DOCKER_IMAGE).cpu(n_threads).storage(storage)
    load_data_task.depends_on(*tasks_to_hold)
    python_command = f"""python3 {SCRIPT_DIR}/load_results.py
    --input_dir {output_root}
    --pheno {pheno}
    --gene_map_ht_raw_path gs://ukbb-pharma-exome-analysis/mt/ukb.exomes.gene_map.raw.ht
    --ukb_vep_ht_path gs://ukbb-pharma-exome-analysis/mt/ukb.exomes.vep.ht
    --overwrite
    --n_threads {n_threads} | tee {load_data_task.stdout}
    ;""".replace('\n', ' ')

    python_command = python_command.replace('\n', '; ').strip()
    command = 'set -o pipefail; PYTHONPATH=$PYTHONPATH:/ PYSPARK_SUBMIT_ARGS="--conf spark.driver.memory=4g ' \
              '--conf spark.executor.memory=24g pyspark-shell" ' + python_command
    load_data_task.command(command)
    p.write_output(load_data_task.stdout, f'{output_root}/{pheno}_loading.log')


def get_tasks_from_pipeline(p):
    return dict(Counter(map(lambda x: x.name, p.select_tasks(""))))


def get_phenos_to_run(phenos_to_run, data_type: str, run_one: bool = False):
    with hl.hadoop_open(f'gs://ukbb-pharma-exome-analysis/misc/phenos_{data_type}.tsv', 'r') as f:
        header = f.readline().strip().split('\t')
        for line in f:
            fields = dict(zip(header, line.strip().split('\t')))
            if data_type == 'icd':
                if int(fields['n_cases_all']) >= MIN_CASES:
                    phenos_to_run.add(fields['icd_code'])
                    if run_one:
                        break
            else:
                if int(fields['score']) >= 1:
                    if data_type == 'categorical' and int(fields['n_cases']) < MIN_CASES:
                        continue
                    if fields['coding']:
                        pheno = f'{fields["pheno"]}-{fields["coding"]}'
                    else:
                        pheno = fields['pheno']
                    phenos_to_run.add(pheno)
                    if run_one:
                        break

    return phenos_to_run


def main(args):
    hl.init(log='/tmp/saige_temp_hail.log')
    trait_type = args.trait_type

    if trait_type == 'icd':
        # phenos_to_run = {'K519', 'K509', 'E109', 'E119', 'J459', 'I251'}
        phenos_to_run = {'K519', 'K509', 'E109', 'E119', 'J459', 'I251',
                         'K51', 'K50', 'E10', 'E11', 'J45', 'I25'}
        # phenos_to_run = set()
    else:
        phenos_to_run = {'50-raw', '699-raw', '23104-raw'}
    if args.run_all_phenos:
        phenos_to_run = set()
        phenos_to_run = get_phenos_to_run(phenos_to_run, trait_type, args.local_test)
    logger.info(f'Got {len(phenos_to_run)} phenotypes...')

    num_pcs = 20
    start_time = time.time()
    covariates = ','.join(['sex', 'age', 'age2', 'age_sex', 'age2_sex'] + [f'pc{x}' for x in range(1, num_pcs + 1)])
    relatedness_cutoff = '0.125'
    num_markers = 2000
    n_threads = 8
    chromosomes = list(range(1, 23)) # + ['X', 'Y']  # 23

    sparse_grm_root = f'{root}/tmp/sparse'
    # sparse_grm_root = f'{old_root}/tmp/sparse'
    sparse_grm_extension = f'_relatednessCutoff_{relatedness_cutoff}_{num_markers}_randomMarkersUsed.sparseGRM.mtx'

    if args.local_test:
        backend = pipeline.LocalBackend(gsa_key_file='/Users/konradk/.hail/ukb_exomes.json')
    else:
        backend = pipeline.BatchBackend(_service='batch2')
    p = pipeline.Pipeline(name='saige', backend=backend, default_image=SAIGE_DOCKER_IMAGE,
                 # default_memory='1Gi',
                 default_storage='500Mi', default_cpu=n_threads)

    sparse_grm = None
    analysis_type = "variant" if args.single_variant_only else "gene"
    if analysis_type == "gene":
        if args.create_sparse_grm:
            sparse_grm = create_sparse_grm(p, sparse_grm_root, ukb_for_grm_plink_path,
                                           relatedness_cutoff, num_markers, n_threads=n_threads)
        else:
            sparse_grm = p.read_input_group(**{ext: f'{sparse_grm_root}.{ext}' for ext in
                                               (sparse_grm_extension, f'{sparse_grm_extension}.sampleIDs.txt')})

    pheno_export_dir = f'{root}/pheno_export_data'
    if not args.overwrite_pheno_data:
        phenos_already_exported = {x['path'] for x in hl.hadoop_ls(pheno_export_dir)}
    pheno_exports = {}
    for pheno in phenos_to_run:
        pheno_export_path = f'{pheno_export_dir}/{pheno}.tsv'
        if not args.overwrite_pheno_data and pheno_export_path in phenos_already_exported:
            pheno_file = p.read_input(pheno_export_path)
        else:
            pheno_file = export_pheno(p, pheno_export_path, pheno, get_ukb_pheno_mt_path(trait_type),
                                      get_ukb_covariates_ht_path(), data_type=trait_type)
        pheno_exports[pheno] = pheno_file
    completed = Counter([isinstance(x, pipeline.resource.InputResourceFile) for x in pheno_exports.values()])
    logger.info(f'Exporting {completed[False]} phenos (already found {completed[True]})...')


    overwrite_null_models = args.create_null_models
    null_model_dir = f'{root}/null_glmm'
    if not overwrite_null_models:
        null_models_already_created = {x['path'] for x in hl.hadoop_ls(null_model_dir)}
    null_models = {}

    for pheno in phenos_to_run:
        sparse_sigma_file = None
        null_glmm_root = f'{null_model_dir}/{pheno}'
        model_file_path = f'{null_glmm_root}.rda'
        variance_ratio_file_path = f'{null_glmm_root}.{analysis_type}.varianceRatio.txt'
        sparse_sigma_file_path = variance_ratio_file_path + sparse_grm_extension.replace("GRM", "Sigma")

        if not overwrite_null_models and model_file_path in null_models_already_created and \
                variance_ratio_file_path in null_models_already_created:
            model_file = p.read_input(model_file_path)
            variance_ratio_file = p.read_input(variance_ratio_file_path)
            if analysis_type == 'gene': sparse_sigma_file = p.read_input(sparse_sigma_file_path)
        else:
            if args.skip_any_null_models: continue
            fit_null_task = fit_null_glmm(p, null_glmm_root, pheno_exports[pheno], pheno, trait_type, covariates,
                                          ukb_for_grm_plink_path,
                                          sparse_grm=sparse_grm, sparse_grm_extension=sparse_grm_extension,
                                          n_threads=n_threads)
            model_file = fit_null_task.null_glmm.rda
            variance_ratio_file = fit_null_task.null_glmm[f'{analysis_type}.varianceRatio.txt']
            if analysis_type == 'gene': sparse_sigma_file = fit_null_task.null_glmm[
                f'{analysis_type}.varianceRatio.txt{sparse_grm_extension.replace("GRM", "Sigma")}']
        null_models[pheno] = (model_file, variance_ratio_file, sparse_sigma_file)

    completed = Counter([type(x[0]) == pipeline.resource.InputResourceFile for x in null_models.values()])
    logger.info(f'Running {completed[False]} null models (already found {completed[True]})...')

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
            # if args.local_test:
            #     break
        if args.local_test:
            break

    completed = Counter([type(x[1]) == pipeline.resource.InputResourceFile for x in vcfs.values()])
    logger.info(f'Creating {completed[False]} VCFs (already found {completed[True]})...')

    result_dir = f'{root}/result'
    overwrite_results = args.overwrite_results
    for i, pheno in enumerate(phenos_to_run):
        if pheno not in null_models: continue
        if not i % 10:
            n_jobs = dict(Counter(map(lambda x: x.name, p.select_tasks("")))).get("run_saige", 0)
            logger.info(f'Read {i} phenotypes ({n_jobs} new to run so far)...')

        pheno_results_dir = f'{result_dir}/{pheno}'
        results_already_created = {}
        if not overwrite_results and not args.skip_saige and hl.hadoop_exists(pheno_results_dir):
            results_already_created = {x['path'] for x in hl.hadoop_ls(pheno_results_dir)}

        model_file, variance_ratio_file, sparse_sigma_file = null_models[pheno]
        saige_tasks = []
        for chrom in chromosomes:
            if args.skip_saige: break
            chromosome = f'chr{chrom}'
            chrom_length = chrom_lengths[chromosome]
            for start_pos in range(1, chrom_length, chunk_size):
                end_pos = chrom_length if start_pos + chunk_size > chrom_length else (start_pos + chunk_size)
                interval = f'{chromosome}:{start_pos}-{end_pos}'
                vcf_file, group_file = vcfs[interval]
                results_path = f'{pheno_results_dir}/result_{pheno}_{chromosome}_{str(start_pos).zfill(9)}'
                if overwrite_results or f'{results_path}.{analysis_type}.txt' not in results_already_created:
                    samples_file = p.read_input(get_ukb_samples_file_path())
                    saige_task = run_saige(p, results_path, model_file, variance_ratio_file, vcf_file, samples_file,
                                           group_file, sparse_sigma_file, trait_type=trait_type, use_bgen=use_bgen,
                                           chrom=chromosome)
                    saige_tasks.append(saige_task)
                # if args.local_test:
                #     break
            if args.local_test:
                break
        load_results_into_hail(p, pheno_results_dir, pheno, saige_tasks)
        if args.limit and n_jobs >= args.limit:
            break

    logger.info(f'Setup took: {time.strftime("%H:%M:%S", time.gmtime(time.time() - start_time))}')
    logger.info(f'Submitting: {get_tasks_from_pipeline(p)}')
    p.run(dry_run=args.dry_run, verbose=True, delete_scratch_on_exit=False)
    logger.info(f'Finished: {get_tasks_from_pipeline(p)}')


if __name__ == '__main__':
    parser = argparse.ArgumentParser()

    parser.add_argument('--create_sparse_grm', help='Run single variant SAIGE', action='store_true')
    parser.add_argument('--overwrite_pheno_data', help='Run single variant SAIGE', action='store_true')
    parser.add_argument('--skip_any_null_models', help='Run single variant SAIGE', action='store_true')
    parser.add_argument('--skip_saige', help='Run single variant SAIGE', action='store_true')
    parser.add_argument('--create_null_models', help='Run single variant SAIGE', action='store_true')
    parser.add_argument('--create_vcfs', help='Run single variant SAIGE', action='store_true')
    parser.add_argument('--overwrite_results', help='Run single variant SAIGE', action='store_true')
    parser.add_argument('--single_variant_only', help='Run single variant SAIGE', action='store_true')
    parser.add_argument('--local_test', help='Run single variant SAIGE', action='store_true')
    parser.add_argument('--use_bgen', help='Run single variant SAIGE', action='store_true')
    parser.add_argument('--trait_type', help='Run single variant SAIGE')
    parser.add_argument('--limit', help='Run single variant SAIGE', type=int)
    parser.add_argument('--run_all_phenos', help='Dry run only', action='store_true')
    parser.add_argument('--dry_run', help='Dry run only', action='store_true')
    parser.add_argument('--send_slack', help='Dry run only', action='store_true')
    args = parser.parse_args()

    if args.local_test:
        try_slack('@konradjk', main, args)
    else:
        main(args)


