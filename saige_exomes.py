#!/usr/bin/env python3

__author__ = 'konradk'

import time

from ukb_common import *
from ukb_common.utils.saige_pipeline import *
from ukb_exomes import *

bucket = 'gs://ukbb-pharma-exome-analysis'
root = f'{bucket}/{CURRENT_TRANCHE}/results'


HAIL_DOCKER_IMAGE = 'gcr.io/ukbb-exome-pharma/hail_utils:3.5'
SAIGE_DOCKER_IMAGE = 'wzhou88/saige:0.36.3'
QQ_DOCKER_IMAGE = 'konradjk/saige_qq:0.2'


def get_phenos_to_run(sex: str, limit: int = None, pilot: bool = False):
    ht = hl.read_matrix_table(get_ukb_pheno_mt_path()).cols()
    ht = ht.filter((ht[f'n_cases_{sex}'] >= MIN_CASES) &
                   ((ht.score >= 1) | (ht.data_type == 'icd'))).select('data_type')
    output = {(x.pheno, x.coding, x.data_type) for x in ht.collect()}

    if pilot:
        output = output.intersection(PILOT_PHENOTYPES)
    if limit:
        output = set(sorted(output)[:limit])
    return output


def main(args):
    hl.init(log='/tmp/saige_temp_hail.log')
    sex = 'both_sexes'

    phenos_to_run = get_phenos_to_run(sex, 1 if args.local_test else None, pilot=not args.run_all_phenos)
    logger.info(f'Got {len(phenos_to_run)} phenotypes...')
    if len(phenos_to_run) <= 20:
        logger.info(phenos_to_run)

    num_pcs = 20
    start_time = time.time()
    basic_covars = ['sex', 'age', 'age2', 'age_sex', 'age2_sex', 'batch_1', 'batch_2']
    covariates = ','.join(basic_covars + [f'pc{x}' for x in range(1, num_pcs + 1)])
    relatedness_cutoff = '0.125'
    num_markers = 2000
    n_threads = 8
    chromosomes = list(range(1, 23)) + ['X', 'Y']

    sparse_grm_root = f'{root}/grm/sparse'
    sparse_grm_extension = f'_relatednessCutoff_{relatedness_cutoff}_{num_markers}_randomMarkersUsed.sparseGRM.mtx'

    if args.local_test:
        backend = pipeline.LocalBackend(gsa_key_file='/Users/konradk/.hail/ukb_exomes.json')
    else:
        backend = pipeline.BatchBackend(billing_project='ukb_pharma')
    p = pipeline.Pipeline(name='saige_exomes', backend=backend, default_image=SAIGE_DOCKER_IMAGE,
                 # default_memory='1Gi',
                 default_storage='500Mi', default_cpu=n_threads)

    sparse_grm = None
    analysis_type = "variant" if args.single_variant_only else "gene"
    if analysis_type == "gene":
        if args.create_sparse_grm:
            sparse_grm = create_sparse_grm(p, sparse_grm_root, get_ukb_grm_plink_path(), SAIGE_DOCKER_IMAGE,
                                           relatedness_cutoff, num_markers, n_threads=n_threads)
        else:
            sparse_grm = p.read_input_group(**{ext: f'{sparse_grm_root}.{ext}' for ext in
                                               (sparse_grm_extension, f'{sparse_grm_extension}.sampleIDs.txt')})

    pheno_export_dir = f'{root}/pheno_export_data'
    phenos_already_exported = {}
    if not args.overwrite_pheno_data and hl.hadoop_exists(pheno_export_dir):
        phenos_already_exported = {x['path'] for x in hl.hadoop_ls(pheno_export_dir)}
    pheno_exports = {}
    for pheno_coding_trait in phenos_to_run:
        pheno, coding, trait_type = pheno_coding_trait
        pheno_export_path = f'{pheno_export_dir}/{trait_type}-{pheno}-{coding}.tsv'
        if not args.overwrite_pheno_data and pheno_export_path in phenos_already_exported:
            pheno_file = p.read_input(pheno_export_path)
        else:
            pheno_task = export_pheno(p, pheno_export_path, pheno, coding, trait_type, 'ukb_exomes', HAIL_DOCKER_IMAGE)
            pheno_file = pheno_task.out
        pheno_exports[pheno_coding_trait] = pheno_file
    completed = Counter([isinstance(x, pipeline.resource.InputResourceFile) for x in pheno_exports.values()])
    logger.info(f'Exporting {completed[False]} phenos (already found {completed[True]})...')


    overwrite_null_models = args.create_null_models
    null_model_dir = f'{root}/null_glmm'
    null_models_already_created = {}
    if not overwrite_null_models and hl.hadoop_exists(null_model_dir):
        null_models_already_created = {x['path'] for x in hl.hadoop_ls(null_model_dir)}
    null_models = {}

    for pheno_coding_trait in phenos_to_run:
        pheno, coding, trait_type = pheno_coding_trait
        sparse_sigma_file = None
        null_glmm_root = f'{null_model_dir}/{trait_type}-{pheno}-{coding}'
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
            fit_null_task = fit_null_glmm(p, null_glmm_root, pheno_exports[pheno_coding_trait], trait_type, covariates,
                                          get_ukb_grm_plink_path(), SAIGE_DOCKER_IMAGE,
                                          sparse_grm=sparse_grm, sparse_grm_extension=sparse_grm_extension,
                                          n_threads=n_threads)
            fit_null_task.attributes['pheno'] = pheno
            model_file = fit_null_task.null_glmm.rda
            variance_ratio_file = fit_null_task.null_glmm[f'{analysis_type}.varianceRatio.txt']
            if analysis_type == 'gene': sparse_sigma_file = fit_null_task.null_glmm[
                f'{analysis_type}.varianceRatio.txt{sparse_grm_extension.replace("GRM", "Sigma")}']
        null_models[pheno_coding_trait] = (model_file, variance_ratio_file, sparse_sigma_file)

    completed = Counter([type(x[0]) == pipeline.resource.InputResourceFile for x in null_models.values()])
    logger.info(f'Running {completed[False]} null models (already found {completed[True]})...')

    chrom_lengths = hl.get_reference('GRCh38').lengths

    use_bgen = args.use_bgen
    vcf_dir = f'{root}/vcf'
    test_extension = 'bgen' if use_bgen else 'vcf.gz'
    overwrite_vcfs = args.create_vcfs
    vcfs_already_created = {}
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
            if f'{vcf_root}.{test_extension}' in vcfs_already_created:
                if use_bgen:
                    vcf_file = p.read_input_group(**{'bgen': f'{vcf_root}.bgen',
                                                     'bgen.bgi': f'{vcf_root}.bgen.bgi',
                                                     'sample': f'{vcf_root}.sample'})
                else:
                    vcf_file = p.read_input_group(**{'vcf.gz': f'{vcf_root}.vcf.gz',
                                                     'vcf.gz.tbi': f'{vcf_root}.vcf.gz.tbi'})
                group_file = p.read_input(f'{vcf_root}.gene.txt')
            else:
                vcf_task = extract_vcf_from_mt(p, vcf_root, HAIL_DOCKER_IMAGE, interval=interval, export_bgen=use_bgen,
                                               gene_map_ht_path=get_ukb_gene_map_ht_path())
                vcf_file = vcf_task.out
                group_file = vcf_task.group_file
            vcfs[interval] = (vcf_file, group_file)
            if args.local_test:
                break
        if args.local_test:
            break

    completed = Counter([type(x[1]) == pipeline.resource.InputResourceFile for x in vcfs.values()])
    logger.info(f'Creating {completed[False]} VCFs (already found {completed[True]})...')

    result_dir = f'{root}/result'
    overwrite_results = args.overwrite_results
    for i, pheno_coding_trait in enumerate(phenos_to_run):
        if pheno_coding_trait not in null_models: continue
        pheno, coding, trait_type = pheno_coding_trait
        if not i % 10:
            n_jobs = dict(Counter(map(lambda x: x.name, p.select_tasks("")))).get("run_saige", 1)
            logger.info(f'Read {i} phenotypes ({n_jobs} new to run so far)...')

        pheno_results_dir = f'{result_dir}/{trait_type}-{pheno}-{coding}'
        results_already_created = {}
        if not overwrite_results and not args.skip_saige and hl.hadoop_exists(pheno_results_dir):
            results_already_created = {x['path'] for x in hl.hadoop_ls(pheno_results_dir)}

        model_file, variance_ratio_file, sparse_sigma_file = null_models[pheno_coding_trait]
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
                                           SAIGE_DOCKER_IMAGE, group_file=group_file,
                                           sparse_sigma_file=sparse_sigma_file, trait_type=trait_type, use_bgen=use_bgen,
                                           chrom=chromosome)
                    saige_task.attributes.update({'interval': interval, 'pheno': pheno, 'coding': coding, 'trait_type': trait_type})
                    saige_tasks.append(saige_task)
                if args.local_test:
                    break
            if args.local_test:
                break
        res_tasks = []
        if overwrite_results or \
                f'{pheno_results_dir}/variant_results.mt' not in results_already_created or \
                not hl.hadoop_exists(f'{pheno_results_dir}/variant_results.mt/_SUCCESS'):
            res_tasks.append(load_results_into_hail(p, pheno_results_dir, pheno, coding, trait_type,
                                                    saige_tasks, get_ukb_vep_path(), HAIL_DOCKER_IMAGE,
                                                    analysis_type=analysis_type,
                                                    gene_map_path=get_ukb_gene_map_ht_path(post_processed=False)))
            qq_plot_results(p, pheno_results_dir, res_tasks, HAIL_DOCKER_IMAGE, QQ_DOCKER_IMAGE)
        if args.limit and n_jobs >= args.limit:
            break

    logger.info(f'Setup took: {time.strftime("%H:%M:%S", time.gmtime(time.time() - start_time))}')
    logger.info(f'Submitting: {get_tasks_from_pipeline(p)}')
    logger.info(f"Total size: {sum([len(x._pretty()) for x in p.select_tasks('')])}")
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
    parser.add_argument('--limit', help='Run single variant SAIGE', type=int)
    parser.add_argument('--run_all_phenos', help='Dry run only', action='store_true')
    parser.add_argument('--dry_run', help='Dry run only', action='store_true')
    parser.add_argument('--send_slack', help='Dry run only', action='store_true')
    args = parser.parse_args()

    if args.local_test:
        try_slack('@konradjk', main, args)
    else:
        main(args)


