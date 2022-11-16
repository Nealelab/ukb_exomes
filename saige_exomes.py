#!/usr/bin/env python3

__author__ = 'konradk'

import time

import re
import argparse
from gnomad.utils import slack
from ukb_common import *
from ukb_common.utils.saige_pipeline import *
from ukb_exomes import *

root = f'{bucket}/{CURRENT_TRANCHE}/results'


HAIL_DOCKER_IMAGE = 'gcr.io/ukbb-exome-pharma/hail_utils:6.5'
SAIGE_DOCKER_IMAGE = 'gcr.io/ukbb-exome-pharma/saige:0.3'
QQ_DOCKER_IMAGE = 'gcr.io/ukbb-exome-pharma/saige_qq:0.1'


def get_exclusions_200k():
    exclusions = set()
    with hl.hadoop_open(f'{root}/exclusions_200k.txt') as f:
        for line in f:
            pheno = line.strip().split('\t')
            if len(pheno) == 2:
                pheno, coding = pheno
            else:
                pheno = pheno[0]
                coding = ''
            exclusions.add((pheno, coding))
    return exclusions


def get_phenos_to_run(limit: int = None, pilot: bool = False, specific_phenos: str = None, min_cases: int = None):
    ht = hl.read_matrix_table(get_ukb_pheno_mt_path()).cols()
    exclusion_ht = hl.import_table(f'gs://ukbb-pharma-exome-analysis/500k/450K_phentoypes_removeflag_vF04132021.txt',
                                   impute=True, key=PHENO_KEY_FIELDS)
    ht = ht.annotate(keep=(exclusion_ht[ht.key]['status.remove'] != 'remove') | (ht.phenocode == 'COVID19') | ht.phenocode.startswith('random_'))
    if not min_cases:
        min_cases = MIN_CASES
    cancer_phenos = hl.literal({"C15", "C16", "C17", "C18", "C187", "C19", "C20", "C21", "C22", "C25", "C32", "C34",
                                "C43", "C44", "C50", "C51", "C53", "C54", "C541", "C56", "C61", "C62", "C64", "C67",
                                "C71", "C73", "C77", "C78", "C79", "C80", "C81", "C82", "C83", "C84", "C85", "C88",
                                "C90", "C91", "C92", "D05", "D06", "D25", "Z08", "Z09", "Z12", "Z80", "Z85"})
    ht = ht.filter(
        (
                hl.set({'biogen', 'abbvie', 'pfizer'}).contains(ht.modifier) |
                ht.phenocode.endswith('_pfe') |
                cancer_phenos.contains(ht.phenocode) |
                (
                        (ht.n_cases_defined >= min_cases) &
                        ((ht.score >= 1) | (ht.trait_type == 'icd_first_occurrence') |
                         (ht.phenocode == 'COVID19') |
                         (ht.phenocode.startswith('30') & hl.is_missing(ht.score))
                         ) &  # biomarkers
                        ht.keep & (ht.modifier != 'raw') &
                        ~ht.description.contains('Source of report')
                )
        ) & hl.is_defined(ht.modifier)
    )
    output = set([tuple(x[field] for field in PHENO_KEY_FIELDS) for x in ht.select().collect()])

    if pilot:
        output = output.intersection(PILOT_PHENOTYPES)
    if specific_phenos:
        specific_phenos = specific_phenos.split(',')
        output = [x for x in output if all(map(lambda y: y is not None, x)) and any([re.match(pcd, '-'.join(x)) for pcd in specific_phenos])]
    if limit:
        output = set(sorted(output)[:limit])
    pheno_key_dict = [dict(zip(PHENO_KEY_FIELDS, x)) for x in output]
    return pheno_key_dict


def main(args):
    hl.init(log='/tmp/saige_temp_hail.log')

    phenos_to_run = get_phenos_to_run(1 if args.local_test else None, pilot=not args.run_all_phenos,
                                      specific_phenos=args.phenos, min_cases=args.min_cases)
    logger.info(f'Got {len(phenos_to_run)} phenotypes...')
    if len(phenos_to_run) <= 20:
        logger.info(phenos_to_run)

    num_pcs = 20
    start_time = time.time()
    basic_covars = ['sex', 'age', 'age2', 'age_sex', 'age2_sex', 'batch_150K', 'batch_200K', 'batch_300K']
    covariates = ','.join(basic_covars + [f'pc{x}' for x in range(1, num_pcs + 1)])
    relatedness_cutoff = '0.125'
    num_markers = 2000
    n_threads = 8
    chromosomes = ['X', 'Y']  # list(range(1, 23)) + ['X', 'Y']

    sparse_grm_root = f'{root}/grm/sparse'
    sparse_grm_extension = f'_relatednessCutoff_{relatedness_cutoff}_{num_markers}_randomMarkersUsed.sparseGRM.mtx'

    if args.local_test:
        backend = hb.LocalBackend(gsa_key_file='/Users/konradk/.hail/ukb_exomes.json')
    else:
        backend = hb.ServiceBackend(billing_project='ukb_pharma', bucket='ukb-pharma-exome-analysis-temp')
    p = hb.Batch(name='saige_exomes', backend=backend, default_image=SAIGE_DOCKER_IMAGE,
                 # default_memory='1Gi',
                 default_cpu=n_threads)

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
    for pheno_key_dict in phenos_to_run:
        pheno_export_path = get_pheno_output_path(pheno_export_dir, pheno_key_dict)
        if not args.overwrite_pheno_data and pheno_export_path in phenos_already_exported:
            pheno_file = p.read_input(pheno_export_path)
        else:
            pheno_task = export_pheno(p, pheno_export_path, pheno_key_dict, 'ukb_exomes', HAIL_DOCKER_IMAGE)
            pheno_file = pheno_task.out
        pheno_exports[stringify_pheno_key_dict(pheno_key_dict)] = pheno_file
    completed = Counter([isinstance(x, InputResourceFile) for x in pheno_exports.values()])
    logger.info(f'Exporting {completed[False]} phenos (already found {completed[True]})...')


    overwrite_null_models = args.create_null_models
    null_model_dir = f'{root}/null_glmm'
    null_models_already_created = {}
    if not overwrite_null_models and hl.hadoop_exists(null_model_dir):
        null_models_already_created = {x['path'] for x in hl.hadoop_ls(null_model_dir)}
    null_models = {}

    for pheno_key_dict in phenos_to_run:
        null_glmm_root = get_pheno_output_path(null_model_dir, pheno_key_dict, '')
        sparse_sigma_file = None
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
            fit_null_task = fit_null_glmm(p, null_glmm_root, pheno_exports[stringify_pheno_key_dict(pheno_key_dict)],
                                          pheno_key_dict['trait_type'], covariates,
                                          get_ukb_grm_plink_path(), SAIGE_DOCKER_IMAGE,
                                          sparse_grm=sparse_grm, sparse_grm_extension=sparse_grm_extension,
                                          n_threads=n_threads)
            fit_null_task.attributes.update(copy.deepcopy(pheno_key_dict))
            model_file = fit_null_task.null_glmm.rda
            variance_ratio_file = fit_null_task.null_glmm[f'{analysis_type}.varianceRatio.txt']
            if analysis_type == 'gene': sparse_sigma_file = fit_null_task.null_glmm[
                f'{analysis_type}.varianceRatio.txt{sparse_grm_extension.replace("GRM", "Sigma")}']
        null_models[stringify_pheno_key_dict(pheno_key_dict)] = (model_file, variance_ratio_file, sparse_sigma_file)

    completed = Counter([type(x[0]) == InputResourceFile for x in null_models.values()])
    logger.info(f'Running {completed[False]} null models (already found {completed[True]})...')

    chrom_lengths = hl.get_reference('GRCh38').lengths

    use_bgen = True  # args.use_bgen
    group_file_only = True
    common_variants_only = True
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
            interval, long_interval = get_intervals(chrom_length, chromosome, chunk_size, start_pos)
            vcf_root = f'{vcf_dir}/test_{chromosome}_{str(start_pos).zfill(9)}'
            if not overwrite_vcfs and f'{vcf_root}.{test_extension}' in vcfs_already_created or args.skip_vcfs:
                if f'{vcf_root}.{test_extension}' not in vcfs_already_created:
                    continue
                vcf_file = get_vcf_file(p, use_bgen, vcf_root)
                group_file = p.read_input(f'{vcf_root}.gene.txt')
            else:
                vcf_task = extract_vcf_from_mt(p, vcf_root, HAIL_DOCKER_IMAGE, interval=long_interval,
                                               gene_ht_interval=interval, export_bgen=use_bgen,
                                               n_threads=n_threads, gene_map_ht_path=get_ukb_gene_map_ht_path(),
                                               additional_args=interval, group_file_only=group_file_only,
                                               common_variants_only=common_variants_only)
                if group_file_only:
                    vcf_file = get_vcf_file(p, use_bgen, vcf_root)
                else:
                    vcf_file = vcf_task.out
                group_file = vcf_task.group_file
            vcfs[interval] = (vcf_file, group_file)
            if args.local_test:
                break
        if args.local_test:
            break

    completed = Counter([type(x[1]) == InputResourceFile for x in vcfs.values()])
    logger.info(f'Creating {completed[False]} VCFs (already found {completed[True]})...')

    result_dir = f'{root}/result'
    overwrite_results = args.overwrite_results
    for i, pheno_key_dict in enumerate(phenos_to_run):
        if stringify_pheno_key_dict(pheno_key_dict) not in null_models: continue
        model_file, variance_ratio_file, sparse_sigma_file = null_models[stringify_pheno_key_dict(pheno_key_dict)]

        if not i % 10:
            n_jobs = dict(Counter(map(lambda x: x.name, p.select_jobs("")))).get("run_saige", 1)
            logger.info(f'Read {i} phenotypes ({n_jobs} new to run so far)...')

        pheno_results_dir = get_pheno_output_path(result_dir, pheno_key_dict, '')
        results_already_created = {}
        if not overwrite_results and not args.skip_saige and hl.hadoop_exists(pheno_results_dir):
            results_already_created = {x['path'] for x in hl.hadoop_ls(pheno_results_dir)}

        saige_tasks = []
        for chrom in chromosomes:
            if args.skip_saige: break
            chromosome = f'chr{chrom}'
            chrom_length = chrom_lengths[chromosome]
            for start_pos in range(1, chrom_length, chunk_size):
                interval, _ = get_intervals(chrom_length, chromosome, chunk_size, start_pos)
                if interval not in vcfs: continue
                vcf_file, group_file = vcfs[interval]
                results_path = get_results_prefix(pheno_results_dir, pheno_key_dict, chromosome, start_pos)
                if overwrite_results or f'{results_path}.{args.add_suffix}{analysis_type}.txt' not in results_already_created:
                    samples_file = p.read_input(get_ukb_samples_file_path())
                    saige_task = run_saige(p, results_path, model_file, variance_ratio_file, vcf_file, samples_file,
                                           SAIGE_DOCKER_IMAGE, group_file=group_file,
                                           sparse_sigma_file=sparse_sigma_file, trait_type=pheno_key_dict['trait_type'], use_bgen=use_bgen,
                                           chrom=chromosome, max_maf=0.999, add_suffix=args.add_suffix, min_maf=0.04)
                    saige_task.attributes.update(copy.deepcopy(pheno_key_dict))
                    saige_tasks.append(saige_task)
                if args.local_test:
                    break
            if args.local_test:
                break
        res_tasks = []
        if overwrite_results or args.overwrite_hail_results or \
                f'{pheno_results_dir}/variant_results.mt' not in results_already_created or \
                not hl.hadoop_exists(f'{pheno_results_dir}/variant_results.mt/_SUCCESS'):
            null_glmm_root = get_pheno_output_path(null_model_dir, pheno_key_dict, f'.{analysis_type}.log')

            prefix = get_results_prefix(pheno_results_dir, pheno_key_dict, 'chr{chrom}', 1)
            saige_log = f'{prefix}.{analysis_type}.log'

            load_task = load_results_into_hail(p, pheno_results_dir, pheno_key_dict,
                                               saige_tasks, get_ukb_vep_path(), HAIL_DOCKER_IMAGE,
                                               saige_log=saige_log,
                                               analysis_type=analysis_type,
                                               gene_map_path=get_ukb_gene_map_ht_path(post_processed=False),
                                               n_threads=n_threads, null_glmm_log=null_glmm_root)
            res_tasks.append(load_task)
            if not args.skip_qq:
                qq_export, qq_plot = qq_plot_results(p, pheno_results_dir, res_tasks, HAIL_DOCKER_IMAGE, QQ_DOCKER_IMAGE,
                                                     n_threads=n_threads)
                qq_export.attributes.update(copy.deepcopy(pheno_key_dict))
                qq_plot.attributes.update(copy.deepcopy(pheno_key_dict))

        if args.limit and n_jobs >= args.limit:
            break

    logger.info(f'Setup took: {time.strftime("%H:%M:%S", time.gmtime(time.time() - start_time))}')
    logger.info(f'Submitting: {get_tasks_from_pipeline(p)}')
    logger.info(f"Total size: {sum([len(x._pretty()) for x in p.select_jobs('')])}")
    p.run(dry_run=args.dry_run, delete_scratch_on_exit=False)
    logger.info(f'Finished: {get_tasks_from_pipeline(p)}')


def get_vcf_file(p, use_bgen, vcf_root):
    if use_bgen:
        vcf_file = p.read_input_group(**{'bgen': f'{vcf_root}.bgen',
                                         'bgen.bgi': f'{vcf_root}.bgen.bgi',
                                         'sample': f'{vcf_root}.sample'})
    else:
        vcf_file = p.read_input_group(**{'vcf.gz': f'{vcf_root}.vcf.gz',
                                         'vcf.gz.tbi': f'{vcf_root}.vcf.gz.tbi'})
    return vcf_file


def get_intervals(chrom_length, chromosome, chunk_size, start_pos):
    end_pos = chrom_length if start_pos + chunk_size > chrom_length else (start_pos + chunk_size)
    short_interval = f'{chromosome}:{start_pos}-{end_pos}'
    padded_chunk_size = chunk_size * 3
    end_pos = chrom_length if start_pos + padded_chunk_size > chrom_length else (start_pos + padded_chunk_size)
    long_interval = f'{chromosome}:{start_pos}-{end_pos}'
    return short_interval, long_interval


if __name__ == '__main__':
    parser = argparse.ArgumentParser()

    parser.add_argument('--create_sparse_grm', help='Run single variant SAIGE', action='store_true')
    parser.add_argument('--overwrite_pheno_data', help='Run single variant SAIGE', action='store_true')
    parser.add_argument('--skip_any_null_models', help='Run single variant SAIGE', action='store_true')
    parser.add_argument('--skip_saige', help='Run single variant SAIGE', action='store_true')
    parser.add_argument('--create_null_models', help='Run single variant SAIGE', action='store_true')
    parser.add_argument('--create_vcfs', help='Run single variant SAIGE', action='store_true')
    parser.add_argument('--skip_vcfs', help='Run single variant SAIGE', action='store_true')
    parser.add_argument('--overwrite_results', help='Run single variant SAIGE', action='store_true')
    parser.add_argument('--overwrite_hail_results', help='Run single variant SAIGE', action='store_true')
    parser.add_argument('--skip_qq', help='Run single variant SAIGE', action='store_true')
    parser.add_argument('--single_variant_only', help='Run single variant SAIGE', action='store_true')
    parser.add_argument('--local_test', help='Run single variant SAIGE', action='store_true')
    parser.add_argument('--use_bgen', help='Run single variant SAIGE', action='store_true')
    parser.add_argument('--limit', help='Run single variant SAIGE', type=int)
    parser.add_argument('--phenos', help='Comma-separated list of trait_type-phenocode-pheno_sex-coding-modifier regexes '
                                         '(e.g. continuous-50-both_sexes--,icd10-E1.*,brain_mri-.* )')
    parser.add_argument('--add_suffix', default='')
    parser.add_argument('--min_cases', help='Dry run only', type=int)
    parser.add_argument('--run_all_phenos', help='Dry run only', action='store_true')
    parser.add_argument('--dry_run', help='Dry run only', action='store_true')
    parser.add_argument('--send_slack', help='Dry run only', action='store_true')
    args = parser.parse_args()

    if args.local_test:
        from slack_token_pkg.slack_creds import slack_token
        with slack.slack_notifications(slack_token, '@konradjk'):
            main(args)
    else:
        main(args)


