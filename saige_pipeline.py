#!/usr/bin/env python3

__author__ = 'konradk'

from gnomad_hail import *
import ukbb_qc.resources as ukb

import sys
sys.path.extend(['/home/hail/hail/pipeline', '/home/hail/hail/batch'])
from pipeline import *

bucket = 'gs://ukbb-pharma-exome-analysis'
root = f'{bucket}/data'
common_high_callrate_pruned_plink_path = f'{root}/ukb.common.high.callrate.pruned.plink'
pheno_path = 'gs://phenotype_pharma/PHESANT_output/100k/January_2019/phesant_output_combined_both_sexes_no_sex_specific.tsv'


def fit_null_glmm(p, pheno_name, trait_type, covariates, n_threads, output_root):
    pheno_file = p.read_input(pheno_path)
    in_bfile = p.read_input_group(
        **{ext: f'{common_high_callrate_pruned_plink_path}.{ext}' for ext in ('bed', 'bim', 'fam')})
    bim_fix_task = p.new_task().label('bim_fix')
    bim_fix_task.command(f'perl -pi -e s/^chr// {in_bfile.bim}')
    fit_null_task = p.new_task().label('fit_null_model')
    fit_null_task.depends_on(bim_fix_task)
    fit_null_task.declare_resource_group(
        null_glmm={ext: f'{{root}}{ext if ext.startswith("_") else "." + ext}' for ext in
                   ('rda', '_30markers.SAIGE.results.txt', 'varianceRatio.txt')})
    command1 = (f'Rscript /usr/local/bin/step1_fitNULLGLMM.R '
                f'--plinkFile={in_bfile} '
                f'--phenoFile={pheno_file} '
                f'--covarColList={covariates} '
                f'--phenoCol={pheno_name} '
                f'--sampleIDColinphenoFile=userID '
                f'--traitType={trait_type} '
                f'--outputPrefix={fit_null_task.null_glmm} '
                f'--nThreads={n_threads} --LOCO=FALSE &> {fit_null_task.stdout}')
    fit_null_task.command(command1)
    p.write_output(fit_null_task.null_glmm, output_root)
    p.write_output(fit_null_task.stdout, f'{output_root}.log')
    # Runtimes: 8 threads: ~5 minutes of 100% CPU (~3G RAM), followed by ~9 minutes of 800% (~6G RAM)
    return fit_null_task


def run_saige(p, model_file, run_saige_task, variance_ratio_file, chrom: str = 'chr1',
              min_mac: int = 1, min_maf: float = 0.0001):
    samples_file = p.read_input(f'{common_high_callrate_pruned_plink_path}.samples')
    vcf_file = p.read_input_group(**{'vcf.gz': f'{common_high_callrate_pruned_plink_path}.vcf.bgz',
                                     'vcf.gz.tbi': f'{common_high_callrate_pruned_plink_path}.vcf.bgz.tbi'})
    touch_index_task = p.new_task()
    touch_index_task.command(f'touch {vcf_file["vcf.gz"]}; touch {vcf_file["vcf.gz.tbi"]}')
    run_saige_task.depends_on(touch_index_task)
    run_saige_task.declare_resource_group(result={'txt': '{root}.txt'})
    command2 = (f'Rscript /usr/local/bin/step2_SPAtests.R '
                f'--vcfFile={vcf_file["vcf.gz"]} '
                f'--vcfFileIndex={vcf_file["vcf.gz.tbi"]} '
                f'--vcfField=GT '
                f'--minMAF={min_maf} '
                f'--minMAC={min_mac} '
                f'--chrom={chrom} '
                f'--sampleFile={samples_file} '
                f'--GMMATmodelFile={model_file} '
                f'--varianceRatioFile={variance_ratio_file} '
                f'--SAIGEOutputFile={run_saige_task.result.txt} '
                f'--IsOutputAFinCaseCtrl=TRUE &> {run_saige_task.stdout}')
    run_saige_task.command(command2)
    p.write_output(run_saige_task.result.txt, f'{root}/tmp/result')
    p.write_output(run_saige_task.stdout, f'{root}/tmp/result.log')


def main(args):
    trait_type = 'quantitative'
    n_threads = 8
    pheno_name = 'X50_irnt'
    covariates = 'X46_irnt'
    null_glmm_root = f'{root}/tmp/test'

    local_backend = LocalBackend()
    p = Pipeline(backend=local_backend, default_image='wzhou88/saige:0.35.8.2')

    # tabix_task = p.new_task().label('tabix')
    # tabix_task.command(f'tabix {vcf_file}')

    step0_done = hl.hadoop_exists(f'{null_glmm_root}.rda')

    if not step0_done:
        fit_null_task = fit_null_glmm(p, pheno_name, trait_type, covariates, n_threads, null_glmm_root)

    run_saige_task = p.new_task().label('run_saige')
    if step0_done:
        model_file = p.read_input(f'{null_glmm_root}.rda')
        variance_ratio_file = p.read_input(f'{null_glmm_root}.varianceRatio.txt')
    else:
        run_saige_task.depends_on(fit_null_task)#, tabix_task)
        model_file = fit_null_task.null_glmm.rda
        variance_ratio_file = fit_null_task.null_glmm["varianceRatio.txt"]
    run_saige(p, model_file, run_saige_task, variance_ratio_file)
    p.run(verbose=True, delete_scratch_on_exit=False)


if __name__ == '__main__':
    parser = argparse.ArgumentParser()

    parser.add_argument('--overwrite', help='Overwrite everything', action='store_true')
    parser.add_argument('--slack_channel', help='Send message to Slack channel/user', default='@konradjk')
    args = parser.parse_args()

    if args.slack_channel:
        try_slack(args.slack_channel, main, args)
    else:
        main(args)


