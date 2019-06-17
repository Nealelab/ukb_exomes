#!/usr/bin/env python3

__author__ = 'konradk'

from gnomad_hail import *
# import ukbb_qc.resources as ukb
from ukb_exomes import *
import shlex

import sys
sys.path.extend(['/Users/konradk/hail/pipeline', '/Users/konradk/hail/batch', '/Users/konradk/hail/hailjwt'])
from pipeline import *
import pipeline.pipeline

bucket = 'gs://ukbb-pharma-exome-analysis'
root = f'{bucket}/data'
ukb_for_grm_plink_path = f'{bucket}/mt/ukb.for_grm.pruned.plink'
pheno_path = 'gs://phenotype_pharma/PHESANT_output/100k/January_2019/phesant_output_combined_both_sexes_no_sex_specific.tsv'

HAIL_DOCKER_IMAGE = 'konradjk/hail_utils:latest'
# SAIGE_DOCKER_IMAGE = 'wzhou88/saige:0.35.8.3'
SAIGE_DOCKER_IMAGE = 'konradjk/saige:0.35.8.2.2'


def create_sparse_grm(p: Pipeline, output_path: str, plink_file_root: str,
                      relatedness_cutoff: str = '0.125', num_markers: int = 2000, n_threads: int = 8):
    in_bfile = p.read_input_group(
        **{ext: f'{plink_file_root}.{ext}' for ext in ('bed', 'bim', 'fam')})
    create_sparse_grm_task: pipeline.pipeline.Task = p.new_task(name='create_sparse_grm')
    create_sparse_grm_task.cpu(n_threads)
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
                        gene: str = None, interval: str = None, groups=None):
    if groups is None:
        groups = {'pLoF', 'missense|LC', 'pLoF|missense|LC', 'synonymous'}
    extract_task: pipeline.pipeline.Task = p.new_task(name='extract_vcf')
    extract_task.image(HAIL_DOCKER_IMAGE)
    extract_task.declare_resource_group(vcf={'vcf.gz': f'{{root}}.vcf.gz',
                                             'vcf.gz.tbi': f'{{root}}.vcf.gz.tbi'})
    python_command = f"""import hail as hl
gene_ht = hl.read_table('{gene_map_ht_path}')"""
    if gene is not None:
        python_command += f"""
gene_ht = gene_ht.filter(gene_ht.gene_symbol == '{gene}')
interval = gene_ht.aggregate(hl.agg.take(gene_ht.interval, 1), _localize=False)[0]"""
    elif interval is not None:
        python_command += f"""
interval = hl.parse_locus_interval('{interval}', 'GRCh38')
gene_ht = gene_ht.filter(interval.contains(gene_ht.interval.start))"""
    python_command += (f"""
gene_ht = gene_ht.filter(hl.set({{'{"', '".join(groups)}'}}).contains(gene_ht.annotation)).key_by()
gene_ht.select(group=gene_ht.gene_id + '_' + gene_ht.gene_symbol + '_' + gene_ht.annotation, variant=hl.delimit(gene_ht.variants, '"\\\\t"')).export('{extract_task.group_file}', header=False)
mt = hl.read_matrix_table('{mt_path}').select_entries('GT')
ht = hl.read_table('gs://broad-ukbb/regeneron.freeze_4/sample_qc/meta_w_pop_adj.ht')
ht = ht.filter(~ht.is_filtered & (ht.hybrid_pop == '12'))
mt = mt.filter_cols(hl.is_defined(ht[mt.col_key]))
mapping = hl.import_table('{sample_mapping_file}', key='eid_sample', delimiter=',')
mt = mt.annotate_cols(exome_id=mapping[mt.s.split('_')[1]].eid_26041)
mt = mt.filter_cols(hl.is_defined(mt.exome_id)).key_cols_by('exome_id')
mt = hl.filter_intervals(mt.annotate_rows(rsid=mt.locus.contig + ':' + hl.str(mt.locus.position) + '_' + mt.alleles[0] + '/' + mt.alleles[1]), [interval])
hl.export_vcf(mt, '{extract_task.bgz}.bgz')""")
    python_command = python_command.replace('\n', '; ').strip()
    command = (f"""
python3 -c "{python_command}";
mv {extract_task.bgz}.bgz {extract_task.vcf["vcf.gz"]};
tabix {extract_task.vcf["vcf.gz"]};
    """)
    extract_task.command(command.replace('\n', ' '))
    p.write_output(extract_task.vcf, output_root)
    p.write_output(extract_task.group_file, f'{output_root}.gene.txt')
    return extract_task


def fit_null_glmm(p: Pipeline, output_root: str, pheno_path: str, pheno_name: str, trait_type: str, covariates: str,
                  plink_file_root: str, sparse_grm: pipeline.pipeline.Resource = None,
                  sparse_grm_extension: str = None, skip_model_fitting: bool = False, n_threads: int = 8):
    analysis_type = "variant" if sparse_grm is None else "gene"

    pheno_file = p.read_input(pheno_path)
    in_bfile = p.read_input_group(**{ext: f'{plink_file_root}.{ext}' for ext in ('bed', 'bim', 'fam')})
    fit_null_task = p.new_task(name=f'fit_null_model_{analysis_type}').cpu(n_threads).storage('1500Mi')
    output_files = {ext: f'{{root}}{ext if ext.startswith("_") else "." + ext}' for ext in
                   ('rda', '_30markers.SAIGE.results.txt', f'{analysis_type}.varianceRatio.txt')}
    if analysis_type == 'gene':
        sparse_sigma_extension = sparse_grm_extension.replace("GRM", "Sigma")
        output_files[f'{analysis_type}.varianceRatio.txt{sparse_sigma_extension}'] = \
            f'{{root}}.{analysis_type}.varianceRatio.txt{sparse_sigma_extension}'
    fit_null_task.declare_resource_group(null_glmm=output_files)
    bim_fix_command = f'perl -pi -e s/^chr// {in_bfile.bim}'

    command = (f'Rscript /usr/local/bin/step1_fitNULLGLMM.R '
               f'--plinkFile={in_bfile} '
               f'--phenoFile={pheno_file} '
               f'--covarColList={covariates} '
               f'--phenoCol={pheno_name} '
               f'--sampleIDColinphenoFile=userID '
               f'--traitType={trait_type} '
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
              vcf_file: pipeline.pipeline.ResourceGroup, samples_file_path: str = None,
              group_file: str = None, sparse_sigma_file: str = None, dependency: pipeline.pipeline.Task = None,
              chrom: str = 'chr1', min_mac: int = 1, min_maf: float = 0, max_maf: float = 0.5):

    analysis_type = "gene" if sparse_sigma_file is not None else "variant"
    samples_file = p.read_input(samples_file_path)

    # Step 2 is single-threaded only
    run_saige_task: pipeline.pipeline.Task = p.new_task(name=f'run_saige_{analysis_type}').cpu(1)
    if dependency is not None:
        run_saige_task.depends_on(dependency)

    command = (f'Rscript /usr/local/bin/step2_SPAtests.R '
               f'--vcfFile={vcf_file["vcf.gz"]} '
               f'--vcfFileIndex={vcf_file["vcf.gz.tbi"]} '
               f'--vcfField=GT '
               f'--minMAF={min_maf} '
               f'--minMAC={min_mac} '
               f'--maxMAFforGroupTest={max_maf} '
               f'--chrom={chrom} '
               f'--sampleFile={samples_file} '
               f'--GMMATmodelFile={model_file} '
               f'--varianceRatioFile={variance_ratio_file} '
               f'--SAIGEOutputFile={run_saige_task.result} ')
    if analysis_type == "gene":
        command += (f'--groupFile={group_file} '
                    f'--sparseSigmaFile={sparse_sigma_file} '
                    f'--IsSingleVarinGroupTest=TRUE ')
    command += f'--IsOutputAFinCaseCtrl=TRUE 2>&1 | tee {run_saige_task.stdout}'
    run_saige_task.command(command)
    p.write_output(run_saige_task.result, f'{output_root}.txt')
    p.write_output(run_saige_task.stdout, f'{output_root}.log')
    # Runtimes: PCSK9 (4 groups, 31v, 289v, 320v, 116v): 22 minutes


def main(args):
    trait_type = 'quantitative'
    phenos = ('X50_irnt', 'X49_irnt',)
    covariates = 'X46_irnt'
    relatedness_cutoff = '0.125'
    num_markers = 2000
    n_threads = 6

    sparse_grm_root = f'{root}/tmp/sparse'
    sparse_grm_extension = f'_relatednessCutoff_{relatedness_cutoff}_{num_markers}_randomMarkersUsed.sparseGRM.mtx'

    # backend = LocalBackend(gsa_key_file='/Users/konradk/.hail/privateKeyData')
    backend = LocalBackend(gsa_key_file='/Users/konradk/.hail/ukb_exomes.json')
    p = Pipeline(default_image=SAIGE_DOCKER_IMAGE, default_storage='500Mi', backend=backend, default_cpu=n_threads)

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

    if args.gene_test:
        overwrite_sparse_grm = False
        if overwrite_sparse_grm:
            sparse_grm = create_sparse_grm(p, sparse_grm_root, ukb_for_grm_plink_path,
                                           relatedness_cutoff, num_markers, n_threads=n_threads)
        else:
            sparse_grm = p.read_input_group(**{ext: f'{sparse_grm_root}.{ext}' for ext in
                                               (sparse_grm_extension, f'{sparse_grm_extension}.sampleIDs.txt')})

        overwrite_null_models = True
        null_model_dir = f'{root}/null_glmm'
        if not overwrite_null_models:
            null_models_already_created = {x['path'] for x in hl.hadoop_ls(null_model_dir)}
        null_models = {}

        for pheno in phenos:
            null_glmm_root = f'{null_model_dir}/{pheno}'
            model_file_path = f'{null_glmm_root}.rda'
            gene_variance_ratio_file_path = f'{null_glmm_root}.gene.varianceRatio.txt'
            sparse_sigma_file_path = gene_variance_ratio_file_path + sparse_grm_extension.replace("GRM", "Sigma")

            if not overwrite_null_models and null_glmm_root in null_models_already_created:
                model_file = p.read_input(model_file_path)
                variance_ratio_file = p.read_input(gene_variance_ratio_file_path)
                sparse_sigma_file = p.read_input(sparse_sigma_file_path)
            else:
                fit_null_task = fit_null_glmm(p, null_glmm_root, pheno_path, pheno, trait_type, covariates,
                                              ukb_for_grm_plink_path,
                                              sparse_grm=sparse_grm, sparse_grm_extension=sparse_grm_extension,
                                              skip_model_fitting=not overwrite_null_models, n_threads=n_threads)
                model_file = fit_null_task.null_glmm.rda
                variance_ratio_file = fit_null_task.null_glmm['gene.varianceRatio.txt']
                sparse_sigma_file = fit_null_task.null_glmm[
                    f'gene.varianceRatio.txt{sparse_grm_extension.replace("GRM", "Sigma")}']
            null_models[pheno] = (model_file, variance_ratio_file, sparse_sigma_file)

        chrom_lengths = {'chr1': 248956422}  # hl.get_reference('GRCh38').lengths

        vcf_dir = f'{root}/vcf'
        overwrite_vcfs = True
        if not overwrite_vcfs:
            vcfs_already_created = {x['path'] for x in hl.hadoop_ls(vcf_dir)}
        chunk_size = int(1e6)
        vcfs = {}
        for chrom in range(1, 3):  # 23
            i = 0
            chromosome = f'chr{chrom}'
            for start_pos in range(1, chrom_lengths[chromosome], chunk_size):
                i += 1
                if i == 6:
                    break
                interval = f'{chromosome}:{start_pos}-{start_pos + chunk_size}'
                vcf_root = f'{vcf_dir}/test_{chromosome}_{start_pos}'
                if not overwrite_vcfs and vcf_root in vcfs_already_created:
                    vcf_file = p.read_input_group(**{'vcf.gz': f'{vcf_root}.vcf.gz',
                                                     'vcf.gz.tbi': f'{vcf_root}.vcf.gz.tbi'})
                    group_file = p.read_input(f'{vcf_root}.gene.txt')
                else:
                    vcf_task = extract_vcf_from_mt(p, vcf_root, get_ukb_gene_map_ht_path(), get_ukb_exomes_mt_path(),
                                                   sample_mapping_file, interval=interval)
                    vcf_file = vcf_task.vcf
                    group_file = vcf_task.group_file
                vcfs[interval] = (vcf_file, group_file)

        result_dir = f'{root}/result'
        for pheno in phenos:
            model_file, variance_ratio_file, sparse_sigma_file = null_models[pheno]
            for chrom in range(1, 3):  # 23
                i = 0
                chromosome = f'chr{chrom}'
                for start_pos in range(1, chrom_lengths[chromosome], chunk_size):
                    i += 1
                    if i == 6:
                        break
                    interval = f'{chromosome}:{start_pos}-{start_pos + chunk_size}'
                    vcf_file, group_file = vcfs[interval]
                    results_path = f'{result_dir}/{pheno}/result_{pheno}_{chromosome}_{start_pos}'
                    run_saige(p, results_path, model_file, variance_ratio_file, vcf_file, get_ukb_samples_file_path(),
                              group_file, sparse_sigma_file)

    p.run(dry_run=False, verbose=True, delete_scratch_on_exit=False)


if __name__ == '__main__':
    parser = argparse.ArgumentParser()

    parser.add_argument('--overwrite', help='Overwrite everything', action='store_true')
    parser.add_argument('--create_sparse_grm', help='Run single variant SAIGE', action='store_true')
    parser.add_argument('--variant_test', help='Run single variant SAIGE', action='store_true')
    parser.add_argument('--gene_test', help='Run SAIGE-GENE', action='store_true')
    parser.add_argument('--dry_run', help='Dry run only', action='store_true')
    parser.add_argument('--slack_channel', help='Send message to Slack channel/user', default='@konradjk')
    args = parser.parse_args()

    if args.slack_channel:
        try_slack(args.slack_channel, main, args)
    else:
        main(args)


