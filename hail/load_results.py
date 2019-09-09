#!/usr/bin/env python3

__author__ = 'konradk'

from gnomad_hail import *
from ukb_exomes.resources import *

temp_bucket = 'gs://ukbb-pharma-exome-analysis-temp'
results_dir = 'gs://ukbb-pharma-exome-analysis/data/result/'
final_gene_results_ht = 'gs://ukbb-pharma-exome-analysis/data/results.ht'
final_variant_results_ht = 'gs://ukbb-pharma-exome-analysis/data/variant_results.ht'


def union_mts_by_tree(all_mts, temp_dir):
    chunk_size = int(len(all_mts) ** 0.5) + 1
    outer_mts = []
    for i in range(chunk_size):
        if i * chunk_size >= len(all_mts): break
        mt = all_mts[i * chunk_size]
        for j in range(1, chunk_size):
            if i * chunk_size + j >= len(all_mts): break
            mt = mt.union_cols(all_mts[i * chunk_size + j])
        outer_mts.append(mt.checkpoint(f'{temp_dir}/temp_output_{i}.mt', overwrite=True))
    mt = outer_mts[0]
    for next_mt in outer_mts[1:]:
        mt = mt.union_cols(next_mt)
    return mt


def get_pheno_annotations():
    cont_ht = hl.read_matrix_table(get_ukb_pheno_mt_path('continuous', 'full')).cols()
    cont_ht = cont_ht.key_by(pheno=hl.str(cont_ht.pheno), coding=cont_ht.coding)
    cont_ht = cont_ht.select(meaning=hl.coalesce(cont_ht.both_sexes_pheno.Field,
                                                 cont_ht.females_pheno.Field,
                                                 cont_ht.males_pheno.Field),
                             path=hl.coalesce(cont_ht.both_sexes_pheno.Path,
                                              cont_ht.females_pheno.Path,
                                              cont_ht.males_pheno.Path)
                             )
    icd_ht = hl.read_table(icd_full_codings_ht_path)
    icd_ht = icd_ht.key_by(pheno=icd_ht.coding, coding='')
    icd_ht = icd_ht.select(meaning=icd_ht.short_meaning, path=icd_ht.meaning)
    pheno_ht = cont_ht.union(icd_ht)
    return pheno_ht


def load_variant_data(directory_root, pheno, overwrite: bool = False):
    if len(pheno) == 1:
        pheno.append('')
    pheno, coding = pheno
    output_ht_path = f'{directory_root}/variant_results.ht'
    ht = hl.import_table(f'{directory_root}/*.single.txt', delimiter=' ', impute=True)
    locus_alleles = ht.markerID.split('_')
    ht = ht.key_by(locus=hl.parse_locus(locus_alleles[0]), alleles=locus_alleles[1].split('/'),
                   pheno=pheno, coding=coding).distinct().naive_coalesce(50)
    ht = ht.transmute(Pvalue=ht['p.value'])
    ht = ht.annotate(**get_vep_formatted_data()[
        hl.struct(locus=ht.locus, alleles=ht.alleles)])  # TODO: fix this for variants that overlap multiple genes
    ht = ht.checkpoint(output_ht_path, overwrite=overwrite, _read_if_exists=not overwrite)
    mt = ht.to_matrix_table(['locus', 'alleles'], ['pheno', 'coding'],
                            ['markerID', 'gene', 'annotation'], [])
    mt.checkpoint(output_ht_path.replace('.ht', '.mt'), overwrite=overwrite, _read_if_exists=not overwrite)


def main(args):
    hl.init(default_reference='GRCh38')
    all_phenos = hl.hadoop_ls(results_dir)

    if args.load_gene_results:
        all_gene_outputs = []
        all_gene_mt_outputs = []
        for directory in all_phenos:
            # if 'E109' in directory['path'] or \
            #         'E119' in directory['path'] or \
            #         'K509' in directory['path'] or \
            #         'J459' in directory['path']: continue
            pheno = directory['path'].split('/')[-1].split('-')
            load_gene_data(directory['path'], pheno)
            all_gene_outputs.append(f'{directory["path"]}/gene_results.ht')
            all_gene_mt_outputs.append(f'{directory["path"]}/gene_results.mt')

        pheno_ht = get_pheno_annotations()
        all_hts = list(map(lambda x: hl.read_table(x), all_gene_outputs))
        ht = all_hts[0].union(*all_hts[1:], unify=True)
        ht = ht.annotate(**pheno_ht[hl.struct(pheno=ht.pheno, coding=ht.coding)])
        ht = ht.checkpoint(final_gene_results_ht, overwrite=args.overwrite, _read_if_exists=not args.overwrite)

        all_mts = list(map(lambda x: hl.read_matrix_table(x), all_gene_mt_outputs))
        mt = union_mts_by_tree(all_mts, temp_bucket)
        mt = mt.annotate_cols(**pheno_ht[hl.struct(pheno=mt.pheno, coding=mt.coding)])
        mt = mt.checkpoint(final_gene_results_ht.replace('.ht', '.mt'), overwrite=args.overwrite, _read_if_exists=not args.overwrite)

    if args.load_variant_results:
        all_variant_outputs = []
        all_variant_mt_outputs = []
        for directory in all_phenos:
            # if 'K519' not in directory['path']: continue
            pheno = directory['path'].split('/')[-1].split('-')
            load_variant_data(directory["path"], pheno)
            all_variant_outputs.append(f'{directory["path"]}/variant_results.ht')
            all_variant_mt_outputs.append(f'{directory["path"]}/variant_results.mt')

        pheno_ht = get_pheno_annotations()
        all_hts = list(map(lambda x: hl.read_table(x), all_variant_outputs))
        ht = all_hts[0].union(*all_hts[1:])
        ht = ht.annotate(**pheno_ht[hl.struct(pheno=ht.pheno, coding=ht.coding)])
        ht = ht.checkpoint(final_variant_results_ht, overwrite=args.overwrite, _read_if_exists=not args.overwrite)

        # all_mts = list(map(lambda x: hl.read_matrix_table(x), all_variant_mt_outputs))
        # mt = union_mts_by_tree(all_mts, temp_bucket)
        # mt = mt.annotate_cols(**pheno_ht[hl.struct(pheno=mt.pheno, coding=mt.coding)])
        # mt = mt.checkpoint(final_variant_results_ht.replace('.ht', '.mt'), overwrite=args.overwrite, _read_if_exists=not args.overwrite)


def load_gene_data(directory, pheno, overwrite: bool = False):
    if len(pheno) == 1:
        pheno.append('')
    pheno, coding = pheno
    output_ht_path = f'{directory}/gene_results.ht'
    print(f'Loading: {directory}/*.gene.txt ...')
    ht = hl.import_table(f'{directory}/*.gene.txt', delimiter=' ', impute=True,
                         types={'Pvalue_SKAT': hl.tfloat64})
    if 'Pvalue_skato_NA' not in list(ht.row):
        ht = ht.annotate(Pvalue_skato_NA=hl.null(hl.tfloat64),
                         Pvalue_burden_NA=hl.null(hl.tfloat64),
                         Pvalue_skat_NA=hl.null(hl.tfloat64))
    fields = ht.Gene.split('_')
    gene_ht = hl.read_table(get_ukb_gene_map_ht_path(False)).select('interval').distinct()
    ht = ht.key_by(gene_id=fields[0], gene_symbol=fields[1], annotation=fields[2],
                   pheno=pheno, coding=coding).drop('Gene').naive_coalesce(10)
    ht = ht.annotate(total_variants=hl.sum([v for k, v in list(ht.row_value.items()) if 'Nmarker' in k]),
                     interval=gene_ht.key_by('gene_id')[ht.gene_id].interval)
    ht = ht.checkpoint(output_ht_path, overwrite=overwrite, _read_if_exists=not overwrite)
    mt = ht.to_matrix_table(['gene_symbol', 'gene_id', 'annotation', 'interval'],
                            ['pheno', 'coding'], [], [])
    mt.checkpoint(output_ht_path.replace('.ht', '.mt'), overwrite=overwrite, _read_if_exists=not overwrite)


def get_vep_formatted_data():
    ht = hl.read_table(get_ukb_vep_path())
    ht = process_consequences(ht)
    ht = ht.explode(ht.vep.worst_csq_by_gene_canonical)
    ht = ht.annotate(
        gene=ht.vep.worst_csq_by_gene_canonical.gene_symbol,
        annotation=hl.case(missing_false=True)
            .when(ht.vep.worst_csq_by_gene_canonical.lof == 'HC', 'pLoF')
            .when(ht.vep.worst_csq_by_gene_canonical.lof == 'LC', 'LC')
            .when((ht.vep.worst_csq_by_gene_canonical.most_severe_consequence == 'missense_variant') |
                  (ht.vep.worst_csq_by_gene_canonical.most_severe_consequence == 'inframe_insertion') |
                  (ht.vep.worst_csq_by_gene_canonical.most_severe_consequence == 'inframe_deletion'),
                  'missense')
            .when(ht.vep.worst_csq_by_gene_canonical.most_severe_consequence == 'synonymous_variant',
                  'synonymous')
            .or_missing())
    return ht.select('gene', 'annotation')


if __name__ == '__main__':
    parser = argparse.ArgumentParser()

    parser.add_argument('--overwrite', help='Overwrite everything', action='store_true')
    parser.add_argument('--load_gene_results', help='Overwrite everything', action='store_true')
    parser.add_argument('--load_variant_results', help='Overwrite everything', action='store_true')
    parser.add_argument('--slack_channel', help='Send message to Slack channel/user', default='@konradjk')
    args = parser.parse_args()

    if args.slack_channel:
        try_slack(args.slack_channel, main, args)
    else:
        main(args)