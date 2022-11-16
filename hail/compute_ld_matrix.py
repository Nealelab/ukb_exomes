#!/usr/bin/env python3

__author__ = 'konradk'

import sys
print(sys.path)
import subprocess
import hailtop.batch as hb
# subprocess.check_output(['unzip', sys.path[1]])
# subprocess.check_output(['ls'])
import argparse
from gnomad.utils.vep import process_consequences
from gnomad.utils import slack
from ukbb_qc.resources.variant_qc import var_annotations_ht_path
from ukbb_common import *
from ukb_exomes import *

ukb_exomes_path = 'gs://ukbb-pharma-exome-analysis/500k_temp/500k.vds'


def ld_matrix_path(gene_id, tranche: str = CURRENT_TRANCHE):
    return f'{public_bucket}/{tranche}/ld/{gene_id}/ld.bm'


def ld_index_path(gene_id, tranche: str = CURRENT_TRANCHE):
    return f'{public_bucket}/{tranche}/ld/{gene_id}/ld.ht'


def get_cohort_index(call_stats_ht):
    freq_index = hl.zip_with_index(call_stats_ht.freq_meta).find(lambda f: (hl.len(f[1]) == 1) &
                                                                           (f[1].get('group') == 'cohort'))[0]
    return freq_index.collect()[0]


async def write_ld_matrix(chrom, overwrite=False, gene_id=None, interval=None, radius=1000000):
    import hail as hl
    # hl.init(default_reference='GRCh38', master='local[8]')
    await hl.init_batch(billing_project='ukb_pharma', remote_tmpdir='gs://ukb-pharma-exome-analysis-temp/ld_temp',
                        default_reference='GRCh38')
    do_write_ld_matrix(chrom, overwrite, interval, gene_id, radius)


def do_write_ld_matrix(chrom, overwrite=False, gene_id=None, interval=None, radius=1000000):
    def read_vds(chrom):
        vds = hl.vds.read_vds(ukb_exomes_path)
        return hl.vds.filter_chromosomes(vds, keep=chrom)

    vds = read_vds(chrom)
    if gene_id is not None and interval is not None:
        vds = hl.vds.filter_intervals(vds, [hl.eval(hl.parse_locus_interval(interval))])
    mt = hl.vds.to_dense_mt(vds)
    out_path = gene_id if gene_id is not None else chrom
    mt.rows().write(ld_index_path(out_path), overwrite)
    # ld = hl.ld_matrix(mt.GT.n_alt_alleles(), mt.locus, radius, block_size=1024)
    # ld = ld.sparsify_triangle()
    # ld.write(ld_matrix_path(out_path), overwrite=overwrite)


def run_async(f, *args, **kwargs):
    import asyncio
    asyncio.run(f(*args, **kwargs))


def main(args):
    hl.init(default_reference='GRCh38', log='/pre_process.log')

    if args.pre_process_data:
        vds = hl.vds.read_vds("gs://gnomad/v4.0/raw/exomes/gnomad_v4.0.vds")
        meta_ht = hl.read_table("gs://broad-ukbb/broad.freeze_7/sample_qc/meta.ht")
        meta_ht = meta_ht.filter(meta_ht.sample_filters.high_quality)
        vds = hl.vds.split_multi(hl.vds.filter_samples(vds, meta_ht, remove_dead_alleles=True))
        call_stats_ht = hl.read_table(ukb.var_annotations_ht_path('ukb_freq', *TRANCHE_DATA[CURRENT_TRANCHE]))
        freq_index = get_cohort_index(call_stats_ht)
        var = vds.variant_data.annotate_rows(call_stats=call_stats_ht[vds.variant_data.row_key].freq[freq_index])
        var = var.filter_rows(var.call_stats.AC >= args.min_ac)
        vds = hl.vds.variant_dataset.VariantDataset(vds.reference_data, var)
        vds.write(ukb_exomes_path, overwrite=args.overwrite)

    if args.run_ld_dataproc:
        chromosomes = [f'chr{x}' for x in range(1, 23)] + ['chrX', 'chrY']
        chromosomes = ['chrY', 'chrX', 'chr22', 'chr21', 'chr20']
        # chromosomes = ['chr22']
        for chrom in chromosomes[::-1]:
            if chrom != 'chr22':
                do_write_ld_matrix(chrom, True)

    if args.run_ld:
        backend = hb.ServiceBackend(billing_project='ukb_pharma', remote_tmpdir='gs://ukb-pharma-exome-analysis-temp/ld_temp')
        b = hb.Batch(default_python_image='hailgenetics/hail:0.2.93', backend=backend)

        if args.by_gene:
            ht = hl.read_matrix_table(get_results_mt_path('gene')).rows()
            intervals = ht.group_by(ht.gene_id, ht.gene_symbol).aggregate(
                interval=hl.agg.take(ht.interval, 1)[0]
            ).collect()
            for gene in intervals:
                j = b.new_python_job(name=gene.gene_id, attributes={'gene_id': gene.gene_id, 'gene_symbol': gene.gene_symbol, 'interval': str(gene.interval)})
                j.call(run_async, write_ld_matrix, gene.interval.start.contig, args.overwrite, gene.gene_id, str(gene.interval), )
                break
        else:
            # chromosomes = [f'chr{x}' for x in range(1, 23)] + ['chrX', 'chrY']
            chromosomes = ['chr1']
            for chrom in chromosomes:
                j = b.new_python_job(name=chrom)
                j.call(run_async, write_ld_matrix, chrom, args.overwrite)
        b.run()

    # confirming
    ht = hl.read_table(ld_index_path('chr22'))
    rel_ht = hl.read_matrix_table(get_results_mt_path('variant')).rows()
    rel_ht = rel_ht.filter(rel_ht.locus.contig == 'chr22')
    rel_ht.filter(hl.is_missing(ht[rel_ht.key])).show()
    # Examples
    # +----------------+-----------------+----------------------------+----------+--------------+---------------+---------------+---------------+
    # | locus          | alleles         | markerID                   | gene     | annotation   | call_stats.AC | call_stats.AF | call_stats.AN |
    # +----------------+-----------------+----------------------------+----------+--------------+---------------+---------------+---------------+
    # | locus<GRCh38>  | array<str>      | str                        | str      | str          |         int32 |       float64 |         int32 |
    # +----------------+-----------------+----------------------------+----------+--------------+---------------+---------------+---------------+
    # | chr22:15528165 | ["CCCTTGA","C"] | "chr22:15528165_CCCTTGA/C" | "OR11H1" | "missense"   |             1 |      1.36e-06 |        736778 |
    # | chr22:15528172 | ["C","T"]       | "chr22:15528172_C/T"       | "OR11H1" | "missense"   |             1 |      1.31e-06 |        761422 |
    # | chr22:15528182 | ["C","T"]       | "chr22:15528182_C/T"       | "OR11H1" | "synonymous" |             5 |      6.60e-06 |        757018 |
    # | chr22:15528188 | ["C","T"]       | "chr22:15528188_C/T"       | "OR11H1" | "synonymous" |           901 |      1.16e-03 |        777558 |
    # | chr22:15528190 | ["T","TGC"]     | "chr22:15528190_T/TGC"     | "OR11H1" | "pLoF"       |             1 |      1.26e-06 |        795620 |
    # | chr22:15528200 | ["C","T"]       | "chr22:15528200_C/T"       | "OR11H1" | "synonymous" |             3 |      3.65e-06 |        821090 |
    # | chr22:15528218 | ["CT","C"]      | "chr22:15528218_CT/C"      | "OR11H1" | "pLoF"       |             1 |      1.21e-06 |        826654 |
    # | chr22:15528223 | ["C","T"]       | "chr22:15528223_C/T"       | "OR11H1" | "missense"   |            55 |      6.65e-05 |        826834 |
    # | chr22:15528250 | ["G","A"]       | "chr22:15528250_G/A"       | "OR11H1" | "missense"   |             2 |      2.42e-06 |        827370 |
    # | chr22:15528264 | ["T","A"]       | "chr22:15528264_T/A"       | "OR11H1" | "missense"   |             1 |      1.21e-06 |        827370 |
    # +----------------+-----------------+----------------------------+----------+--------------+---------------+---------------+---------------+
    # TODO: figure out why these are still missing
    # TODO: Got here
    vds = hl.vds.read_vds("gs://gnomad/v4.0/raw/exomes/gnomad_v4.0.vds")
    meta_ht = hl.read_table("gs://broad-ukbb/broad.freeze_7/sample_qc/meta.ht")
    meta_ht = meta_ht.filter(meta_ht.sample_filters.high_quality)
    vds = hl.vds.filter_samples(vds, meta_ht, remove_dead_alleles=True)
    call_stats_ht = hl.read_table(ukb.var_annotations_ht_path('ukb_freq', *TRANCHE_DATA[CURRENT_TRANCHE]))
    freq_index = get_cohort_index(call_stats_ht)
    var = vds.variant_data.annotate_rows(call_stats=call_stats_ht[vds.variant_data.row_key].freq[freq_index])
    var = var.filter_rows(var.call_stats.AC >= args.min_ac)
    vds = hl.vds.variant_dataset.VariantDataset(vds.reference_data, var)
    mt = vds.variant_data.filter_rows(vds.variant_data.locus == hl.parse_locus("chr22:15528188"))
    mt.rows().show()

    # temp stuff for figuring out AN differences
    # vds = hl.vds.read_vds("gs://gnomad/v4.0/raw/exomes/gnomad_v4.0.vds")
    hl.read_table("gs://broad-ukbb/broad.freeze_6/sample_qc/meta.ht").naive_coalesce(100).write('/tmp/300k.ht')
    hl.read_table("gs://broad-ukbb/broad.freeze_7/sample_qc/meta.ht").naive_coalesce(100).write('/tmp/450k.ht')
    meta_300k_ht = hl.read_table('/tmp/300k.ht')
    meta_450k_ht = hl.read_table('/tmp/450k.ht')
    meta_300k_ht.aggregate(hl.agg.count_where(meta_300k_ht.sample_filters.high_quality))
    meta_450k_ht.aggregate(hl.agg.count_where(meta_450k_ht.sample_filters.high_quality))
    meta_300k_ht.aggregate(hl.agg.count_where(meta_300k_ht.sample_filters.high_quality
                                              & ~meta_300k_ht.sample_filters.duplicate &
                                              (meta_300k_ht.gnomad_pc_project_pop_data.pop == 'nfe')))
    meta_450k_ht.aggregate(hl.agg.count_where(meta_450k_ht.sample_filters.high_quality
                                              & ~meta_450k_ht.sample_filters.duplicate &
                                              (meta_450k_ht.pan_ancestry_meta.pop == 'EUR')))


if __name__ == '__main__':
    parser = argparse.ArgumentParser()

    parser.add_argument('--overwrite', help='Overwrite everything', action='store_true')
    parser.add_argument('--pre_process_data', help='Overwrite everything', action='store_true')
    parser.add_argument('--run_ld', help='Overwrite everything', action='store_true')
    parser.add_argument('--run_ld_dataproc', help='Overwrite everything', action='store_true')
    parser.add_argument('--by_gene', help='Overwrite everything', action='store_true')
    parser.add_argument('--radius', help='Radius at which to calculate LD information (bp; default 1e6)', default=1000000, type=int)
    parser.add_argument('--min_ac', help='Minimum allele count to compute LD scores (default 1)', default=1, type=float)
    parser.add_argument('--slack_channel', help='Send message to Slack channel/user')
    args = parser.parse_args()

    if args.slack_channel:
        from slack_token_pkg.slack_creds import slack_token
        with slack.slack_notifications(slack_token, args.slack_channel):
            main(args)
    else:
        main(args)
