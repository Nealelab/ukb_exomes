#!/usr/bin/env python3

__author__ = 'konradk'

from ukb_exomes import *
import tempfile


def main(args):
    hl.init(master=f'local[{args.n_threads}]',
            log=hl.utils.timestamp_path(os.path.join(tempfile.gettempdir(), 'load_results'), suffix='.log'),
            default_reference='GRCh38')

    cases, controls = get_cases_and_controls_from_log(f'{args.input_dir}/result_{args.pheno}')

    load_gene_data(args.input_dir, args.pheno, args.coding, args.gene_map_ht_raw_path, cases, controls, args.overwrite)
    load_variant_data(args.input_dir, args.pheno, args.coding, args.ukb_vep_ht_path, cases, controls, args.overwrite)


if __name__ == '__main__':
    parser = argparse.ArgumentParser()

    parser.add_argument('--input_dir', help='Input directory', required=True)
    parser.add_argument('--pheno', help='Phenotype ID', required=True)
    parser.add_argument('--coding', help='Phenotype coding', default='')
    parser.add_argument('--gene_map_ht_raw_path', help='Path to raw gene map', required=True)
    parser.add_argument('--ukb_vep_ht_path', help='Path to UKB VEP data', required=True)
    parser.add_argument('--n_threads', help='Number of threads to run', type=int, default=8)
    parser.add_argument('--overwrite', help='Overwrite everything', action='store_true')
    args = parser.parse_args()

    main(args)