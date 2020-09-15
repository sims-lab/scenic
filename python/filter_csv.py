# Script to filter csv file for use in pySCENIC

import argparse
import pandas as pd

parser = argparse.ArgumentParser()
parser.add_argument('--input',
                    default = 'data.dir/*_raw-expression.csv',
                    help = 'path to h5ad-formatted hdf5 file containing AnnData object')
parser.add_argument('--umi_counts', default = 3, 
                    help = 'UMI counts to use for gene filtering, genes must have at least umi_counts * min_percent * numberofcells counts to be retained')
parser.add_argument('--min_percent', default = 0.01, 
                    help = 'percentage of cells for gene filtering, genes must have at least min_percent * numberofcells counts to be retained')
parser.add_argument('--output', default = "filtered-raw-expression.csv", 
                    help = "path to output file")
args = parser.parse_args()

exp_mtx = pd.read_csv(args.input, index_col = 0)
ncells = exp_mtx.shape[1]

# pySCENIC thresholds
min_counts_per_gene = int(args.umi_counts) * float(args.min_percent) * ncells    # minimum total counts per gene
print("min counts per gene: ", min_counts_per_gene)

min_samples = float(args.min_percent) * ncells    # minimum number of cells in which gene is expressed
print("min_samples: ", min_samples)

# Rows to filter
rows_to_filter = exp_mtx.sum(axis = 1) > min_counts_per_gene
exp_mtx = exp_mtx[rows_to_filter]

rows_to_filter = exp_mtx.astype(bool).sum(axis = 1) > min_samples
exp_mtx = exp_mtx[rows_to_filter]

exp_mtx.to_csv(args.output)