# Script to extract raw expression matrix from h5ad-formatted hdf5 file containing 
# AnnData object

# Script also filters genes to remove those with very low expression

import scanpy as sc
import pandas as pd
import numpy as np
import argparse

parser = argparse.ArgumentParser()
parser.add_argument('--input',
                    default = '../../data/anndata.h5ad',
                    help = 'path to h5ad-formatted hdf5 file containing AnnData object')
parser.add_argument('--umi_counts', default = 3, 
                    help = 'UMI counts to use for gene filtering, genes must have at least umi_counts * min_percent * numberofcells counts to be retained')
parser.add_argument('--min_percent', default = 0.01, 
                    help = 'percentage of cells for gene filtering, genes must have at least min_percent * numberofcells counts to be retained')
parser.add_argument('--downsample', action = 'store_true', 
                    help = 'include if you want to downsample your expression matrix, False by default')
parser.add_argument('--cell_number', default = 10000, 
                    help = "cell number to downsample to, default is 10000, must be used in combination with --downsample")
parser.add_argument('--output', default = "filtered-raw-expression.csv", 
                    help = "path to output file")
args = parser.parse_args()

# Read in h5ad-formatted hdf5 file
adata = sc.read_h5ad(args.input)

# pySCENIC thresholds

ncells = adata.shape[0]

min_counts_per_gene = int(args.umi_counts) * float(args.min_percent) * ncells    # minimum total counts per gene
print("min counts per gene: ", min_counts_per_gene)

min_samples = 0.01 * float(args.min_percent)    # minimum number of cells in which gene is expressed
print("min_samples: ", min_samples)

# Genes to filter
genes_to_filter = adata.layers["counts"].sum(axis = 0) > min_counts_per_gene
genes_to_filter = genes_to_filter.T.tolist()

adata = adata[:, np.where(genes_to_filter)[0]]

# Further filtering of genes
total_nonzero = []
for row in adata.layers["counts"].T:
    nonzero = row.count_nonzero()
    total_nonzero.append(nonzero)

genes_to_filter = [item > min_samples for item in total_nonzero]

adata = adata[:, np.where(genes_to_filter)[0]]

# Filter AnnData object to keep random subset of cells

if args.downsample:
    # Downsample by taking a random sample of cells from raw_counts
    keep = np.random.randint(adata.layers["counts"].shape[0], size = int(args.cell_number))
    adata = adata[keep, :]

# Write out dense matrix

raw_counts_dense = pd.DataFrame(adata.layers["counts"].todense())

gene_ids = adata.var_names
cell_names = adata.obs.index

raw_counts_dense.index = cell_names
raw_counts_dense.columns = gene_ids

raw_counts_dense.to_csv(args.output)
