"""
Author: Lucy Garner

This script works with outputs from pipeline_pyscenic.py:
    1. Checks whether the default threshold for pyscenic_aucell is appropriate
    2. Generates a heatmap of AUCell scores
"""

import os
import argparse
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import logging

os.environ['NUMEXPR_MAX_THREADS'] = '15'

parser = argparse.ArgumentParser()
parser.add_argument('--sample', default = 'merged-all',
                    help = 'sample name')
parser.add_argument('--exp_mtx', default = 'pyscenic_results.dir/normalised.dir/merged-all_filtered-expression.csv',    # genes x cells
                    help = 'path to csv file containing filtered expression values used for pySCENIC')
parser.add_argument('--aucell_output', default = 'pyscenic_results.dir/normalised.dir/merged-all_aucell.csv',
                    help = 'output from pyscenic aucell, matrix of AUCell scores')
parser.add_argument('-t', action = 'store_true',
                    help = "whether to transpose the AUCell matrix - yes if started with csv file")
args = parser.parse_args()

logging.basicConfig(filename = 'aucell_heatmap.log', level = logging.DEBUG)
logging.info(args)

sample = args.sample

if "raw" in args.exp_mtx:
    datatype = "raw"
elif "normalised" in args.exp_mtx:
    datatype = "normalised"

# Check how many genes are detected per cell
# Decide whether default AUCell threshold of 0.05 is appropriate
if args.t:
    exp_mtx = pd.read_csv(args.exp_mtx, index_col = 0).T
else:
    exp_mtx = pd.read_csv(args.exp_mtx, index_col = 0)

nGenesDetectedPerCell = np.sum(exp_mtx > 0, axis = 1)
percentiles = nGenesDetectedPerCell.quantile([.01, .05, .10, .50, 1])
logging.info(percentiles)

fig, ax = plt.subplots(1, 1, figsize = (8, 5), dpi = 150)
sns.distplot(nGenesDetectedPerCell, norm_hist = False, kde = False, bins = 'fd')
for i,x in enumerate(percentiles):
    fig.gca().axvline(x = x, ymin = 0,ymax = 1, color = 'red')
    ax.text(x = x, y = ax.get_ylim()[1], s = f'{int(x)} ({percentiles.index.values[i]*100}%)',
            color= 'red', rotation = 30, size = 'x-small', rotation_mode = 'anchor')
ax.set_xlabel('# of genes')
ax.set_ylabel('# of cells')
fig.tight_layout()
plt.savefig("plots.dir/" + datatype + ".dir/" + sample + "_number_of_genes_expressed_above_AUC_thresholds.png")
plt.close()

# Heatmap of AUCell scores
if args.t:
    auc_mtx = pd.read_csv(args.aucell_output, index_col = 0).T
else:
    auc_mtx = pd.read_csv(args.aucell_output, index_col = 0)

logging.info("Read in AUCell matrix")
sns_plot = sns.clustermap(auc_mtx, figsize = (12, 12))
sns_plot.savefig("plots.dir/" + datatype + ".dir/" + sample + "_aucell_heatmap.png")

logging.info("Plotted heatmap of AUCell scores")
