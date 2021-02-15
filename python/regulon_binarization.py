"""
Author: Lucy Garner

This script works with outputs from pipeline_pyscenic.py:
    1. Binarizes the AUCell matrix
    2. Plots thresholds used for binarization
"""

import os
import argparse
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import logging
from pyscenic.binarization import binarize

os.environ['NUMEXPR_MAX_THREADS'] = '25'

parser = argparse.ArgumentParser()
parser.add_argument('--sample', default = 'merged-all',
                    help = 'sample name')
parser.add_argument('--aucell_output', default = 'pyscenic_results.dir/normalised.dir/all-samples.dir/aucell.csv',
                    help = 'Output from pyscenic aucell, matrix of AUCell scores')
parser.add_argument('-t', action = 'store_true',
                    help = "whether to transpose the AUCell matrix - yes if started with csv file")
parser.add_argument('--binarize_threads', default = 25,
                    help = "number of threads to use for generating binary AUCell matrix")
parser.add_argument('--custom_auc_thresholds',
                    help = "custom thresholds for generating binarizing the AUCell matrix")
args = parser.parse_args()

logging.basicConfig(filename = 'regulon_binarization.log', level = logging.DEBUG)
logging.info(args)

sample = args.sample

if "raw" in args.aucell_output:
    datatype = "raw"
elif "normalised" in args.aucell_output:
    datatype = "normalised"

# Functions

# Plot "binarization" process for the given regulon
def plot_binarization(auc_mtx: pd.DataFrame, regulon_name: str, threshold: float, bins: int = 200, ax = None) -> None:
    """
    Plot the "binarization" process for the given regulon.
    :param auc_mtx: The dataframe with the AUC values for all cells and regulons (n_cells x n_regulons).
    :param regulon_name: The name of the regulon.
    :param bins: The number of bins to use in the AUC histogram.
    :param threshold: The threshold to use for binarization.
    """
    if ax is None:
        ax = plt.gca()

    sns.distplot(auc_mtx[regulon_name], ax = ax, norm_hist = True, bins = bins)

    ylim = ax.get_ylim()
    ax.plot([threshold] * 2, ylim, 'r:')
    ax.set_ylim(ylim)
    ax.set_xlabel('AUC')
    ax.set_ylabel('#')
    ax.set_title(regulon_name)

    plt.savefig("plots.dir/" + datatype + ".dir/" + sample + ".dir/" + regulon_name + "_auc_threshold.png")
    plt.close()

# Analysis

# Import AUCell scores
if args.t:
    auc_mtx = pd.read_csv(args.aucell_output, index_col = 0).T
else:
    auc_mtx = pd.read_csv(args.aucell_output, index_col = 0)

# "Binarize" the supplied AUCell matrix, i.e. decide if for each cell in the matrix
# a regulon is active or not based on the bimodal distribution of the AUC values for that regulon
logging.info("About to perform binarization")

if args.custom_auc_thresholds is None:
    logging.info("No custom AUCell thresholds provided")
    logging.info(auc_mtx.head())
    binary_mtx, auc_thresholds = binarize(auc_mtx, num_workers = int(args.binarize_threads))    # binarized dataframe and AUC threshold used for each regulon
else:
    logging.info("Custom AUCell thresholds provided")
    logging.info(auc_mtx.head())
    custom_auc_thresholds = pd.read_csv(args.custom_auc_thresholds, names = ['regulon', 'threshold'])
    custom_auc_thresholds = dict(zip(custom_auc_thresholds.regulon, custom_auc_thresholds.threshold))
    binary_mtx, auc_thresholds = binarize(auc_mtx, threshold_overides = custom_auc_thresholds, num_workers = int(args.binarize_threads))    # binarized dataframe and AUC threshold used for each regulon

logging.info(binary_mtx.head())

if args.custom_auc_thresholds is None:
    binary_mtx.to_csv("pyscenic_results.dir/" + datatype + ".dir/" + sample + ".dir/binary_matrix.csv")
    auc_thresholds.to_csv("pyscenic_results.dir/" + datatype + ".dir/" + sample + ".dir/aucell_thresholds.csv", header = False, index = True)
else:
    binary_mtx.to_csv("pyscenic_results.dir/" + datatype + ".dir/" + sample + ".dir/binary_matrix_custom_thresholds.csv")
    auc_thresholds.to_csv("pyscenic_results.dir/" + datatype + ".dir/" + sample + ".dir/aucell_thresholds.csv", header = False, index = True)

logging.info("Finished binarization")

# Loop over all regulons and plot the thresholds
for regulon in binary_mtx.columns:
    plot_binarization(auc_mtx, regulon_name = regulon, threshold = auc_thresholds[regulon])
logging.info("Plotted thresholds for binarization")
