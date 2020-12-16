"""
Author: Lucy Garner

This script works with outputs from pipeline_pyscenic.py:
    1. Calculates regulon specificity scores
    2. Calculates AUCell z-score per cell
    3. Calculates AUCell z-score per annotation group (e.g. cell type)
"""

import os
import argparse
import pandas as pd
from pyscenic.rss import regulon_specificity_scores
import logging

os.environ['NUMEXPR_MAX_THREADS'] = '5'

parser = argparse.ArgumentParser()
parser.add_argument('--sample', default = 'merged-all',
                    help = 'sample name')
parser.add_argument('--exp_mtx', default = 'pyscenic_results.dir/normalised.dir/merged-all_filtered-expression.csv',    # genes x cells
                    help = 'path to h5ad-formatted hdf5 file containing AnnData object')
parser.add_argument('--aucell_output', default = 'pyscenic_results.dir/normalised.dir/merged-all_aucell.csv',
                    help = 'Output from pyscenic aucell, matrix of AUCell scores')
parser.add_argument('--annotation_input', default = 'stimulation-annotation.csv,cluster-annotation.csv',
                    help = 'files containing mapping between cell barcodes (1st column) and annotations of interest e.g. clusters (2nd column)')
parser.add_argument('-t', action = 'store_true',
                    help = 'whether to transpose the AUCell matrix - yes if started with csv file')
args = parser.parse_args()

logging.basicConfig(filename = 'rss_zscore.log', level = logging.DEBUG)
logging.info(args)

sample = args.sample
annotation_input = args.annotation_input.split(',')

annotation_files = []
annotation_names = []
for i, annotation in enumerate(annotation_input):
    annotation_name = annotation.split('.')[0]
    annotation = sample + '_' + annotation
    annotation_files.append(annotation)
    annotation_names.append(annotation_name)

if 'raw' in args.exp_mtx:
    datatype = 'raw'
elif 'normalised' in args.exp_mtx:
    datatype = 'normalised'

# Regulon specificity scores
logging.info("About to calculate regulon specificity scores")

# Make a dictionary of annotations that you want regulon specificity scores for
if args.t:
    exp_mtx = pd.read_csv(args.exp_mtx, index_col = 0).T
else:
    exp_mtx = pd.read_csv(args.exp_mtx, index_col = 0)

annotation_dict = {}
for i, annotation in enumerate(annotation_files):
    annotation = pd.read_csv('data.dir/' + annotation, squeeze = True, index_col = 0, header = 0)
    annotation = annotation[annotation.index.isin(exp_mtx.index)]
    annotation_dict[annotation_names[i]] = annotation

# Calculate RSS for each item (i.e. each annotation) in dictionary
if args.t:
    auc_mtx = pd.read_csv(args.aucell_output, index_col = 0).T
else:
    auc_mtx = pd.read_csv(args.aucell_output, index_col = 0)

for annotation_name in annotation_dict.keys():
    rss_louvain = regulon_specificity_scores(auc_mtx, annotation_dict[annotation_name])
    rss_louvain.insert(0, annotation_name, rss_louvain.index.values)
    rss_louvain.to_csv('pyscenic_results.dir/' + datatype + '.dir/' + sample + '.dir/' + annotation_name + '_RSS.csv', index = False)

logging.info("Calculated regulon specificity scores")

# Generate a Z-score for each regulon to enable comparison between regulons

# Score per cell
auc_mtx_Z = pd.DataFrame(index = auc_mtx.index)
for col in list(auc_mtx.columns):
    auc_mtx_Z[col] = (auc_mtx[col] - auc_mtx[col].mean()) / auc_mtx[col].std(ddof = 0)
auc_mtx_Z.to_csv('pyscenic_results.dir/' + datatype + '.dir/' + sample + '.dir/aucell_zscores.csv')

# Z-score per cluster/annotation group
# To find cell type specific regulators we use a Z score (i.e. the average AUCell
# score for the cells of a given type are standardized using the overall average AUCell scores
# and its standard deviation)
auc_mtx_long = auc_mtx.T
auc_mtx_long['regulon'] = auc_mtx_long.index.values
auc_mtx_long = pd.melt(auc_mtx_long, id_vars = 'regulon', var_name = 'cellid', value_name = 'score')

for annotation_name in annotation_dict.keys():
    annotation = annotation_dict[annotation_name]
    aucell_score = auc_mtx_long.merge(annotation, left_on = 'cellid', right_index = True).drop('cellid', axis = 1)
    aucell_score.columns = ['regulon', 'score', 'annotation']
    aucell_score = aucell_score.merge(aucell_score.groupby(['regulon', 'annotation'])['score'].mean().rename('annotation_mean').reset_index())
    aucell_score = aucell_score.merge(aucell_score.groupby(['regulon'])['score'].mean().rename('regulon_mean').reset_index())
    aucell_score = aucell_score.merge(aucell_score.groupby(['regulon'])['score'].std().rename('regulon_sd').reset_index())
    aucell_score['annotation_zscore'] = (aucell_score.annotation_mean - aucell_score.regulon_mean)/aucell_score.regulon_sd
    aucell_score = aucell_score.drop_duplicates(subset = ['regulon', 'annotation']).loc[:, ['regulon', 'annotation', 'annotation_zscore']]
    aucell_score.to_csv('pyscenic_results.dir/' + datatype + '.dir/' + sample + '.dir/aucell_zscores_' + annotation_name + '.csv', index = False)

    # Cell specific regulators have a Z-score greater than 3
    aucell_score_3 = aucell_score[(aucell_score.annotation_zscore >= 3.0)].sort_values('annotation_zscore', ascending = False)
    aucell_score_3.to_csv('pyscenic_results.dir/' + datatype + '.dir/' + sample + '.dir/aucell_zscores_' + annotation_name + '_specific_regulators_3.csv', index = False)

    # Cell specific regulators have a Z-score greater than 2
    aucell_score_2 = aucell_score[(aucell_score.annotation_zscore >= 2.0)].sort_values('annotation_zscore', ascending = False)
    aucell_score_2.to_csv('pyscenic_results.dir/' + datatype + '.dir/' + sample + '.dir/aucell_zscores_' + annotation_name + '_specific_regulators_2.csv', index = False)
