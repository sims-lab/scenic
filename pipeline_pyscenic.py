"""
=================
Pipeline pySCENIC
=================

Authors: Lucy Garner and Devika Agarwal

Overview
========
This pipeline performs pySCENIC analysis (https://pyscenic.readthedocs.io/en/latest/index.html)
including:
1. Deriving co-expression modules
2. Finding enriched motifs and corresponding target genes for modules
3. Quantifying activity of gene signatures/regulons across single cells

Usage
=====

See https://cgat-core.readthedocs.io/en/latest/ for how to use CGAT-core pipelines

Configuration
-------------
The pipeline requires a pipeline.yml configuration file

Default configuration files can be generated by executing:
    python pipeline_pyscenic.py config

Input files
-----------
Four files as follows:

1. Raw or normalised expression matrix:
NB: Only currently works for CSV

Can be supplied as h5ad-formatted hdf5 file containing an AnnData object or as a csv (genes as rows, cells as columns) file
Should have been filtered to remove low quality cells and doublets
AnnData object can be downsampled within pipeline, whereas csv files should be downsampled before input

h5ad file should be named *_anndata.h5ad and csv file should be named *_raw-expression.csv/*_normalised-expression.csv
Files should be inside a data.dir directory

2. A list of TFs from the genome under study
(inside resources.dir directory; tfs_list in pipeline.yml)

3. Databases ranking the whole genome of your species of interest based on regulatory features (i.e. transcription factors)
Ranking databases are typically stored in the feather format and can be downloaded from cisTargetDBs
(inside resources.dir directory; database_fname in pipeline.yml)

4. Motif annotation database providing the missing link between an enriched motif and
the transcription factor that binds this motif
This pipeline needs a TSV text file where every line represents a particular annotation
(inside resources.dir directory; annotations in pipeline.yml)

Dependencies
------------
This pipeline requires: cgatcore, pyscenic and dependencies

Pipeline output
===============

adjacencies.csv    # gives strength of evidence for association between TF and targets
adjacencies_wCor.csv    # gives correlation to the associations between TF and targets
reg.csv    # modules of TFs and their target genes (not been refined to regulons)
aucell.csv   # activity of identified gene regulatory modules within each cell

"""

from ruffus import *
import sys
import os
from cgatcore import pipeline as P

PARAMS = P.get_parameters(
    ["%s/pipeline.yml" % os.path.splitext(__file__)[0],
      "../pipeline.yml",
      "pipeline.yml"])

#if not os.path.exists("data.dir/" + PARAMS["datatype"] + ".dir"):
#    os.makedirs("data.dir/" + PARAMS["datatype"] + ".dir")

@active_if(PARAMS["filtering_input_format"] == "csv")
@follows(mkdir("pyscenic_results.dir"))
@transform("data.dir/*_*-expression.csv", regex(r"data.dir/([^_]+)_(r.*|n.*)-expression.csv"), r"pyscenic_results.dir/\2.dir/\1.dir/filtered-expression.csv")
def gene_filtering(infile, outfile):
    '''
    Filtering of the input raw UMI expression matrix to remove genes expressed in few cells
    '''

    sample = infile.split('/')[1]
    sample = sample.split('_')[0]

    PY_PATH = os.path.join(os.path.dirname(os.path.dirname(__file__)), "python")

    statement = """python %(PY_PATH)s/filter_csv.py
                --input %(infile)s --umi_counts %(filtering_UMI_counts)s
                --min_percent %(filtering_min_percent)s
                --output %(outfile)s"""

    P.run(statement, job_threads = PARAMS["filtering_threads"], job_memory = '2G', job_queue = PARAMS["cluster_queue"])

# TODO: Fix script to work with AnnData objects
# @active_if(PARAMS["filtering_input_format"] == "h5ad")
# @follows(mkdir("pyscenic_results.dir"))
# @transform("data.dir/*_anndata.h5ad", regex(r"data.dir/([^_]+)_anndata.h5ad"), r"data.dir/\1_filtered-raw-expression.csv")
# def extract_anndata(infile, outfile):
#     '''
#     Extracts raw expression matrix from AnnData object
#     Filters to remove genes expressed in few cells
#     '''
#
#     if PARAMS["filtering_other_options"] == None:
#         filtering_other_options = ""
#     else:
#         filtering_other_options = PARAMS["filtering_other_options"]
#
#     PY_PATH = os.path.join(os.path.dirname(os.path.dirname(__file__)), "python")
#
#     statement = """python %(PY_PATH)s/extract_from_anndata.py
#                 --input %(infile)s --umi_counts %(filtering_UMI_counts)s
#                 --min_percent %(filtering_min_percent)s
#                 %(filtering_other_options)s
#                 --output %(outfile)s"""
#
#     P.run(statement, job_threads = PARAMS["filtering_threads"], job_memory = '10G', job_queue = PARAMS["cluster_queue"])

# @transform([gene_filtering, extract_anndata], regex(r"pyscenic_results.dir/(r.*|n.*).dir/([^_]+)_filtered-expression.csv"), r"pyscenic_results.dir/\1.dir/\2_adjacencies.tsv")
@transform(gene_filtering, regex(r"pyscenic_results.dir/(r.*|n.*).dir/([^_]+).dir/filtered-expression.csv"), r"pyscenic_results.dir/\1.dir/\2.dir/adjacencies.tsv")
def arboreto_with_multiprocessing(infile, outfile):
    '''
    Gene regulatory network inference
    Output contains co-expression modules - list of adjacencies between a TF and its targets
    '''

    if PARAMS["grn_other_options"] == None:
        grn_other_options = ""
    else:
        grn_other_options = PARAMS["grn_other_options"]

    PY_PATH = os.path.join(os.path.dirname(os.path.dirname(__file__)), "python")

    statement = """python %(PY_PATH)s/arboreto_with_multiprocessing.py
                %(infile)s %(grn_tfs_list)s
                --output %(outfile)s
                --num_workers %(grn_threads)s
                %(grn_other_options)s"""

    P.run(statement, job_threads = PARAMS["grn_threads"], job_memory = '10G', job_queue = PARAMS["cluster_queue"])

@transform(arboreto_with_multiprocessing, regex(r"pyscenic_results.dir/(r.*|n.*).dir/([^_]+).dir/adjacencies.tsv"), add_inputs(r"pyscenic_results.dir/\1.dir/\2.dir/filtered-expression.csv"), r"pyscenic_results.dir/\1.dir/\2.dir/adjacencies_wCor.csv")
def pyscenic_add_cor(infiles, outfile):
    '''
    Add Pearson correlations based on TF-gene expression to the network adjacencies output from the
    GRN step, and output these to a new adjacencies file
    '''

    adjacencies = infiles[0]
    expression_matrix = infiles[1]

    if PARAMS["add_cor_other_options"] == None:
        add_cor_other_options = ""
    else:
        add_cor_other_options = PARAMS["add_cor_other_options"]

    statement = """pyscenic add_cor
                %(adjacencies)s
                %(expression_matrix)s
                %(add_cor_other_options)s
                -o %(outfile)s"""

    P.run(statement, job_threads = PARAMS["add_cor_threads"], job_memory = '10G', job_queue = PARAMS["cluster_queue"])

@transform(pyscenic_add_cor, regex(r"pyscenic_results.dir/(r.*|n.*).dir/([^_]+).dir/adjacencies_wCor.csv"), add_inputs(r"pyscenic_results.dir/\1.dir/\2.dir/filtered-expression.csv"), r"pyscenic_results.dir/\1.dir/\2.dir/reg.csv")
def pyscenic_ctx(infiles, outfile):
    '''
    Regulons are derived from adjacencies
    '''

    adjacencies = infiles[0]
    expression_matrix = infiles[1]
    database_fname_1 = PARAMS["ctx_database_fname"][0]

    if len(PARAMS["ctx_database_fname"]) == 1:
        database_fname_2 = ""
    elif PARAMS["ctx_database_fname"][1] == None:
        database_fname_2 = ""
    else:
        database_fname_2 = PARAMS["ctx_database_fname"][1]

    if PARAMS["ctx_other_options"] == None:
        ctx_other_options = ""
    else:
        ctx_other_options = PARAMS["ctx_other_options"]

    statement = """pyscenic ctx
                -o %(outfile)s
                --annotations_fname %(ctx_annotations)s
                --num_workers %(ctx_threads)s
                --expression_mtx_fname %(expression_matrix)s
                %(ctx_other_options)s
                %(adjacencies)s
                %(database_fname_1)s
                %(database_fname_2)s"""

    P.run(statement, job_threads = PARAMS["ctx_threads"], job_memory = '10G', job_queue = PARAMS["cluster_queue"])

@transform(pyscenic_ctx, regex(r"pyscenic_results.dir/(r.*|n.*).dir/([^_]+).dir/reg.csv"), add_inputs(r"pyscenic_results.dir/\1.dir/\2.dir/filtered-expression.csv"), r"pyscenic_results.dir/\1.dir/\2.dir/aucell.csv")
def pyscenic_aucell(infiles, outfile):
    '''
    Characterise cells by enrichment of previously discovered regulons
    Enrichment of a regulon is measured as the Area Under the recovery Curve (AUC) of the genes that define this regulon
    '''

    reg = infiles[0]
    expression_matrix = infiles[1]

    if PARAMS["aucell_other_options"] == None:
        aucell_other_options = ""
    else:
        aucell_other_options = PARAMS["aucell_other_options"]

    statement = """pyscenic aucell
                -o %(outfile)s
                --num_workers %(aucell_threads)s
                --auc_threshold %(aucell_auc_threshold)s
                %(aucell_other_options)s
                %(expression_matrix)s
                %(reg)s"""

    P.run(statement, job_threads = PARAMS["aucell_threads"], job_memory = '10G', job_queue = PARAMS["cluster_queue"])

@follows(pyscenic_aucell)
def full():
    pass

def main(argv = None):
    if argv is None:
        argv = sys.argv
    P.main(argv)

if __name__ == "__main__":
    sys.exit(P.main(sys.argv))
