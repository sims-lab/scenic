"""
=============================
Pipeline pySCENIC downstream
=============================

Author: Lucy Garner


Overview
========
This pipeline performs pySCENIC downstream analysis steps (https://pyscenic.readthedocs.io/en/latest/index.html)
including:
1. AUCell score heatmap
2. Generation of regulons
3. Binarization of regulons
4. Calculation of regulon specificity scores
5. Calculation of AUCell z-scores for individual cells and cell groups (e.g. clusters)

"""

from ruffus import *
import sys
import os
from cgatcore import pipeline as P

PARAMS = P.get_parameters(
    ["%s/pipeline.yml" % os.path.splitext(__file__)[0],
      "../pipeline.yml",
      "pipeline.yml"])

@follows(mkdir("plots.dir"))
@transform("pyscenic_results.dir/*_aucell.csv", regex(r"pyscenic_results.dir/([^_]+)_aucell.csv"), r"plots.dir/aucell_heatmap.png")
def aucell_heatmap(infile, outfile):

    exp_mtx = infile.replace("pyscenic_results.dir", "data.dir")
    if PARAMS["datatype"] == "raw":
        exp_mtx = exp_mtx.replace("aucell.csv", "filtered-raw-expression.csv")
    elif PARAMS["datatype"] == "normalised":
        exp_mtx = exp_mtx.replace("aucell.csv", "filtered-normalised-expression.csv")

    PY_PATH = os.path.join(os.path.dirname(os.path.dirname(__file__)), "python")

    statement = """python %(PY_PATH)s/aucell_heatmap.py
                --exp_mtx %(exp_mtx)s
                --aucell_output %(infile)s
                %(aucell_tab)s
                """

    P.run(statement, job_threads = PARAMS["aucell_threads"], job_memory = '10G', job_queue = PARAMS["cluster_queue"])

@transform("pyscenic_results.dir/*_reg.csv", regex(r"pyscenic_results.dir/([^_]+)_reg.csv"), r"pyscenic_results.dir/regulons.csv")
def generate_regulons(infile, outfile):

    PY_PATH = os.path.join(os.path.dirname(os.path.dirname(__file__)), "python")

    statement = """python %(PY_PATH)s/generate_regulons.py
                --ctx_output %(infile)s
                """

    P.run(statement, job_threads = PARAMS["regulons_threads"], job_memory = '10G', job_queue = PARAMS["cluster_queue"])

@transform("pyscenic_results.dir/*_aucell.csv", regex(r"pyscenic_results.dir/([^_]+)_aucell.csv"), r"pyscenic_results.dir/aucell_thresholds.csv")
def regulon_binarization(infile, outfile):

    if PARAMS["binarize_tab"] == None:
        binarize_tab = ""
    else:
        binarize_tab = PARAMS["binarize_tab"]

    if PARAMS["binarize_custom_aucell_thresholds"] == None:
        custom_aucell_thresholds = ""
    else:
        custom_aucell_thresholds = '--custom_auc_thresholds ' + PARAMS["binarize_custom_aucell_thresholds"]

    PY_PATH = os.path.join(os.path.dirname(os.path.dirname(__file__)), "python")

    statement = """python %(PY_PATH)s/regulon_binarization.py
                --aucell_output %(infile)s
                %(binarize_tab)s
                %(custom_aucell_thresholds)s
                """

    P.run(statement, job_threads = PARAMS["binarize_threads"], job_memory = '10G', job_queue = PARAMS["cluster_queue"])

@transform("pyscenic_results.dir/*_aucell.csv", regex(r"pyscenic_results.dir/([^_]+)_aucell.csv"), r"pyscenic_results.dir/aucell_zscores.csv")
def rss_zscore(infile, outfile):

    exp_mtx = infile.replace("pyscenic_results.dir", "data.dir")
    if PARAMS["datatype"] == "raw":
        exp_mtx = exp_mtx.replace("aucell.csv", "filtered-raw-expression.csv")
    elif PARAMS["datatype"] == "normalised":
        exp_mtx = exp_mtx.replace("aucell.csv", "filtered-normalised-expression.csv")

    if PARAMS["rss_zscore_tab"] == None:
        rss_zscore_tab = ""
    else:
        rss_zscore_tab = PARAMS["rss_zscore_tab"]

    PY_PATH = os.path.join(os.path.dirname(os.path.dirname(__file__)), "python")

    statement = """python %(PY_PATH)s/rss_zscore.py
                --exp_mtx %(exp_mtx)s
                --aucell_output %(infile)s
                --annotation_files %(rss_zscore_annotation_files)s
                %(rss_zscore_tab)s
                """

    P.run(statement, job_threads = PARAMS["rss_zscore_threads"], job_memory = '2G', job_queue = PARAMS["cluster_queue"])

@follows(aucell_heatmap, generate_regulons, regulon_binarization, rss_zscore)
def full():
    pass

def main(argv=None):
    if argv is None:
        argv = sys.argv
    P.main(argv)

if __name__ == "__main__":
    sys.exit(P.main(sys.argv))
