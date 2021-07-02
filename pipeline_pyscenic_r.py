"""
=============================
Pipeline pySCENIC R analysis
=============================

Authors: Devika Agarwal and Lucy Garner

Overview
========
This pipeline performs pySCENIC downstream analysis steps (https://pyscenic.readthedocs.io/en/latest/index.html)
including:
1. AUCell score based UMAP projection and clustering in seurat
2. AUCell score based heatmaps and violin plots
3. Binarized score projection on UMAP for each regulon
4. Wilcoxon (or chosen test) and KS statistical test for each regulon for celltype and condition
5. Average AUCell score heatmaps based on wilcoxon (ot chosen test)/KS test results
6. Pathway analyses and plots for top regulons per condition/celltype for regulon target genes

"""

from ruffus import *
import sys
import os
import re
from cgatcore import pipeline as P

PARAMS = P.get_parameters(
    ["%s/pipeline.yml" % os.path.splitext(__file__)[0],
      "../pipeline.yml",
      "pipeline.yml"])

@follows(mkdir("reports.dir"))
@transform("pyscenic_results.dir/*.dir/*.dir/aucell.csv",
           regex(r"pyscenic_results.dir/(r.*|n.*).dir/([^_]+).dir/aucell.csv"),
           r"reports.dir/\1.dir/\2.dir/scenic_seurat.html")
def scenic_seurat(infile, outfile):
    ''' Seurat based analyses for scenic results'''

    if "normalised" in infile:
        datatype = "normalised"
    elif "raw" in infile:
        datatype = "raw"

    sample = infile.split("/")[2]
    sample = sample.replace(".dir", "")

    R_PATH = os.path.join(os.getcwd(), "R")
    outbase = P.snip(outfile, ".html")
    dir = re.sub("aucell.csv", "", infile)

    aucell_zscores = dir + PARAMS["aucell_zscores"]
    aucell_zscores_celltype = dir + "aucell_zscores_" + PARAMS["annotation_celltype"] + "-annotation.csv"
    aucell_zscores_condition = dir + "aucell_zscores_" + PARAMS["annotation_condition"] + "-annotation.csv"
    binary_matrix = dir + PARAMS["binary_matrix"]
    working_dir = PARAMS["working_dir"]
    seurat_object = "data.dir/" + sample + "_" + PARAMS["rseurat_seurat_object"]
    results_directory = "pyscenic_results.dir/" + datatype + ".dir/" + sample + ".dir/pyscenic_r.dir"
    plots_directory = "plots.dir/" + datatype + ".dir/" + sample + ".dir/pyscenic_r.dir"
    umap_pcs = PARAMS["rseurat_umap_pcs"]
    celltype = PARAMS["rseurat_celltype"]
    clustering_resolution = PARAMS["rseurat_clustering_resolution"]
    stacked_vln_function = PARAMS["rseurat_stacked_vln_function"]
    dotplot_function = PARAMS["rseurat_dotplot_function"]
    diff_exp_test = PARAMS["rseurat_diff_exp_test"]
    latent_variables = PARAMS["rseurat_latent_variables"]
    reference_condition = PARAMS["rseurat_reference_condition"]
    celltype_condition = PARAMS["rseurat_celltype_condition"]
    FDR_threshold = PARAMS["rseurat_FDR_threshold"]

    condition = PARAMS["rseurat_condition"]
    condition_exclusion = PARAMS["rseurat_condition_exclusion"]
    if condition_exclusion != "None" and sample in condition_exclusion:
        condition = "None"

    os.makedirs(results_directory)

    statement = """Rscript -e "rmarkdown::render('%(R_PATH)s/scenic_seurat.Rmd',
                                                 params = list(working_dir = '%(working_dir)s',
                                                 seurat_object = '%(working_dir)s/%(seurat_object)s',
                                                 aucell_scores = '%(working_dir)s/%(infile)s',
                                                 aucell_zscores = '%(working_dir)s/%(aucell_zscores)s',
                                                 aucell_zscores_celltype = '%(working_dir)s/%(aucell_zscores_celltype)s',
                                                 aucell_zscores_condition = '%(working_dir)s/%(aucell_zscores_condition)s',
                                                 binary_matrix = '%(working_dir)s/%(binary_matrix)s',
                                                 results_directory = '%(working_dir)s/%(results_directory)s',
                                                 plots_directory = '%(working_dir)s/%(plots_directory)s',
                                                 datatype = '%(datatype)s',
                                                 umap_pcs = '%(umap_pcs)s',
                                                 celltype = '%(celltype)s',
                                                 condition = '%(condition)s',
                                                 clustering_resolution = '%(clustering_resolution)s',
                                                 stacked_vln_function = '%(working_dir)s/%(stacked_vln_function)s',
                                                 dotplot_function = '%(working_dir)s/%(dotplot_function)s',
                                                 diff_exp_test = '%(diff_exp_test)s',
                                                 latent_variables = '%(latent_variables)s',
                                                 reference_condition = '%(reference_condition)s',
                                                 celltype_condition = '%(celltype_condition)s',
                                                 FDR_threshold = '%(FDR_threshold)s'),
                                                 output_file = '%(working_dir)s/%(outfile)s')"
                   > %(working_dir)s/%(results_directory)s/scenic_seurat.log
                   2> %(working_dir)s/%(results_directory)s/scenic_seurat.err"""

    P.run(statement, job_threads = PARAMS["rseurat_threads"], job_memory = '10G',
          job_queue = PARAMS["cluster_queue"], job_condaenv = PARAMS["conda_env"])

@follows(scenic_seurat)
@transform("pyscenic_results.dir/*.dir/*.dir/binary_matrix.csv",
           regex(r"pyscenic_results.dir/(r.*|n.*).dir/([^_]+).dir/binary_matrix.csv"),
           r"reports.dir/\1.dir/\2.dir/scenic_analysis_R.html")
def rscenic(infile, outfile):
    ''' R based analyses for scenic results'''

    if "normalised" in infile:
        datatype = "normalised"
    elif "raw" in infile:
        datatype = "raw"

    sample = infile.split("/")[2]
    sample = sample.replace(".dir", "")

    R_PATH = os.path.join(os.getcwd(), "R")
    outbase = P.snip(outfile, ".html")
    dir = re.sub("binary_matrix.csv", "", infile)

    exp_matrix = dir + "filtered-expression.csv"
    celltype = PARAMS["rscenic_celltype"]
    condition = PARAMS["rscenic_condition"]
    condition_exclusion = PARAMS["rscenic_condition_exclusion"]
    if condition_exclusion != "None" and sample in condition_exclusion:
        condition = "None"
    celltype_condition = PARAMS["rscenic_celltype_condition"]

    annotation_celltype = PARAMS["rscenic_annotation_celltype"]
    annotation_condition = PARAMS["rscenic_annotation_condition"]
    annotation_celltype_condition = "data.dir/" + sample + "_" + PARAMS["rscenic_annotation_celltype_condition"]

    zscores_celltype = dir + "aucell_zscores_" + PARAMS["rscenic_annotation_celltype"]
    rss_celltype = dir + PARAMS["rscenic_annotation_celltype"].split(".")[0] + "_RSS.csv"
    zscore_filter_threshold_celltype = PARAMS["rscenic_zscore_filter_threshold_celltype"]

    celltype_top = PARAMS["rseurat_diff_exp_test"] + '_celltype_top10.csv'
    if condition != "None":
        zscores_condition = dir + "aucell_zscores_" + PARAMS["rscenic_annotation_condition"]
        rss_condition = dir + PARAMS["rscenic_annotation_condition"].split(".")[0] + "_RSS.csv"
        zscore_filter_threshold_condition = PARAMS["rscenic_zscore_filter_threshold_condition"]
        condition_reference_top = PARAMS["rseurat_diff_exp_test"] + '_condition_reference_top10.csv'
        condition_pairwise_top = PARAMS["rseurat_diff_exp_test"] + '_condition_pairwise_top10.csv'
        ks_condition_top = 'ks_condition_top10.csv'
    else:
        zscores_condition = "None"
        rss_condition = "None"
        zscore_filter_threshold_condition = "None"
        condition_reference_top = "None"
        condition_pairwise_top = "None"
        ks_condition_top = "None"

    zscores_celltype_condition = dir + PARAMS["rscenic_zscores_celltype_condition"]
    rss_celltype_condition = dir + PARAMS["rscenic_rss_celltype_condition"]
    zscore_filter_threshold_celltype_condition = PARAMS["rscenic_zscore_filter_threshold_celltype_condition"]

    regulon_genes = dir + "regulons.csv"
    species = PARAMS["rscenic_species"]
    working_dir = PARAMS["working_dir"]
    num_workers = PARAMS["rscenic_num_workers"]

    results_directory = "pyscenic_results.dir/" + datatype + ".dir/" + sample + ".dir/pyscenic_r.dir"
    plots_directory = "plots.dir/" + datatype + ".dir/" + sample + ".dir/pyscenic_r.dir"

    binary_heatmap_cell_annotations = sample + "_" + annotation_celltype
    go_ont_option = PARAMS["rscenic_go_ont_option"]
    pvaluecutoff = PARAMS["rscenic_pvaluecutoff"]
    qvaluecutoff = PARAMS["rscenic_qvaluecutoff"]
    maxGSSize = PARAMS["rscenic_maxGSSize"]
    msigdb_geneset = PARAMS["rscenic_msigdb_geneset"]

    statement = """Rscript -e "rmarkdown::render('%(R_PATH)s/scenic_analysis_R.Rmd',
                                                 params = list(working_dir = '%(working_dir)s',
                                                 num_workers = '%(num_workers)s',
                                                 datatype = '%(datatype)s',
                                                 results_directory = '%(working_dir)s/%(results_directory)s',
                                                 plots_directory = '%(working_dir)s/%(plots_directory)s',
                                                 celltype  = '%(celltype)s',
                                                 condition = '%(condition)s',
                                                 celltype_condition = '%(celltype_condition)s',
                                                 zscores_celltype = '%(working_dir)s/%(zscores_celltype)s',
                                                 zscores_condition = '%(working_dir)s/%(zscores_condition)s',
                                                 zscores_celltype_condition = '%(working_dir)s/%(zscores_celltype_condition)s',
                                                 zscore_filter_threshold_celltype = '%(zscore_filter_threshold_celltype)s',
                                                 zscore_filter_threshold_condition = '%(zscore_filter_threshold_condition)s',
                                                 zscore_filter_threshold_celltype_condition = '%(zscore_filter_threshold_celltype_condition)s',
                                                 binary_mtx = '%(working_dir)s/%(infile)s',
                                                 annotation_celltype = '%(working_dir)s/data.dir/%(sample)s_%(annotation_celltype)s',
                                                 annotation_condition = '%(working_dir)s/data.dir/%(sample)s_%(annotation_condition)s',
                                                 annotation_celltype_condition = '%(working_dir)s/%(annotation_celltype_condition)s',
                                                 binary_heatmap_cell_annotations ='%(working_dir)s/data.dir/%(binary_heatmap_cell_annotations)s',
                                                 rss_celltype = '%(working_dir)s/%(rss_celltype)s',
                                                 rss_condition = '%(working_dir)s/%(rss_condition)s',
                                                 rss_celltype_condition = '%(working_dir)s/%(rss_celltype_condition)s',
                                                 species = '%(species)s',
                                                 regulon_genes = '%(working_dir)s/%(regulon_genes)s',
                                                 exp_matrix = '%(working_dir)s/%(exp_matrix)s',
                                                 go_ont_option = '%(go_ont_option)s',
                                                 diff_exp_test = '%(rseurat_diff_exp_test)s',
                                                 celltype_top = '%(celltype_top)s',
                                                 condition_reference_top = '%(condition_reference_top)s',
                                                 condition_pairwise_top = '%(condition_pairwise_top)s',
                                                 ks_celltype_top = 'ks_celltype_top10.csv',
                                                 ks_condition_top = '%(ks_condition_top)s',
                                                 pvaluecutoff = '%(pvaluecutoff)s',
                                                 qvaluecutoff = '%(qvaluecutoff)s',
                                                 maxGSSize = '%(maxGSSize)s',
                                                 msigdb_geneset = '%(msigdb_geneset)s'),
                                                 output_file = '%(working_dir)s/%(outfile)s')"
                   > %(working_dir)s/%(results_directory)s/scenic_analysis_R.log
                   2> %(working_dir)s/%(results_directory)s/scenic_analysis_R.err"""

    P.run(statement, job_threads = PARAMS["rscenic_num_workers"], job_memory = '15G',
          job_queue = PARAMS["cluster_queue"], job_condaenv = PARAMS["conda_env"])

@follows(scenic_seurat, rscenic)
def full():
    pass

def main(argv = None):
    if argv is None:
        argv = sys.argv
    P.main(argv)

if __name__ == "__main__":
    sys.exit(P.main(sys.argv))
