"""
=============================
Pipeline pySCENIC R analysis
=============================
Author: Devika Agarwal
Overview
========
This pipeline performs pySCENIC downstream analysis steps (https://pyscenic.readthedocs.io/en/latest/index.html)
including:
1. AUCell score based UMAP projection and clustering in seurat
2. AUCell score based heatmaps and violin plots
3. Binarized score projection on UMAP for each regulon
4. Wilcoxon and KS statistical test for each regulon for celltype and condition
5. Average AUCell score heatmaps based on wilcoxon /KS test results
6. Pathway analyses and plots for top regulons per condition/celltype for regulon target genes
"""

from ruffus import *
import sys
import os
from cgatcore import pipeline as P

PARAMS = P.get_parameters(
    ["%s/pipeline.yml" % os.path.splitext(__file__)[0],
      "../pipeline.yml",
      "pipeline.yml"])

@follows(mkdir("reports.dir"))
@transform("pyscenic_results.dir/*_aucell.csv", regex(r"pyscenic_results.dir/([^_]+)_aucell.csv"), r"reports.dir/scenic_seurat.html")
def scenic_seurat(infile, outfile):
    ''' Seurat based analyses for scenic results'''

    if PARAMS["datatype"] == "raw":
        outfile = outfile.replace("scenic_seurat.html", "scenic_seurat_raw.html")
    elif PARAMS["datatype"] == "normalised":
        outfile = outfile.replace("scenic_seurat.html", "scenic_seurat_normalised.html")

    R_PATH = os.path.join(os.path.dirname(os.path.dirname(__file__)), "R")
    outbase = P.snip(outfile, ".html")
    working_dir = PARAMS["working_dir"]
    seurat_object = PARAMS["rseurat_seurat_object"]
    aucell_zscores = PARAMS["aucell_zscores"]
    aucell_zscores_celltype = PARAMS["aucell_zscores_celltype"]
    aucell_zscores_condition = PARAMS["aucell_zscores_condition"]
    binary_matrix = PARAMS["binary_matrix"]
    results_directory = PARAMS["results_directory"]
    plots_directory = PARAMS["plots_directory"]
    datatype = PARAMS["datatype"]
    umap_pcs = PARAMS["rseurat_umap_pcs"]
    cell_type = PARAMS["rseurat_cell_type"]
    condition = PARAMS["rseurat_condition"]
    clustering_resolution = PARAMS["rseurat_clustering_resolution"]
    stacked_vln_function = PARAMS["rseurat_stacked_vln_function"]
    reference_condition = PARAMS["rseurat_reference_condition"]
    FDR_threshold = PARAMS["rseurat_FDR_threshold"]

    statement = """ Rscript -e "rmarkdown::render('%(R_PATH)s/scenic_seurat.Rmd',
                        params = list(  working_dir = '%(working_dir)s',
                                        seurat_object = '%(working_dir)s/%(seurat_object)s',
                                        aucell_scores = '%(working_dir)s/%(infile)s',
                                        aucell_zscores = '%(working_dir)s/%(aucell_zscores)s',
                                        aucell_zscores_celltype = '%(working_dir)s/%(aucell_zscores_celltype)s',
                                        aucell_zscores_condition = '%(working_dir)s/%(aucell_zscores_condition)s',
                                        binary_matrix = '%(working_dir)s/%(binary_matrix)s',
                                        results_directory = '%(results_directory)s',
                                        plots_directory = '%(plots_directory)s',
                                        datatype = '%(datatype)s',
                                        umap_pcs = '%(umap_pcs)s',
                                        cell_type = '%(cell_type)s',
                                        condition = '%(condition)s',
                                        clustering_resolution = '%(clustering_resolution)s',
                                        stacked_vln_function = '%(working_dir)s/%(stacked_vln_function)s',
                                        reference_condition = '%(reference_condition)s' ,
                                        FDR_threshold = '%(FDR_threshold)s'),
                                        output_file = '%(working_dir)s/%(outfile)s')"
                        > %(outbase)s.log
                        2> %(outbase)s.err """

    P.run(statement, job_threads = PARAMS["rseurat_threads"], job_memory = '10G', job_queue = PARAMS["cluster_queue"])

@follows(scenic_seurat)
@transform("pyscenic_results.dir/binary_matrix.csv", regex(r"pyscenic_results.dir/binary_matrix.csv"), r"reports.dir/Scenic_analysis_R.html")
def rscenic(infile, outfile):
    ''' R based analyses for scenic results'''
    if PARAMS["datatype"] == "raw":
        outfile = outfile.replace("Scenic_analysis_R.html", "Scenic_analysis_R_raw.html")
    elif PARAMS["datatype"] == "normalised":
        outfile = outfile.replace("Scenic_analysis_R.html", "Scenic_analysis_R_normalised.html")

    R_PATH = os.path.join(os.path.dirname(os.path.dirname(__file__)), "R")
    outbase = P.snip(outfile, ".html")
    working_dir = PARAMS["working_dir"]
    num_workers = PARAMS["rscenic_num_workers"]
    results_directory = PARAMS["results_directory"]
    plots_directory = PARAMS["plots_directory"]
    datatype = PARAMS["datatype"]
    celltype = PARAMS["rscenic_celltype"]
    condition = PARAMS["rscenic_condition"]
    celltype_condition = PARAMS["rscenic_celltype_condition"]
    celltype_zscores = PARAMS["aucell_zscores_celltype"]
    celltype_stim_zscores = PARAMS["rscenic_celltype_stim_zscores"]
    stim_zscores = PARAMS["aucell_zscores_condition"]
    celltype_zscore_filter_threshold = PARAMS["rscenic_celltype_zscore_filter_threshold"]
    stim_zscore_filter_threshold = PARAMS["rscenic_stim_zscore_filter_threshold"]
    celltype_stim_zscore_filter_threshold = PARAMS["rscenic_celltype_stim_zscore_filter_threshold"]
    binary_mtx = PARAMS["binary_matrix"]
    annotation_celltype = PARAMS["rscenic_annotation_celltype"]
    annotation_stim = PARAMS["rscenic_annotation_stim"]
    annotation_celltype_stim = PARAMS["rscenic_annotation_celltype_stim"]
    binary_heatmap_cell_annotations = PARAMS["rscenic_binary_heatmap_cell_annotations"]
    rss_celltype = PARAMS["rscenic_rss_celltype"]
    rss_stim = PARAMS["rscenic_rss_stim"]
    rss_celltype_stim = PARAMS["rscenic_rss_celltype_stim"]
    species = PARAMS["rscenic_species"]
    regulon_genes = PARAMS["rscenic_regulon_genes"]
    exp_matrix = PARAMS["rscenic_exp_matrix"]
    go_ont_option = PARAMS["rscenic_go_ont_option"]
    wilcoxon_celltype_top = PARAMS["rscenic_wilcoxon_celltype_top"]
    wilcoxon_condition_top = PARAMS["rscenic_wilcoxon_condition_top"]
    ks_celltype_top = PARAMS["rscenic_ks_celltype_top"]
    ks_condition_top = PARAMS["rscenic_ks_condition_top"]
    pvaluecutoff = PARAMS["rscenic_pvaluecutoff"]
    qvaluecutoff = PARAMS["rscenic_qvaluecutoff"]
    maxGSSize = PARAMS["rscenic_maxGSSize"]
    msigdb_geneset = PARAMS["rscenic_msigdb_geneset"]

    statement = """ Rscript -e "rmarkdown::render('%(R_PATH)s/scenic_seurat.Rmd',
                            params = list(  working_dir = '%(working_dir)s',
                                            num_workers = '%(num_workers)s',
                                            datatype = '%(datatype)s',
                                            results_directory= '%(results_directory)s',
                                            plots_directory = '%(plots_directory)s',
                                            celltype  = '%(celltype)s',
                                            condition = '%(condition)s',
                                            celltype_condition = '%(celltype_condition)s',
                                            celltype_zscores = '%(working_dir)s/%(celltype_zscores)s',
                                            celltype_stim_zscores = '%(working_dir)s/%(celltype_stim_zscores)s',
                                            stim_zscores = '%(working_dir)s/%(stim_zscores)s',
                                            celltype_zscore_filter_threshold = '%(celltype_zscore_filter_threshold)s',
                                            stim_zscore_filter_threshold = '%(stim_zscore_filter_threshold)s',
                                            celltype_stim_zscore_filter_threshold = '%(celltype_stim_zscore_filter_threshold)s',
                                            binary_mtx = '%(working_dir)s/%(binary_mtx)s',
                                            annotation_celltype = '%(working_dir)s/%(annotation_celltype)s',
                                            annotation_stim = '%(working_dir)s/%(annotation_stim)s',
                                            annotation_celltype_stim = '%(working_dir)s/%(annotation_celltype_stim)s',
                                            binary_heatmap_cell_annotations ='%(working_dir)s/%(binary_heatmap_cell_annotations)s',
                                            rss_celltype = '%(working_dir)s/%(rss_celltype)s',
                                            rss_stim = '%(working_dir)s/%(rss_stim)s',
                                            rss_celltype_stim = '%(working_dir)s/%(rss_celltype_stim)s',
                                            species = '%(species)s',
                                            regulon_genes = '%(working_dir)s/%(regulon_genes)s',
                                            exp_matrix = '%(working_dir)s/%(exp_matrix)s',
                                            go_ont_option ='%(go_ont_option)s',
                                            wilcoxon_celltype_top = '%(wilcoxon_celltype_top)s',
                                            wilcoxon_condition_top = '%(wilcoxon_condition_top)s',
                                            ks_celltype_top = '%(ks_celltype_top)s',
                                            ks_condition_top = '%(ks_condition_top)s',
                                            pvaluecutoff = '%(pvaluecutoff)s',
                                            qvaluecutoff = '%(qvaluecutoff)s',
                                            maxGSSize = '%(maxGSSize)s',
                                            msigdb_geneset = '%(msigdb_geneset)s'),
                                            output_file = '%(working_dir)s/%(outfile)s')"
                            > %(outbase)s.log
                            2> %(outbase)s.err """

    P.run(statement, job_threads = PARAMS["rscenic_threads"], job_memory = '10G', job_queue = PARAMS["cluster_queue"])

@follows(scenic_seurat,rscenic)
def full():
    pass

def main(argv=None):
    if argv is None:
        argv = sys.argv
    P.main(argv)

if __name__ == "__main__":
    sys.exit(P.main(sys.argv))
