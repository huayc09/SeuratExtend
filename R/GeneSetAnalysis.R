#' @title Gene Set Enrichment Analysis
#' @description Calculate the GSEA score of gene sets at the single-cell level using the 'AUCell' package:
#'
#' Aibar et al. (2017) SCENIC: Single-cell regulatory network inference and clustering.
#' Nature Methods. doi: 10.1038/nmeth.4463
#'
#' Aibar et al. (2016) AUCell: Analysis of 'gene set' activity in single-cell RNA-seq data.
#' R/Bioconductor package.
#'
#' This function can utilize either customized gene sets or the pre-built GO or Reactome database.
#' @param seu Seurat object.
#' @param genesets List of gene sets (or signatures) to test in the cells.
#' The gene sets should be provided as a character list.
#' @param title Name of the slot where the data will be saved.
#' @param ratio Minimum ratio of genes in the gene set detected in the datasets. Default: 0.4.
#' @param n.min Minimum number of genes in the gene set. Default: 1.
#' @param n.max Maximum number of genes in the gene set. Default: Inf.
#' @param slot Slot from which to pull feature data. Default: 'counts'.
#' @param assay Name of the assay to use. Default: 'RNA'.
#' @param nCores Number of cores to use for computation. Default: 1.
#' @param aucMaxRank Threshold for calculating the AUC.
#' @param export_to_matrix If TRUE, the function will return an AUCell matrix instead of a Seurat object.
#' @param verbose Should the function display progress messages? Default: TRUE.
#' @param n.items.part If the datasets/gene sets are too large, the function can split the gene sets into n parts
#' to reduce RAM usage.
#' @return Seurat object or matrix.
#' @details If returning a Seurat object, the AUCell matrix is saved in seu@misc[["AUCell"]][[title]].
#' @examples
#' library(SeuratExtend)
#' options(spe = "human")
#'
#' # Perform GSEA using the Gene Ontology (GO) database. Given the extensive size of the
#' # entire database, this example only evaluates pathways under the "immune_system_process"
#' # category. The results will be saved in: seu@misc$AUCell$GO[[title]]
#' pbmc <- GeneSetAnalysisGO(pbmc, parent = "immune_system_process")
#' matr <- pbmc@misc$AUCell$GO$immune_system_process
#' matr <- RenameGO(matr)
#' head(matr, 4:3)
#'
#' # For the "parent" argument, you can use any term from the GO database, be it a GO ID or
#' # pathway name. Using `GeneSetAnalysisGO()` without arguments will show commonly used
#' # GO categories:
#' GeneSetAnalysisGO()
#'
#' # To visualize the data, consider using a heatmap (to compare multiple groups with more
#' # features, albeit with less detailed representation), a violin plot (to compare multiple
#' # groups with fewer features, but presenting more details for individual data points), or a
#' # waterfall plot (to contrast only two groups). Below is an example of a heatmap:
#' Heatmap(CalcStats(matr, f = pbmc$cluster, order = "p", n = 4), lab_fill = "zscore")
#'
#' # Violin plot example:
#' VlnPlot2(matr[1:3,], f = pbmc$cluster)
#'
#' # Waterfall plot example:
#' WaterfallPlot(matr, f = pbmc$cluster, ident.1 = "B cell", ident.2 = "CD8 T cell", top.n = 20)
#'
#' # Conduct GSEA using the Reactome database. This example will only assess pathways under
#' # the "Immune System" category. Results will be stored in: seu@misc$AUCell$Reactome[[title]]
#' pbmc <- GeneSetAnalysisReactome(pbmc, parent = "Immune System")
#' matr <- pbmc@misc$AUCell$Reactome$`Immune System`
#' matr <- RenameReactome(matr)
#' Heatmap(CalcStats(matr, f = pbmc$cluster, order = "p", n = 4), lab_fill = "zscore")
#'
#' # As with GO, you can run `GeneSetAnalysisReactome()` without arguments to view
#' # commonly used categories in the Reactome database:
#' GeneSetAnalysisReactome()
#'
#' # For GSEA using custom gene sets, the output AUCell matrix will be saved under:
#' # seu@misc$AUCell[[title]]
#' pbmc <- GeneSetAnalysis(pbmc, genesets = hall50$human)
#' matr <- pbmc@misc$AUCell$genesets
#' Heatmap(CalcStats(matr, f = pbmc$cluster), lab_fill = "zscore")
#' @rdname GeneSetAnalysis
#' @export

GeneSetAnalysis <- function(
    seu = NULL,
    genesets,
    title = "genesets",
    ratio = 0.4, n.min = 1, n.max = Inf,
    slot = "counts", assay = "RNA", nCores = 1,
    aucMaxRank = NULL,
    export_to_matrix = F, verbose = TRUE,
    n.items.part = NULL) {

  if(!require(SeuratObject)) library(Seurat)
  DefaultAssay(seu) <- assay
  genesets <- FilterGenesets(
    genes_in_data = rownames(seu),
    genesets = genesets,
    ratio = ratio,
    n.min = n.min,
    n.max = n.max,
    verbose = verbose)

  seu <- BuildAUCRank(
    seu,
    slot = slot,
    assay = assay,
    verbose = verbose)

  if(is.null(aucMaxRank)) aucMaxRank <- ceiling(0.05 * nrow(seu))
  AUC_matrix <- calcAUC_matrix(
    GenesetList = genesets,
    rankings = seu@misc$AUCell$cells_rankings,
    nCores = nCores,
    aucMaxRank = aucMaxRank,
    verbose = verbose,
    n.items.part = n.items.part)

  if(export_to_matrix) return(AUC_matrix)
  seu@misc[["AUCell"]][[title]] <- AUC_matrix
  return(seu)
}

# Internal ----------------------------------------------------------------

FilterGenesets <- function(
    genes_in_data,
    genesets,
    ratio = 0.4,
    n.min = 1,
    n.max = Inf,
    verbose = TRUE) {
  if(verbose) {
    message(paste0(
      Sys.time(),
      " Start filtering ",length(genesets)," gene set(s): "),
      "n(Genes) >= ", n.min, ", n(Genes) <= ", n.max,
      ", at least ", ratio * 100, "% of genes found in the datasets")
  }
  filter <- sapply(
    genesets,
    function(x){
      length(x) >= n.min &
        length(x) <= n.max &
        sum(genes_in_data %in% x)/length(x) > ratio
    })
  if(!any(filter)) stop("No gene set(s) pass the filter. Please check inputs.")
  genesets <- genesets[filter]
  if(verbose) {
    message(paste(
      Sys.time(),
      length(genesets),"gene set(s) passed the filter"))
  }
  return(genesets)
}
