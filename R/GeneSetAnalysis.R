#' @title Gene Set Enrichment Analysis
#' @description Calculate GSEA score of gene sets at single cell level, By 'AUCell' package:
#'
#' Aibar et al. (2017) SCENIC: single-cell regulatory network inference and clustering.
#' Nature Methods. doi: 10.1038/nmeth.4463
#'
#' Aibar. et al. (2016) AUCell: Analysis of 'gene set' activity in single-cell RNA-seq data.
#' R/Bioconductor package.
#'
#' Can use either customized genesets, or pre-built GO or Reactome database
#' @param seu Seurat object
#' @param genesets List of gene-sets (or signatures) to test in the cells.
#' The gene-sets should be provided as character list.
#' @param title Name of slot where the data will be saved
#' @param ratio Minimum ratio of genes in the geneset detected in the datasets, Default: 0.4
#' @param n.min Min number of genes in the geneset, Default: 1
#' @param n.max Max number of genes in the geneset, Default: Inf
#' @param slot Slot to pull feature data for, Default: 'counts'
#' @param assay Name of assay to use, Default: 'RNA'
#' @param nCores Number of cores to use for computation. Default: 1
#' @param aucMaxRank Threshold to calculate the AUC
#' @param export_to_matrix If TRUE, then return a AUCell matrix instead of Seurat object
#' @param verbose Should the function show progress messages? Default: TRUE
#' @param n.items.part If the datasets/genesets are huge, split the genesets into n parts to
#' let it run with lower RAM
#' @return Seurat object or Matrix
#' @details If return Seurat object, AUCell matrix is saved in seu@misc[["AUCell"]][[title]]
#' @examples
#' options(spe = "human")
#'
#' # Geneset enrichment analysis (GSEA) using customized genesets
#' pbmc <- GeneSetAnalysis(pbmc, genesets = hall50$human)
#' matr <- pbmc@misc$AUCell$genesets
#' Heatmap(CalcStats(matr, f = pbmc$cluster), lab_fill = "zscore")
#'
#' # GSEA using GO database
#' pbmc <- GeneSetAnalysisGO(pbmc, parent = "immune_system_process")
#' matr <- pbmc@misc$AUCell$GO$immune_system_process
#' matr <- RenameGO(matr)
#' Heatmap(CalcStats(matr, f = pbmc$cluster, order = "p", n = 5), lab_fill = "zscore")
#'
#' # GSEA using Reactome database
#' pbmc <- GeneSetAnalysisReactome(pbmc, parent = "Immune System")
#' matr <- pbmc@misc$AUCell$Reactome$`Immune System`
#' matr <- RenameReactome(matr)
#' Heatmap(CalcStats(matr, f = pbmc$cluster, order = "p", n = 5), lab_fill = "zscore")
#'
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
