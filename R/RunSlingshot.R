#' @title Calculate Pseudotime and Map Trajectories Using Slingshot
#' @description This function integrates Slingshot for pseudotime analysis directly within a Seurat workflow, enabling the mapping of cellular trajectories based on user-defined cluster assignments and starting clusters.
#' @param Seu A Seurat object containing single-cell RNA-seq data with precomputed clusters and necessary dimensional reductions.
#' @param group.by Specifies the metadata column in the Seurat object used to define groups or clusters for trajectory analysis.
#' @param reducedDim The name of the dimensionality reduction to use for trajectory inference. Default: 'PCA'.
#' @param assay Name of the assay to employ. Defaults to the active assay.
#' @param start.clus Optional; specifies the starting cluster to initiate trajectory inference, which can be crucial for directing the trajectory analysis in a biologically meaningful way. Default: NULL.
#' @return Modifies the Seurat object by adding Slingshot-derived pseudotime and trajectory data to the `@misc$slingshot` slot.
#' @details Slingshot is a flexible tool for trajectory analysis that uses cluster-based minimum spanning trees to infer developmental pathways. This function facilitates the integration of Slingshot with Seurat objects, allowing for pseudotime calculations that are sensitive to the underlying data structure and cluster dynamics. It is especially useful for complex datasets where multiple trajectories might exist, helping to uncover hidden patterns in cellular development.
#' @examples
#' library(Seurat)
#' library(SeuratExtend)
#'
#' # Load an example Seurat Object
#' mye_small <- readRDS(url("https://zenodo.org/records/10944066/files/pbmc10k_mye_small_velocyto.rds", "rb"))
#'
#' # Run Slingshot
#' mye_small <- RunSlingshot(mye_small, group.by = "cluster", start.clus = "Mono CD14")
#'
#' # Access and visualize the Slingshot output
#' sling <- mye_small@misc$slingshot$PCA$SlingPseudotime
#' print(head(sling))
#' mye_small@meta.data[,colnames(sling)] <- as.data.frame(sling)
#' DimPlot2(mye_small, features = colnames(sling), cols = "C")
#'
#' @rdname RunSlingshot
#' @export


RunSlingshot <- function(
    Seu, group.by,
    reducedDim = "PCA",
    assay = NULL,
    start.clus = NULL
){
  import("slingshot")
  library(Seurat)
  library(dplyr)
  assay <- assay %||% DefaultAssay(Seu)
  cols_to_remove <- grep(pattern = "slingPseudotime_\\d+", colnames(Seu@meta.data), value = TRUE)
  Seu@meta.data <- Seu@meta.data[, !colnames(Seu@meta.data) %in% cols_to_remove]

  sce <- as.SingleCellExperiment(Seu, assay = assay)
  Time1 <- Sys.time()
  sce <- slingshot(sce, clusterLabels = group.by, reducedDim = reducedDim, start.clus = start.clus)
  Time2 <- Sys.time()
  message <- paste0("Time usage of slingshot: ", round(difftime(Time2, Time1, units='mins'), digits = 2), " mins")
  message(message)
  Seu@misc[["slingshot"]][[reducedDim]] <- list()
  Seu@misc[["slingshot"]][[reducedDim]][["SlingshotDataSet"]] <- SlingshotDataSet(sce)
  Seu@misc[["slingshot"]][[reducedDim]][["f"]] <- Seu@meta.data[,group.by]
  slingPseudotimeName <- colnames(sce@colData) %>% .[grepl("slingPseudotime", .)]
  Seu@misc[["slingshot"]][[reducedDim]][["SlingPseudotime"]] <- sce@colData[slingPseudotimeName]
  return(Seu)
}

