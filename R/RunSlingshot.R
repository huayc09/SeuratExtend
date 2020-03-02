#' @title FUNCTION_TITLE
#' @description FUNCTION_DESCRIPTION
#' @param Seu PARAM_DESCRIPTION
#' @param group.by PARAM_DESCRIPTION
#' @param reducedDim PARAM_DESCRIPTION, Default: 'DM'
#' @param assay PARAM_DESCRIPTION, Default: NULL
#' @param start.clus PARAM_DESCRIPTION, Default: NULL
#' @return OUTPUT_DESCRIPTION
#' @details DETAILS
#' @examples
#' \dontrun{
#' if(interactive()){
#'  #EXAMPLE1
#'  }
#' }
#' @rdname RunSlingshot
#' @export
RunSlingshot <- function(Seu, group.by, reducedDim = "DM", assay = NULL,
                         start.clus = NULL){
  library(slingshot)
  library(Seurat)
  library(rlang)
  assay <- assay %||% DefaultAssay(Seu)
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
  Seu@meta.data <- cbind(Seu@meta.data, sce@colData[slingPseudotimeName]) %>% as.data.frame()
  return(Seu)
}
# Seu <- RunSlingshot(Seu, group.by = "cluster", reducedDim = "PCA")
