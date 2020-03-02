#' @title FUNCTION_TITLE
#' @description FUNCTION_DESCRIPTION
#' @param seu PARAM_DESCRIPTION
#' @param slot PARAM_DESCRIPTION, Default: 'counts'
#' @param assay PARAM_DESCRIPTION, Default: 'RNA'
#' @return OUTPUT_DESCRIPTION
#' @details DETAILS
#' @examples
#' \dontrun{
#' if(interactive()){
#'  #EXAMPLE1
#'  }
#' }
#' @rdname BuildAUCRank
#' @export
BuildAUCRank <- function(seu, slot = "counts", assay = "RNA"){
  library(AUCell)
  library(Seurat)
  matr <- GetAssayData(seu, slot = slot, assay = assay)
  seu@misc$AUCell<-list()
  seu@misc$AUCell[["cells_rankings"]] <- AUCell_buildRankings(matr, plotStats=TRUE)
  return(seu)
}
# seu <- BuildAUCRank(seu)
