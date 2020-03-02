#' @title FUNCTION_TITLE
#' @description FUNCTION_DESCRIPTION
#' @param Seu PARAM_DESCRIPTION
#' @param features PARAM_DESCRIPTION
#' @param group.by PARAM_DESCRIPTION
#' @param method PARAM_DESCRIPTION
#' @param assay PARAM_DESCRIPTION, Default: 'RNA'
#' @return OUTPUT_DESCRIPTION
#' @details DETAILS
#' @examples
#' \dontrun{
#' if(interactive()){
#'  #EXAMPLE1
#'  }
#' }
#' @rdname CalcScoreGeneral_v3
#' @export
CalcScoreGeneral_v3<-function(Seu, features, group.by, method, assay = "RNA"){
  library(Seurat)
  if(any(!features %in% rownames(Seu))) {
    WrongGenes <- setdiff(features, rownames(Seu))
    warning(paste0(paste(WrongGenes, collapse = ", "), " are not detected in the matrix"))
    features <- intersect(features, rownames(Seu))
  }
  matr <- as.matrix(GetAssayData(Seu, assay = assay)[features,])
  f <- factor(Seu@meta.data[, group.by])
  return(CalcScoreGeneral(matr, f, method))
}
