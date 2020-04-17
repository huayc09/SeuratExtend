#' @title FUNCTION_TITLE
#' @description FUNCTION_DESCRIPTION
#' @param Seu PARAM_DESCRIPTION
#' @param assay PARAM_DESCRIPTION, Default: NULL
#' @return OUTPUT_DESCRIPTION
#' @details DETAILS
#' @examples
#' \dontrun{
#' if(interactive()){
#'  #EXAMPLE1
#'  }
#' }
#' @rdname RunDiffusionMap
#' @export
RunDiffusionMap <- function(Seu, assay = NULL){
  library(Seurat)
  library(rlang)
  library(destiny)
  library(magrittr)
  DefaultAssay(Seu) <- assay %||% DefaultAssay(Seu)
  dm <- DiffusionMap(t(as.matrix(GetAssayData(Seu))))
  Seu[["dm"]] <- CreateDimReducObject(embeddings = dm@eigenvectors %>% set_rownames(colnames(Seu)),
                                      key = "DM_", assay = DefaultAssay(Seu))
  Seu@misc[["DiffusionMap"]] <- dm
  return(Seu)
}
