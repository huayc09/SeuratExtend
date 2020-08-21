#' @title FUNCTION_TITLE
#' @description FUNCTION_DESCRIPTION
#' @param loom.path PARAM_DESCRIPTION
#' @param seu PARAM_DESCRIPTION, Default: NULL
#' @return OUTPUT_DESCRIPTION
#' @details DETAILS
#' @examples
#' \dontrun{
#' if(interactive()){
#'  #EXAMPLE1
#'  }
#' }
#' @rdname ImportPyscenicLoom
#' @export

ImportPyscenicLoom <- function(loom.path, seu = NULL) {
  library(Seurat)
  import("loomR")

  lfile <- connect(loom.path, skip.validate = T)
  if(is.null(seu)) {} else if(length(lfile$col.attrs$CellID[]) != length(colnames(seu))) {
    stop("Loom file and Seurat object have different cell numbers")
  }else if(!identical(lfile$col.attrs$CellID[], colnames(seu))){
    warning("Loom file and Seurat object have the same cell numbers but different cell names")
  }

  RegulonsAUC <- lfile$col.attrs$RegulonsAUC[]
  colnames(RegulonsAUC) <- sub("_(+)", "", colnames(RegulonsAUC), fixed = T)
  rownames(RegulonsAUC) <- lfile$col.attrs$CellID[]

  Regulons <- lfile$row.attrs$Regulons[]
  colnames(Regulons) <- sub("_(+)", "", colnames(Regulons), fixed = T)
  genes <- lfile$row.attrs$Gene
  Regulons <- apply(Regulons, 2, function(x) genes[which(x==1)])

  lfile$close_all()
  results <- list(RegulonsAUC = RegulonsAUC,
                  Regulons = Regulons)

  if(is.null(seu)) {
    return(results)
  } else {
    seu@misc[["SCENIC"]] <- results
    return(seu)
  }
}
