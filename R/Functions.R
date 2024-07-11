#' @title Proportion of features
#' @description Check the Proportion of positive cells (default: expression above 0)
#' in certain clusters
#' @param seu Seurat object
#' @param feature Features to plot (gene expression, metrics, PC scores, anything that can be
#' retreived by FetchData)
#' @param ident cluster name, Default: all clusters
#' @param group.by A variable name in meta.data to group by, if you don't want to use
#' default Idents of Seurat object
#' @param DefaultAssay Assay to use, Default: 'RNA'
#' @param above Above which value will be considered as positive cell, Default: 0
#' @param total If to calculate proportion in total cells of selected clusters, Default: F
#' @param if.expressed If only return logical value, Default: F
#' @param min.pct Minimum fraction of min.pct cells in the cluster will be considered "expressed",
#' only use this parameter when "if.expressed" is set to TRUE. Default: 0.1
#' @return A matrix
#' @details None.
#' @examples
#' genes <- VariableFeatures(pbmc)[1:10]
#'
#' feature_percent(pbmc, feature = genes)
#'
#' # Count one cell as positive only the expression is above 2
#' feature_percent(pbmc, feature = genes, above = 2)
#'
#' # Only check subset of clusters
#' feature_percent(pbmc, feature = genes, ident = c("B cell", "CD8 T cell"))
#'
#' # Group by a different variable
#' feature_percent(pbmc, feature = genes, group.by = "orig.ident")
#'
#' # Also check the proportion of expressed cells in total clusters
#' feature_percent(pbmc, feature = genes, total = T)
#'
#' # only show logical value; 'TRUE' if more than 20% cells are positive
#' feature_percent(pbmc, feature = genes, if.expressed = T, min.pct = 0.2)
#'
#' @rdname feature_percent
#' @export

feature_percent <- function(
  seu,
  feature,
  ident = NULL,
  group.by = NULL,
  DefaultAssay = "RNA",
  above = 0,
  total = F,
  if.expressed = F,
  min.pct = 0.1
) {
  library(Seurat)
  library(rlist)
  DefaultAssay(seu) <- DefaultAssay
  if(is.null(group.by)){
    f <- Idents(seu)
  } else {
    f <- factor(seu@meta.data[,group.by])
  }
  ident <- ident %||% levels(f)
  cells <- colnames(seu)[f %in% ident]
  m <- FetchData(seu, vars = feature, cells = cells)
  f2 <- f[f %in% ident]
  f2 <- factor(f2)
  m2 <- split(m, f2)
  if(total) m2[["total"]] <- m
  m2 <- lapply(m2, function(x){
    apply(x, 2, function(y) sum(y > above)) / nrow(x)
  })
  m2 <- list.cbind(m2)
  if(if.expressed) m2 <- m2 > min.pct
  return(m2)
}



