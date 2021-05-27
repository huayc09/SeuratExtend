#' @title Modified Violin plot
#' @description Modified violin plots of Seurat style. More compact for multiple variables.
#' Input can be Seurat object or data frame.
#' @param ... Arguments passed to other methods
#' @return ggplot object
#' @details See above
#' @examples
#' VlnPlot2(pbmc, features = c("CD3D","CD79A"))
#'
#' VlnPlot2(pbmc, features = c("CD3D","CD79A"), split.by = "orig.ident")
#'
#' VlnPlot2(pbmc, features = c("CD3D","CD79A"), violin = F)
#'
#' @rdname VlnPlot2
#' @export

VlnPlot2 <- function(object, ...) {
  UseMethod(generic = "VlnPlot2", object = object)
}
