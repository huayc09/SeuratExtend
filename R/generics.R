#' @title Modified violin plot
#' @description Modified violin plots of Seurat style. More compact for multiple variables.
#' Input can be Seurat object or data frame.
#' @param object An object; either Seurat object or matrix
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

#' @title Calculate statistics of matrix
#' @description Calculate mean, median, tscore or zscore of matrix, group by clusters.
#' @param object An object; either Seurat object or matrix
#' @param ... Arguments passed to other methods
#' @return matrix
#' @details Calculate statistis of each feature of each cell type. "p" value according
#' to t.test.
#'
#' For calculating LogFC of log normalized data, it is recommended to set
#' \code{exp.transform} to \code{TRUE}
#' @examples
#' # Use Seurat object as input
#'
#' genes <- VariableFeatures(pbmc)[1:10]
#'
#' CalcStats(pbmc, genes, method = "zscore")
#'
#' CalcStats(pbmc, genes, method = "tscore", group.by = "orig.ident")
#'
#' # Use matrix as input
#'
#' matr <- FetchData(pbmc, genes)
#' matr
#'
#' CalcStats(matr, f = pbmc$cluster, method = "zscore", t = TRUE)
#'
#' # Heatmap
#' CalcStats(pbmc, genes) %>% Heatmap(lab_fill = "zscore")
#'
#' # Order rows
#' CalcStats(pbmc, VariableFeatures(pbmc), method = "zscore", order = "p", n = 4) %>%
#'   Heatmap(lab_fill = "zscore")
#'
#' @rdname CalcStats
#' @export

CalcStats <- function(object, ...) {
  UseMethod(generic = "CalcStats", object = object)
}

#' @title Waterfall Plot
#' @description Generates a waterfall plot.
#' @param object Either a Seurat object or a matrix.
#' @param ... Additional arguments passed to other methods.
#' @return A ggplot object.
#' @details For more detailed usage, see the examples provided.
#' @examples
#' # First, create a matrix using the GeneSetAnalysis() function.
#' # Rows represent the Hallmark 50 genesets, and columns represent cells.
#' pbmc <- GeneSetAnalysis(pbmc, genesets = hall50$human)
#' matr <- pbmc@misc$AUCell$genesets
#'
#' # Generate a waterfall plot comparing CD14+ Mono with CD8 T cells.
#' WaterfallPlot(matr, f = pbmc$cluster, ident.1 = "CD14+ Mono", ident.2 = "CD8 T")
#'
#' # Keep only bars with a length (tscore, in this instance) greater than 1.
#' WaterfallPlot(matr, f = pbmc$cluster, ident.1 = "CD14+ Mono", ident.2 = "CD8 T", len.threshold = 1)
#'
#' # Use a Seurat object as input and compare 100 genes, using LogFC as bar length.
#' genes <- VariableFeatures(pbmc)[1:100]
#' WaterfallPlot(
#'   pbmc, group.by = "cluster", features = genes,
#'   ident.1 = "CD14+ Mono", ident.2 = "CD8 T", length = "logFC")
#'
#' # Keep only the top and bottom 20 genes.
#' WaterfallPlot(
#'   pbmc, group.by = "cluster", features = genes,
#'   ident.1 = "CD14+ Mono", ident.2 = "CD8 T", length = "logFC",
#'   top.n = 20)
#' @rdname WaterfallPlot
#' @export

WaterfallPlot <- function(object, ...) {
  UseMethod(generic = "WaterfallPlot", object = object)
}
