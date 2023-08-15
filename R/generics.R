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

#' @title Calculate Matrix Statistics Grouped by Clusters
#' @description Computes various statistics (mean, median, zscore, tscore, etc.) on a matrix, grouped by clusters.
#' @param object An object. Can either be a Seurat object or a matrix.
#' @param ... Arguments passed to other methods
#' @return A matrix.
#' @details Computes statistics for each feature across cell types. "p" value is determined by t-test.
#' For log-normalized data's LogFC computation, it's advised to set `exp.transform` to TRUE.
#' @examples
#' library(Seurat)
#' library(SeuratExtend)
#'
#' # Using a Seurat object as input. First, select some genes for calculation:
#' genes <- VariableFeatures(pbmc)[1:20]
#'
#' # Calculate zscore, grouping by the default 'ident' (cluster):
#' genes.zscore <- CalcStats(pbmc, features = genes, method = "zscore")
#' head(genes.zscore)
#'
#' # Visualize with a heatmap:
#' Heatmap(genes.zscore, lab_fill = "zscore")
#'
#' # Select more genes and retain the top 4 genes of each cluster, sorted by p-value.
#' # This can be a convenient method to display the top marker genes of each cluster:
#' genes <- VariableFeatures(pbmc)
#' genes.zscore <- CalcStats(
#'   pbmc, features = genes, method = "zscore", group.by = "cluster",
#'   order = "p", n = 4)
#' Heatmap(genes.zscore, lab_fill = "zscore")
#'
#' # It's also possible to use a matrix as input. For instance, we initially perform
#' # Geneset Enrichment Analysis (GSEA) using the Hallmark 50 geneset and obtain the AUCell matrix
#' # (rows represent pathways, columns represent cells):
#' pbmc <- GeneSetAnalysis(pbmc, genesets = hall50$human)
#' matr <- pbmc@misc$AUCell$genesets
#'
#' # Next, we calculate the zscore, grouped by 'cluster':
#' gsea.zscore <- CalcStats(matr, f = pbmc$cluster, method = "zscore")
#' Heatmap(gsea.zscore, lab_fill = "zscore")
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
