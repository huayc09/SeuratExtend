#' @title Enhanced Violin Plot
#' @description Generates advanced violin plots distinct from Seurat's VlnPlot. This improved version offers a more compact design for efficient space utilization, the ability to overlay a boxplot, and convenient inclusion of statistical annotations. The function accommodates input in the form of either a Seurat object or a data frame.
#' @param object An object, either a Seurat object or matrix.
#' @return ggplot object
#' @return A ggplot object.
#' @details This function provides a range of customization options to generate violin plots, with or without additional graphical elements such as boxplots and points. You can specify features to plot, control the appearance of violin plots, boxplots, and points, adjust point position, group data, split violin plots, selectively use cells, control the column layout of multiple plots, and add statistical annotations. For a detailed overview, refer to the provided examples.
#'
#' @examples
#' library(Seurat)
#' library(SeuratExtend)
#'
#' # Using Seurat object as input:
#'
#' # Basic violin plot with box plot and points:
#' genes <- c("CD3D","CD14","CD79A")
#' VlnPlot2(pbmc, features = genes, ncol = 1)
#'
#' # Without violin plot, only box plot, and use quasirandom style for point position adjustment:
#' VlnPlot2(pbmc, features = genes, violin = F, pt.style = "quasirandom", ncol = 1)
#'
#' # Hide points but display outlier points of the box plot:
#' VlnPlot2(pbmc, features = genes, pt = FALSE, ncol = 1)
#'
#' # When hiding points, outliers are shown by default (recommended). However, outliers can be hidden for a cleaner appearance:
#' VlnPlot2(pbmc, features = genes, pt = FALSE, hide.outlier = T, ncol = 1)
#'
#' # Group by cluster and split each cluster by samples:
#' VlnPlot2(pbmc, features = genes, group.by = "cluster", split.by = "orig.ident")
#'
#' # Display only cells from certain subtypes (e.g. B cell, CD14+ Mono, and CD8 T cell), and arrange plots in 3 columns:
#' cells <- colnames(pbmc)[pbmc$cluster %in% c("B cell", "Mono CD14", "CD8 T cell")]
#' VlnPlot2(pbmc, features = genes, group.by = "cluster", cells = cells)
#'
#' # Add statistical annotations using the Wilcoxon test for pairwise comparisons and hide non-significant results:
#' VlnPlot2(pbmc, features = genes, group.by = "cluster", cell = cells,
#'          stat.method = "wilcox.test", hide.ns = TRUE)
#'
#' # Display statistics comparing only the specified clusters using a t-test:
#' VlnPlot2(pbmc, features = genes, group.by = "cluster", cell = cells,
#'          stat.method = "t.test", comparisons = list(c(1,2), c(1,3)), hide.ns = FALSE)
#'
#' # Using a matrix as input:
#' # For instance, after performing Geneset Enrichment Analysis (GSEA) using the Hallmark 50 geneset to obtain the AUCell matrix:
#' pbmc <- GeneSetAnalysis(pbmc, genesets = hall50$human)
#' matr <- pbmc@misc$AUCell$genesets
#'
#' # Create violin plots for the first three pathways:
#' VlnPlot2(matr[1:3,], f = pbmc$cluster, ncol = 1)
#' @rdname VlnPlot2
#' @export

VlnPlot2 <- function(object, ...) {
  UseMethod(generic = "VlnPlot2", object = object)
}

#' @title Calculate Matrix Statistics Grouped by Clusters
#' @description Computes various statistics (mean, median, zscore, tscore, etc.) on a matrix, grouped by clusters.
#' @param object An object. Can either be a Seurat object or a matrix.
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
#' WaterfallPlot(matr, f = pbmc$cluster, ident.1 = "Mono CD14", ident.2 = "CD8 T cell")
#'
#' # Keep only bars with a length (tscore, in this instance) greater than 1.
#' WaterfallPlot(matr, f = pbmc$cluster, ident.1 = "Mono CD14", ident.2 = "CD8 T cell", len.threshold = 1)
#'
#' # Use a Seurat object as input and compare 100 genes, using LogFC as bar length.
#' genes <- VariableFeatures(pbmc)[1:80]
#' WaterfallPlot(
#'   pbmc, group.by = "cluster", features = genes,
#'   ident.1 = "Mono CD14", ident.2 = "CD8 T cell", length = "logFC")
#'
#' # Keep only the top and bottom 20 genes.
#' WaterfallPlot(
#'   pbmc, group.by = "cluster", features = genes,
#'   ident.1 = "Mono CD14", ident.2 = "CD8 T cell", length = "logFC",
#'   top.n = 20)
#' @rdname WaterfallPlot
#' @export

WaterfallPlot <- function(object, ...) {
  UseMethod(generic = "WaterfallPlot", object = object)
}
