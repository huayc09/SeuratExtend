% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/generics.R, R/WaterfallPlot.R
\name{WaterfallPlot}
\alias{WaterfallPlot}
\alias{WaterfallPlot.Seurat}
\alias{WaterfallPlot.default}
\title{Waterfall Plot}
\usage{
WaterfallPlot(object, ...)

\method{WaterfallPlot}{Seurat}(
  seu,
  features,
  group.by = NULL,
  cell = NULL,
  slot = "data",
  assay = NULL,
  priority = c("expr", "none"),
  ident.1 = NULL,
  ident.2 = NULL,
  exp.transform = (length == "logFC"),
  order = TRUE,
  length = "logFC",
  color = "p",
  len.threshold = 0,
  col.threshold = 0,
  color_theme = c(low = muted("blue"), mid = "white", high = muted("red")),
  top.n = NULL,
  flip = FALSE,
  y.label = NULL,
  angle = NULL,
  hjust = NULL,
  vjust = NULL,
  title = NULL
)

\method{WaterfallPlot}{default}(
  matr,
  f,
  ident.1 = NULL,
  ident.2 = NULL,
  exp.transform = FALSE,
  order = TRUE,
  length = "tscore",
  color = "tscore",
  len.threshold = 0,
  col.threshold = 0,
  color_theme = c(low = muted("blue"), mid = "white", high = muted("red")),
  top.n = NULL,
  flip = TRUE,
  y.label = NULL,
  angle = NULL,
  hjust = NULL,
  vjust = NULL,
  title = NULL
)
}
\arguments{
\item{object}{Either a Seurat object or a matrix.}

\item{...}{Additional arguments passed to other methods.}

\item{seu}{A Seurat object. Only applicable when using the Seurat method.}

\item{features}{Features to be plotted, which can include gene expression, metrics, PC scores, or any other data that can be retrieved using the `FetchData()` function. Only applicable for the Seurat method.}

\item{group.by}{A variable from `meta.data` for grouping, or a character vector of the same length as the number of cells. Only applicable for the Seurat method.}

\item{cell}{Cell identifiers to be used in the plot. Defaults to all cells. Only applicable for the Seurat method.}

\item{slot}{Slot from which to retrieve feature data. Only applicable for the Seurat method.}

\item{assay}{Name of the assay to use. If not specified, the active assay will be used. Only applicable for the Seurat method.}

\item{priority}{If set to "expr", the function will fetch data from the expression matrix rather than `meta.data`. Only applicable for the Seurat method.}

\item{ident.1}{The primary identity class. If not specified, the first class will be used.}

\item{ident.2}{An optional secondary identity class for comparison. If NULL, comparisons will be made against all other cells.}

\item{exp.transform}{Indicates whether to transform data using `expm1`. This is particularly useful when calculating the log fold change of normalized gene counts. Defaults to FALSE for matrix input and TRUE for Seurat object input.}

\item{order}{Determines whether the features should be ordered. Defaults to TRUE.}

\item{length}{Specifies the statistic to determine the length of the bar. Possible values are "tscore" (default for matrix input), "p", or "logFC" (default for Seurat object input).}

\item{color}{Specifies the statistic to determine the color of the bar. Possible values are "tscore" (default), "p" (default for Seurat object input), or "logFC".}

\item{len.threshold}{Excludes features with a value for the `length` parameter below this threshold. Defaults to 0.}

\item{col.threshold}{Excludes features with a value for the `color` parameter below this threshold. Defaults to 0.}

\item{top.n}{Retains only the top `n` bars in both positive and negative directions. If `length(top.n)` is 1, the function retains `top.n` bars for both positive and negative directions. If `length(top.n)` is 2, it retains `top.n[1]` positive bars and `top.n[2]` negative bars. Defaults to NULL.}

\item{flip}{Determines whether the plot should be flipped. Defaults to TRUE for matrix input and FALSE for Seurat object input.}

\item{y.label}{Label for the y-axis. Defaults to "length".}

\item{angle}{Angle of the x-axis labels. This argument is passed to `element_text()`.}

\item{hjust}{Horizontal justification for the x-axis labels. This argument is passed to `element_text()`.}

\item{vjust}{Vertical justification for the x-axis labels. This argument is passed to `element_text()`.}

\item{title}{Title of the plot. Defaults to NULL.}

\item{matr}{A matrix or data frame where rows represent features and columns represent cells.}

\item{f}{A factor or vector indicating the identity of each cell. The length should match the number of columns in `matr`.}

\item{color_scheme}{Specifies the color gradient for the heatmap visualization.
  This parameter accepts multiple input formats to provide flexibility in defining color schemes:

  - Predefined color schemes: Users can specify "A", "B", "C", "D", or "E" to use color schemes from the `viridis` package.

  - Named vector for three-point gradients: Provide a named vector with keys "low", "mid", and "high" to define colors at these specific data points. The "mid" value is typically centered at zero, allowing for a diverging color scheme.
    Example: `c(low = "blue", mid = "white", high = "red")`

  - Two-point gradient: Provide a named vector with keys "low" and "high" to create a simple two-color gradient. Example: `c(low = "blue", high = "red")`

  - Custom color gradient: Users can provide a vector of colors to generate a custom gradient across multiple values. This is suitable for more complex data ranges and visual preferences.}
}
\value{
A ggplot object.
}
\description{
Generates a waterfall plot.
}
\details{
For more detailed usage, see the examples provided.
}
\examples{
# First, create a matrix using the GeneSetAnalysis() function.
# Rows represent the Hallmark 50 genesets, and columns represent cells.
pbmc <- GeneSetAnalysis(pbmc, genesets = hall50$human)
matr <- pbmc@misc$AUCell$genesets

# Generate a waterfall plot comparing CD14+ Mono with CD8 T cells.
WaterfallPlot(matr, f = pbmc$cluster, ident.1 = "Mono CD14", ident.2 = "CD8 T cell")

# Keep only bars with a length (tscore, in this instance) greater than 1.
WaterfallPlot(matr, f = pbmc$cluster, ident.1 = "Mono CD14", ident.2 = "CD8 T cell", len.threshold = 1)

# Use a Seurat object as input and compare 100 genes, using LogFC as bar length.
genes <- VariableFeatures(pbmc)[1:80]
WaterfallPlot(
  pbmc, group.by = "cluster", features = genes,
  ident.1 = "Mono CD14", ident.2 = "CD8 T cell", length = "logFC")

# Keep only the top and bottom 20 genes.
WaterfallPlot(
  pbmc, group.by = "cluster", features = genes,
  ident.1 = "Mono CD14", ident.2 = "CD8 T cell", length = "logFC",
  top.n = 20)
}
