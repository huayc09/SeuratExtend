% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/generics.R, R/CalcStats.R
\name{CalcStats}
\alias{CalcStats}
\alias{CalcStats.Seurat}
\alias{CalcStats.default}
\title{Calculate statistics of matrix}
\usage{
CalcStats(object, ...)

\method{CalcStats}{Seurat}(
  seu,
  features,
  group.by = NULL,
  cells = NULL,
  slot = "data",
  method = "zscore",
  exp.transform = F,
  order = NULL,
  n = Inf,
  p.threshold = 0.05
)

\method{CalcStats}{default}(
  matr,
  f,
  method = "zscore",
  exp.transform = F,
  t = F,
  order = NULL,
  n = Inf,
  p.threshold = 0.05
)
}
\arguments{
\item{object}{An object; either Seurat object or matrix}

\item{...}{Arguments passed to other methods}

\item{seu}{(Seurat version) Seurat object}

\item{features}{(Seurat version) Features to plot (gene expression, metrics, PC scores,
anything that can be retreived by FetchData), Default: NULL (All features
in matrix)}

\item{group.by}{(Seurat version) A variable name in meta.data to
group the violin plots by}

\item{cells}{(Seurat version) Cell names to use, Default: all cells}

\item{slot}{Slot to pull feature data for}

\item{method}{Should be either "mean", "median", "zscore", "tscore",
"p", or "logFC". Default: 'zscore'}

\item{exp.transform}{Whether to transform the data with
\code{\link[base:log]{expm1}}, Default: F}

\item{order}{Re-order rows by "value" or "p" (in t.test)}

\item{n}{Top n rows of each cluster. Ignored when
\code{order} is \code{NULL}. Default: \code{Inf}}

\item{p.threshold}{When \code{order = "p"}, rows with pvalue greater
than \code{p.threshold} are removed.}

\item{matr}{Matrix or data frame.Row - features; columns - cells}

\item{f}{Factor or vector. Identity of each cell. Should be the
same length of cells}

\item{t}{If your matrix has features in columns and cells in rows,
you should transpose the matrix first. Default: F}
}
\value{
matrix
}
\description{
Calculate mean, median, tscore or zscore of matrix, group by clusters.
}
\details{
Calculate statistis of each feature of each cell type. "p" value according
to t.test.

For calculating LogFC of log normalized data, it is recommended to set
\code{exp.transform} to \code{TRUE}
}
\examples{
# Use Seurat object as input

genes <- VariableFeatures(pbmc)[1:10]

CalcStats(pbmc, genes, method = "zscore")

CalcStats(pbmc, genes, method = "tscore", group.by = "orig.ident")

# Use matrix as input

matr <- FetchData(pbmc, genes)
matr

CalcStats(matr, f = pbmc$cluster, method = "zscore", t = TRUE)

# Heatmap
CalcStats(pbmc, genes) \%>\% Heatmap(lab_fill = "zscore")

# Order rows
CalcStats(pbmc, VariableFeatures(pbmc), method = "zscore", order = "p", n = 4) \%>\%
  Heatmap(lab_fill = "zscore")

}
