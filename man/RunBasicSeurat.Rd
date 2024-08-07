% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/RunBasicSeurat.R
\name{RunBasicSeurat}
\alias{RunBasicSeurat}
\title{Run Standard Seurat Pipeline}
\usage{
RunBasicSeurat(
  seu,
  spe = NULL,
  nFeature_RNA.min = 500,
  nFeature_RNA.max = Inf,
  percent.mt.max = 20,
  dims = 1:10,
  resolution = 0.5,
  reduction = NULL,
  harmony.by = NULL,
  force.Normalize = FALSE,
  force.RunPCA = FALSE,
  vars.to.regress = NULL,
  force.RunHarmony = FALSE,
  force.FindNeighbors = FALSE,
  force.FindClusters = FALSE,
  force.RunUMAP = FALSE
)
}
\arguments{
\item{seu}{A Seurat object that will be processed.}

\item{spe}{Either "human" or "mouse", used in calculating percent.mt. Default: NULL.}

\item{nFeature_RNA.min}{The minimum number of RNA features to include in the subset. Default: 500.}

\item{nFeature_RNA.max}{The maximum number of RNA features to include in the subset. Default: Inf.}

\item{percent.mt.max}{The maximum percentage of mitochondrial content allowed for cells in the subset. Default: 20.}

\item{dims}{The dimensions to be used in RunPCA and FindNeighbors. Default: 1:10.}

\item{resolution}{The resolution parameter for the clustering algorithm. Default: 0.5.}

\item{reduction}{The dimensional reduction technique to be used ('pca' or 'harmony'). If NULL, the method is determined based on available data in the Seurat object. Default: NULL.}

\item{harmony.by}{A string representing a metadata column in the Seurat object that categorizes data, which might need transformation for numeric-only issues before running Harmony. Default: NULL.}

\item{force.Normalize}{Forces the normalization step regardless of previous command history. Default: FALSE.}

\item{force.RunPCA}{Forces the PCA step regardless of previous command history. Default: FALSE.}

\item{vars.to.regress}{Variables to regress out during scaling. Default: NULL.}

\item{force.RunHarmony}{Forces the Harmony integration step regardless of previous command history or detection. Default: FALSE.}

\item{force.FindNeighbors}{Forces the neighbor finding step regardless of previous command history. Default: FALSE.}

\item{force.FindClusters}{Forces the clustering step regardless of previous command history. Default: FALSE.}

\item{force.RunUMAP}{Forces the UMAP step regardless of previous command history. Default: FALSE.}
}
\value{
Returns the modified Seurat object after performing the requested operations and updates.
}
\description{
This function processes a Seurat object through various steps including normalization, PCA, and clustering based on specified parameters. It allows for conditional execution of each step based on prior executions and parameter changes.
}
\details{
This function provides flexible and conditional execution of various data processing steps in a Seurat pipeline. Depending on the parameters set, the function intelligently determines whether to rerun specific steps based on changes in parameters or previous runs recorded in the Seurat object's command history. This allows for efficient reuse of existing results while ensuring updates are applied when necessary.
}
\examples{
library(SeuratExtend)
pbmc <- RunBasicSeurat(pbmc, force.Normalize = TRUE)
DimPlot2(pbmc, group.by = "cluster")
}
