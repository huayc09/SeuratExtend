% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/scVelo.R
\name{scVelo.SeuratToLoom}
\alias{scVelo.SeuratToLoom}
\alias{scVelo.RunBasic}
\title{(Deprecated) Export Seurat Object and Velocyto Data to Loom for scVelo Analysis}
\usage{
scVelo.SeuratToLoom(
  seu,
  export.name,
  velocyto.loom.path = NULL,
  velocyto.loom.filenames,
  fxn.convert.loomtoseurat.cellname = NULL,
  cell.id.match.table = NULL,
  prefix = NULL,
  postfix = NULL
)

scVelo.RunBasic(loom, save.adata = "adata.obj")
}
\arguments{
\item{seu}{Seurat object containing single-cell RNA sequencing data.}

\item{export.name}{Name of the resulting loom file to be created.}

\item{velocyto.loom.path}{Directory containing the Velocyto loom files; defaults to the current working directory. Default: getwd()}

\item{velocyto.loom.filenames}{Filenames of Velocyto output (loom files) that contain spliced and unspliced RNA data.}

\item{fxn.convert.loomtoseurat.cellname}{(Experimental) Function to convert cell names in Velocyto loom files to Seurat object cell names, if they differ.}

\item{prefix}{Optional prefix added to Seurat cell IDs, typically used to include sample names. Default: NULL}

\item{postfix}{Optional postfix added to Seurat cell IDs, typically used to denote unique identifiers such as '-1'. Default: NULL}

\item{loom}{Path to the loom file created by `scVelo.SeuratToLoom`.}

\item{save.adata}{Name and path where the AnnData object should be saved, typically ending in '.h5ad'. Default: 'adata.obj'}
}
\value{
Generates a loom file that integrates Seurat and Velocyto data, stored in the specified export path.
}
\description{
This function converts a Seurat object and associated Velocyto loom files into a single loom file that integrates both spliced and unspliced RNA matrices. It is part of an older workflow and is soft deprecated in favor of `scVelo.SeuratToAnndata`. However, it is preserved for reproducibility of earlier work that cited its use.
}
\details{
This set of functions supports the transition from Seurat objects to an integrated loom file via `scVelo.SeuratToLoom`, and from loom files to AnnData objects via `scVelo.RunBasic`, ready for basic to advanced scVelo analysis in Python. This older method is preserved for compatibility with previous studies and publications, ensuring reproducibility of earlier work. This older method is preserved for compatibility with previous studies and publications, but users are encouraged to transition to `scVelo.SeuratToAnndata` for a more streamlined and updated approach that directly integrates Seurat objects with scVelo analysis in Python's Scanpy library.
}
\seealso{
\code{\link[SeuratExtend:scVelo.SeuratToAnndata]{scVelo.SeuratToAnndata}} for the recommended method of preparing data for scVelo analysis.
}
