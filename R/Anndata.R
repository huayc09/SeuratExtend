#' @title Seurat to AnnData Conversion and Management Tools
#' @description This documentation encompasses a suite of functions designed to facilitate the conversion and management of data between Seurat objects and AnnData structures. The functions cover the entire workflow from initial conversion of Seurat objects to Loom or AnnData formats, adding supplementary data like dimension reductions and metadata, and saving or loading these objects for further analysis.
#' @param seu The Seurat object from which data is being converted or managed.
#' @param add.normdata Boolean indicating whether to also include the normalized gene expression matrix in the loom file.
#' @param add.metadata Boolean indicating whether to also include the meta.data in the loom file.
#' @param adata.path Path to the AnnData (.h5ad) file where the data is to be saved or from which data is to be loaded.
#' @param save.adata Optional; specifies the file path where the AnnData object should be saved after processing. The file will be saved in the .h5ad format. This parameter facilitates easy storage of processed data.
#' @param load.adata Optional; path to a previously saved AnnData object to be loaded. Default: NULL. Useful for resuming analyses or integrating additional data without repeating preprocessing steps.
#' @param filename Filename for the loom file to be created or accessed. It specifies where the loom file will be stored or is stored, which is essential for data conversion tasks.
#' @param layers Optional; a list of additional matrices to include as layers in the loom file, such as velocyto's spliced/unspliced matrices. Each matrix must match the transposed dimensions of the Seurat object (row for cells and columns for genes). Default: NULL.
#' @param loompath Path to the loom file that needs to be loaded into an AnnData object.
#' @param dr Specifies which dimension reductions from the Seurat object should be transferred to the AnnData object. If NULL, all available reductions are included. Default: NULL.
#' @param prefix A prefix to add to the dimension reduction names when importing them into the AnnData object's 'obsm' attribute. This helps in differentiating various data types and maintaining order within the dataset. Default: 'X_'.
#' @param postfix A postfix to add to the dimension reduction names when importing them into the AnnData object's 'obsm' attribute.
#' @param scv.graph Optional; a boolean indicating whether to calculate the velocity graph after importing dimension reduction data. This is useful for kinetic analyses in single-cell studies. Default: FALSE.
#' @param col The name of the metadata column in the Seurat object that needs to be transferred to the AnnData object.
#' @param conda_env Name of the Conda environment where the necessary Python libraries are installed. This environment is used to run Python code from R, bridging Seurat and AnnData functionalities. Default: 'seuratextend'.
#' @return Depending on the function called, this suite returns either a modified Seurat object or an AnnData object, or it may perform save/load operations without returning an object. Specifically:
#'
#' - `Seu2Adata` and `Seu2Loom` convert Seurat objects to AnnData and Loom formats, respectively.
#'
#' - `adata.LoadLoom` loads a Loom file into an AnnData object.
#'
#' - `adata.AddDR` and `adata.AddMetadata` enrich an AnnData object with additional dimension reduction data and metadata from a Seurat object.
#'
#' - `adata.Save` and `adata.Load` handle saving to and loading from disk operations for AnnData objects.
#' @details These functions are designed to ensure seamless interoperability between R, used primarily for Seurat, and Python, used for AnnData via the reticulate package:
#'
#' - `Seu2Adata` directly converts a Seurat object into an AnnData object, optionally adding normalized data and saving the result.
#'
#' - `Seu2Loom` converts a Seurat object into a Loom file, allowing for intermediate storage or further conversion.
#'
#' - `adata.LoadLoom` initiates an AnnData object from a Loom file, useful for integrating with Python-based analyses.
#'
#' - `adata.AddDR` appends dimension reduction results from Seurat to an AnnData object, essential for preserving comprehensive analytical context.
#'
#' - `adata.AddMetadata` transfers additional metadata from Seurat to AnnData, facilitating enriched data analyses.
#'
#' - `adata.Save` provides functionality to save AnnData objects to disk, securing data for future use.
#'
#' - `adata.Load` retrieves AnnData objects from disk, ensuring continuity of analysis across sessions.
#'
#' These tools are integral for researchers employing hybrid workflows that leverage the strengths of both Seurat and Scanpy/AnnData environments, supporting a wide range of computational tasks from data preprocessing to advanced analyses.
#' @examples
#' library(Seurat)
#' library(SeuratExtend)
#' library(reticulate)
#'
#' # Convert Seurat object to loom file and save it locally
#' pbmc_loom_path <- file.path(tempdir(), "pbmc3k_small.loom")
#' Seu2Loom(pbmc, filename = pbmc_loom_path, add.normdata = TRUE)
#'
#' # Load the loom file into an AnnData object
#' adata.LoadLoom(loompath = pbmc_loom_path)
#'
#' # Use reticulate to execute Python code and print the AnnData object details
#' py_run_string("print(adata)")
#'
#' # Add dimension reduction data from Seurat to AnnData
#' adata.AddDR(pbmc)
#' py_run_string("print(adata)")
#'
#' # Directly convert Seurat object to AnnData object including all processing steps
#' Seu2Adata(pbmc)
#'
#' # Update AnnData with new metadata from Seurat
#' pbmc$cluster2 <- pbmc$cluster
#' adata.AddMetadata(pbmc, col = "cluster2")
#' py_run_string("print(adata)")
#'
#' # Save the AnnData object to a local file
#' pbmc_adata_path <- file.path(tempdir(), "pbmc3k_small.h5ad")
#' adata.Save(pbmc_adata_path)
#'
#' # Load an existing AnnData object from file
#' adata.Load(pbmc_adata_path)
#'
#' @rdname adata.Load
#' @export


adata.Load <- function(
    adata.path,
    conda_env = "seuratextend"
) {
  library(reticulate)
  load_condaenv(conda_env = conda_env)
  reticulate::py_run_string(paste0("
import scanpy as sc
adata = sc.read('", adata.path, "')
"))
}

#' @rdname adata.Load
#' @export

adata.Save <- function(
    adata.path,
    conda_env = "seuratextend"
) {
  library(reticulate)
  load_condaenv(conda_env = conda_env)
  reticulate::py_run_string(paste0("
import scanpy as sc
adata.write('", adata.path, "')
"))
}

#' @rdname adata.Load
#' @export

adata.AddDR <- function(
    seu,
    dr = NULL,
    prefix = "X_",
    postfix = "",
    scv.graph = FALSE,
    load.adata = NULL,
    save.adata = NULL,
    conda_env = "seuratextend"
) {
  library(reticulate)
  load_condaenv(conda_env = conda_env)
  # Initialize Python interface
  py <- reticulate::import_builtins()

  # Check dr input
  dr_list <- Reductions(seu)
  dr <- dr %||% dr_list
  invalid_dr <- dr[!dr %in% dr_list]
  if(length(invalid_dr) > 0){
    stop("Dimension reduction(s) '", paste0(invalid_dr, collapse = ", "), "' not found.\n",
         "Available option(s): ", paste0(dr_list, collapse = ", "))
  }

  # Load adata if load.adata is provided
  if (!is.null(load.adata)) {
    adata.Load(load.adata, conda_env = conda_env)
  }

  # Export all the dr to the adata object
  for (current_dr in dr) {
    embeddings <- Embeddings(seu[[current_dr]])
    py_var_name <- paste0("r_", current_dr)
    py[[py_var_name]] <- embeddings
    reticulate::py_run_string(paste0("
import numpy as np
adata.obsm['", prefix, current_dr, postfix, "'] = np.array(", py_var_name, ")
"))
  }

  # Run velocity_graph if scv.graph is TRUE
  if (scv.graph) {
    reticulate::py_run_string("
import scvelo as scv
scv.tl.velocity(adata)
scv.tl.velocity_graph(adata)
")
  }

  # Save adata if save.adata is provided
  if (!is.null(save.adata)) {
    adata.Save(save.adata, conda_env = conda_env)
  }

  invisible(NULL)
}

#' @rdname adata.Load
#' @export

adata.AddMetadata <- function(
    seu,
    col,
    load.adata = NULL,
    save.adata = NULL,
    conda_env = "seuratextend"
) {
  library(reticulate)
  load_condaenv(conda_env = conda_env)
  # Initialize Python interface
  py <- reticulate::import_builtins()

  # Load adata if load.adata is provided
  if (!is.null(load.adata)) {
    adata.Load(load.adata, conda_env = conda_env)
  }

  # Check if the specified columns exist in the Seurat metadata
  available_cols <- colnames(seu@meta.data)
  if(!all(col %in% available_cols)){
    stop("Metadata column(s) '", paste0(col[!col %in% available_cols], collapse = ", "), "' not found.\n",
         "Available option(s): ", paste0(available_cols, collapse = ", "))
  }

  # Export the specified columns to the adata object
  for (current_col in col) {
    metadata_col <- seu@meta.data[[current_col]]
    py_var_name <- paste0("r_meta_", current_col)
    py[[py_var_name]] <- metadata_col
    reticulate::py_run_string(paste0("
adata.obs['", current_col, "'] = ", py_var_name, "
"))
  }

  # Save adata if save.adata is provided
  if (!is.null(save.adata)) {
    adata.Save(save.adata, conda_env = conda_env)
  }

  invisible(NULL)
}

#' @rdname adata.Load
#' @export

adata.LoadLoom <- function(
    loompath,
    save.adata = NULL,
    conda_env = "seuratextend"
) {
  library(reticulate)
  load_condaenv(conda_env = conda_env)

  # Build Python code string to read loom file using Scanpy
  python_code <- glue::glue("
import scanpy as sc
adata = sc.read_loom('{loompath}')
")
  py_run_string(python_code)

  # Save adata if save.adata is provided
  if (!is.null(save.adata)) {
    adata.Save(save.adata, conda_env = conda_env)
  }

  invisible(NULL)
}

#' @rdname adata.Load
#' @export

Seu2Adata <- function(
    seu,
    add.normdata = TRUE,
    save.adata = NULL,
    conda_env = "seuratextend"
) {
  library(reticulate)
  load_condaenv(conda_env = conda_env)

  loompath <- file.path(tempdir(), "seu.loom")
  Seu2Loom(seu, add.normdata = add.normdata, filename = loompath, overwrite = TRUE)
  adata.LoadLoom(loompath = loompath, conda_env = conda_env)
  adata.AddDR(seu, dr = Reductions(seu), conda_env = conda_env)

  # Save adata if save.adata is provided
  if (!is.null(save.adata)) {
    adata.Save(save.adata, conda_env = conda_env)
  }

  invisible(NULL)
}
