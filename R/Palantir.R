#' @title Run Palantir Diffusion Map and Calculate Pseudotime
#' @description This function suite uses the Palantir algorithm to first calculate the diffusion map based on pre-calculated dimension reductions in Seurat (e.g., PCA, harmony), adding the diffusion map (dm) and multiscale space (ms) embeddings back to the original Seurat object. Subsequently, it calculates pseudotime to determine the developmental trajectory and fate decisions of cells based on specified starting points.
#' @param seu A Seurat object.
#' @param reduction The name the dimension reduction to use. Default: 'pca'.
#' @param n_components An integer specifying how many dimensions to use as input. Default: 50.
#' @param conda_env The name of the Conda environment that contains the necessary Python libraries to execute Palantir. Default: 'seuratextend'
#' @return The Seurat object is updated to include the diffusion map and multiscale space embeddings in the respective slots. Additionally, it adds the calculated pseudotime information to the `SeuratObj@misc$Palantir` slot of the Seurat object.
#' @details These functions integrate the powerful trajectory analysis capabilities of Palantir with the data handling and visualization strengths of Seurat. Initially, `Palantir.RunDM` calculates diffusion maps that represent cellular states and transitions in a high-dimensional space, which are then applied to visually explore cell trajectories in reduced dimensional space. `Palantir.Pseudotime` further leverages these embeddings to model developmental trajectories, allowing for a detailed understanding of cell fate decisions. This comprehensive analysis is crucial for identifying key transitions and stages in cellular development, making it an invaluable tool for developmental biology studies.
#' @examples
#' library(Seurat)
#' library(SeuratExtend)
#'
#' # Download the example Seurat Object
#' mye_small <- readRDS(url("https://zenodo.org/records/10944066/files/pbmc10k_mye_small_velocyto.rds", "rb"))
#'
#' # Running Diffusion Map with Palantir
#' mye_small <- Palantir.RunDM(mye_small)
#' # View the new dimensional reductions "dm" and "ms" in the Seurat Object:
#' print(mye_small@reductions)
#' # View the first 2 ms dimensions:
#' DimPlot2(mye_small, reduction = "ms", group.by = "cluster", label = TRUE)
#'
#' # Optional: Running UMAP based on "ms" if having more than 2 dimensions
#' # Assuming 'ms' has 10 dimensions:
#' # mye_small <- RunUMAP(mye_small, reduction = "ms", dims = 1:10)
#' # DimPlot2(mye_small, group.by = "cluster", label = TRUE, reduction = "umap")
#'
#' # Determining Cell Fates and Calculating Pseudotime. You can manually select "start cells" using the `CellSelector` function.
#' # p <- DimPlot(mye_small, reduction = "ms", group.by = "cluster")
#' # cells <- CellSelector(p)
#' # print(head(cells, 3))
#'
#' # Calculating pseudotime and updating the Seurat object
#' mye_small <- Palantir.Pseudotime(mye_small, start_cell = "sample1_GAGAGGTAGCAGTACG-1")
#' ps <- mye_small@misc$Palantir$Pseudotime
#' print(head(ps))
#' # Visualize cell fate on UMAP
#' colnames(ps)[3:4] <- c("fate1", "fate2")
#' mye_small@meta.data[,colnames(ps)] <- ps
#' DimPlot2(mye_small, features = colnames(ps), reduction = "ms",
#'          cols = list(Entropy = "D"))
#' @rdname Palantir-Methods
#' @export

Palantir.RunDM <- function(
    seu,
    reduction = "pca",
    n_components = 50,
    conda_env = "seuratextend"
) {
  library(Seurat)

  load_condaenv(conda_env = conda_env)

  # Create a temporary directory with a timestamp
  subfolder <- create_temp_dir()

  # Write the CSV file to the temporary subfolder
  csv_file_path <- file.path(subfolder, "dr.csv")
  write.csv(Embeddings(seu, reduction = reduction), csv_file_path)

  code <- glue::glue("
import palantir
import pandas as pd

# Read the PCA projections from the CSV file
pca_projections = pd.read_csv('{csv_file_path}', index_col=0)

# Run diffusion maps and save the results
dm_res = palantir.utils.run_diffusion_maps(pca_projections, n_components={n_components})
dm_res['EigenVectors'].to_csv('{subfolder}/dm_res.csv')

# Determine multiscale space and save the results
ms_data = palantir.utils.determine_multiscale_space(dm_res)
ms_data.to_csv('{subfolder}/ms_data.csv')
")
  py_run_string(code)

  # Define the paths to the CSV files
  dm_res_path <- file.path(subfolder, "dm_res.csv")
  ms_data_path <- file.path(subfolder, "ms_data.csv")

  # Read the CSV files and create DimReduc objects
  dm_data <- read.csv(dm_res_path, row.names = 1)
  colnames(dm_data) <- paste0("DM_", 1:ncol(dm_data))
  seu[["dm"]] <- CreateDimReducObject(
    embeddings = as.matrix(dm_data),
    key = "DM_",
    assay = DefaultAssay(seu)
  )

  ms_data <- read.csv(ms_data_path, row.names = 1)
  colnames(ms_data) <- paste0("MS_", 1:ncol(ms_data))
  seu[["ms"]] <- CreateDimReducObject(
    embeddings = as.matrix(ms_data),
    key = "MS_",
    assay = DefaultAssay(seu)
  )
  return(seu)
}

#' @title (Deprecated) Run Palantir Diffusion Map
#' @description
#' This is a deprecated version of the Palantir Diffusion Map function.
#' Please refer to the updated function, \code{\link{Palantir.RunDM}}, for the latest features and improvements.
#' This older version, in addition to calculating the diffusion map (dm) and multiscale_space (ms),
#' also calculates t-SNE. To use the t-SNE feature, the \code{multicore-tsne} package needs to be installed in the conda environment.
#' @param seu A Seurat object.
#' @param reduction A character string specifying the dimension reduction to use. Default: 'pca'.
#' @param n_components An integer specifying how many dimensions to use as input. Default: 20.
#' @return A Seurat object with added dm, ms, and t-SNE embeddings.
#' @examples
#' \dontrun{
#' library(Seurat)
#' library(SeuratExtend)
#' pbmc <- RunPalantirDiffusionMap(pbmc)
#' DimPlot(pbmc, reduction = "tsne", label = TRUE)
#' }
#' @rdname RunPalantirDiffusionMap
#' @export
#' @seealso \code{\link{Palantir.RunDM}} for the updated version of this function.

RunPalantirDiffusionMap <- function(seu, reduction = "pca", n_components = 20) {
  library(Seurat)
  library(reticulate)
  library(dplyr)
  library(magrittr)

  if(!"tmp" %in% list.files()) dir.create("tmp")
  write.csv(Embeddings(seu, reduction = reduction), "tmp/dr.csv")

  code <- paste0("import palantir\n",
                 "import os\n",
                 "import matplotlib\n",
                 "import matplotlib.pyplot as plt\n",
                 "import pandas as pd\n",
                 "import numpy as np\n",
                 "import random\n\n",
                 "pca_projections = pd.read_csv('tmp/dr.csv', index_col=0)\n",
                 "dm_res = palantir.utils.run_diffusion_maps(pca_projections, n_components=",
                 n_components,
                 ")\n",
                 "dm_res['EigenVectors'].to_csv('tmp/dm_res.csv')\n",
                 "ms_data = palantir.utils.determine_multiscale_space(dm_res)\n",
                 "ms_data.to_csv('tmp/ms_data.csv')\n",
                 "tsne = palantir.utils.run_tsne(ms_data)\n",
                 "tsne.to_csv('tmp/tsne.csv')\n",
                 "quit")
  py_run_string(code)
  seu[["dm"]] <-
    read.csv("tmp/dm_res.csv", row.names = 1) %>%
    set_colnames(paste0("DM_", 1:ncol(.))) %>%
    as.matrix() %>%
    CreateDimReducObject(key = "DM_", assay = DefaultAssay(seu))
  seu[["ms"]] <-
    read.csv("tmp/ms_data.csv", row.names = 1) %>%
    set_colnames(paste0("MS_", 1:ncol(.))) %>%
    as.matrix() %>%
    CreateDimReducObject(key = "MS_", assay = DefaultAssay(seu))
  seu[["tsne"]] <-
    read.csv("tmp/tsne.csv", row.names = 1) %>%
    set_colnames(paste0("TSNE_", 1:ncol(.))) %>%
    as.matrix() %>%
    CreateDimReducObject(key = "TSNE_", assay = DefaultAssay(seu))
  return(seu)
}

#' @param start_cell A vector of cell identifiers from the Seurat object that marks the starting cells for the trajectory analysis. These cells are typically identified as progenitor or early developmental states in the dataset.
#' @param terminal_states Optional; a named vector of cell identifiers that represent terminal states in the developmental trajectory, with names reflecting the fate labels, such as c("fate1" = "sample1_AAACCCAAGGCCCAAA-1", "fate2" = "sample1_CCTCTAGGTGAACGGT-1"). These terminal states help Palantir more accurately model the trajectory towards cellular differentiation. Default: NULL
#' @param title A label used to store and retrieve pseudotime information within the `SeuratObj@misc$Palantir` slot of the Seurat object. Default: 'Pseudotime'
#' @param n_jobs Number of jobs for parallel computation. Default: 1
#' @rdname Palantir-Methods
#' @export

Palantir.Pseudotime <- function(
    seu,
    start_cell,
    terminal_states = NULL,
    title = "Pseudotime",
    conda_env = "seuratextend",
    n_jobs = 1
) {
  library(Seurat)

  load_condaenv(conda_env = conda_env)

  # Check if "ms" is in the dimension reduction of the Seurat object
  if(!"ms" %in% names(seu@reductions)) {
    stop("Multiscale_space (ms) embeddings not found in the Seurat object. Please run Palantir.RunDM() to calculate them before proceeding.")
  }

  # Create a temporary directory with a timestamp
  subfolder <- create_temp_dir()

  # Write the CSV file to the temporary subfolder
  csv_file_path <- file.path(subfolder, "ms.csv")
  write.csv(Embeddings(seu, reduction = "ms"), csv_file_path)

  # Define the common part of the Python code string
  common_code <- glue::glue("
import palantir
import pandas as pd

# Read the MS from the CSV file
ms_data = pd.read_csv('{csv_file_path}', index_col=0)

# Define start_cell
start_cell = '{start_cell[1]}'
")

  # Decide which specific Python code to append based on the terminal_states parameter
  if(is.null(terminal_states)) {
    specific_code <- glue::glue("pr_res = palantir.core.run_palantir(ms_data, start_cell, num_waypoints=500, n_jobs={n_jobs})")
  } else {
    # Convert the terminal_states R named vector to a Python dictionary string
    terminal_states_dict_str <- toString(
      sapply(names(terminal_states),
             function(name)
               sprintf("'%s': '%s'", terminal_states[name], name),
             USE.NAMES = FALSE),
      collapse = ", ")

    # Construct the specific Python code string with the terminal_states pandas Series and run_palantir function
    specific_code <- glue::glue("
# Create the terminal_states pandas Series
terminal_states_dict = {{ {terminal_states_dict_str} }}
terminal_states = pd.Series(terminal_states_dict)

# Run diffusion maps with start_cell and terminal_states and save the results
pr_res = palantir.core.run_palantir(ms_data, start_cell, num_waypoints=500, terminal_states=terminal_states.index, n_jobs={n_jobs})
pr_res.branch_probs.columns = terminal_states[pr_res.branch_probs.columns]
")
  }

  code2 <- glue::glue("
# Combine the results into one DataFrame
result_df = pd.DataFrame({{
    'Pseudotime': pr_res.pseudotime,
    'Entropy': pr_res.entropy
}})
result_df = pd.concat([result_df, pr_res.branch_probs], axis=1)

# Write the combined DataFrame to a CSV file
result_df.to_csv('{subfolder}/pr_res.csv')
")
  # Combine the common and specific Python code strings and run the Python code
  code <- paste(common_code, specific_code, code2, sep = "\n")

  py_run_string(code)

  # Define the paths to the CSV file
  pr_res_path <- file.path(subfolder, "pr_res.csv")

  # Read the CSV files and saved to Seurat Object
  pr_data <- read.csv(pr_res_path, row.names = 1)
  seu@misc[["Palantir"]][[title]] <- pr_data

  return(seu)
}

#' @title Run MAGIC for Gene Expression Smoothing
#' @description Applies the MAGIC (Markov Affinity-based Graph Imputation of Cells) algorithm to a Seurat object to denoise and smooth gene expression data. This technique enhances data interpretability by imputing missing data and reducing technical noise.
#' @param seu A Seurat object containing single-cell RNA-seq data to which the MAGIC algorithm will be applied.
#' @param n_top_genes The number of top variable genes to consider in the MAGIC algorithm, which helps in focusing the smoothing on the most informative genes. Default: 2000
#' @param n_components The number of principal components to use in dimensionality reduction before applying MAGIC. Useful for preprocessing the data to enhance the effects of MAGIC. Default: 20
#' @param conda_env The name of the Conda environment where the necessary Python dependencies for running MAGIC are installed. This environment is used to run Python code from R, ensuring smooth integration of the two platforms. Default: 'seuratextend'
#' @param n_jobs Number of jobs for parallel computation. Default: 1
#' @return Updates the provided Seurat object by adding a new assay named 'magic', which contains the denoised and smoothed gene expression data.
#' @details MAGIC uses a graph-based approach to infer and smooth gene expression across similar cells, effectively filling in gaps in the data where measurements are sparse or noisy. This process is especially beneficial in datasets with high levels of technical noise or when trying to resolve subtle biological signals that might be obscured by this noise.
#' @examples
#' library(Seurat)
#' library(SeuratExtend)
#'
#' # Download the example Seurat Object
#' mye_small <- readRDS(url("https://zenodo.org/records/10944066/files/pbmc10k_mye_small_velocyto.rds", "rb"))
#'
#' # Run MAGIC
#' mye_small <- Palantir.Magic(mye_small)
#'
#' # Visualizing the effects of MAGIC on selected genes
#' DimPlot2(mye_small, features = c("CD14", "magic_CD14", "FLT3", "magic_FLT3"))
#'
#' @rdname Palantir-Magic
#' @export

Palantir.Magic <- function(
    seu,
    n_top_genes = 2000,
    n_components = 20,
    conda_env = "seuratextend",
    n_jobs = 1
) {
  library(Seurat)
  library(glue)

  load_condaenv(conda_env = conda_env)

  # Create a temporary directory with a timestamp
  subfolder <- create_temp_dir()

  # Export loom to tmp folder
  tmp.loom.path <- file.path(subfolder, "seu.loom")
  Seu2Loom(seu, add.normdata = FALSE, add.metadata = FALSE, filename = tmp.loom.path)

  # Run Python code using reticulate
  py_run_string(glue("
import scanpy as sc
import palantir

ad_magic = sc.read_loom('{tmp.loom.path}')
sc.pp.normalize_per_cell(ad_magic)
palantir.preprocess.log_transform(ad_magic)
sc.pp.highly_variable_genes(ad_magic, n_top_genes={n_top_genes})
sc.pp.pca(ad_magic)
palantir.utils.run_diffusion_maps(ad_magic, n_components={n_components})
imputed_X = palantir.utils.run_magic_imputation(ad_magic, n_jobs={n_jobs})
"))

  # Add magic_imputation to "magic" assay
  magic_matr <- py$imputed_X
  rownames(magic_matr) <- colnames(seu)
  colnames(magic_matr) <- rownames(seu)

  seu[["magic"]] <- CreateAssayObject(data = t(as.matrix(magic_matr)))

  return(seu)
}
