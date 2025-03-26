#' @title Convert Seurat Object to AnnData and Generate scVelo Plots for Single-Cell RNA Velocity Analysis
#' @description This set of functions converts a Seurat object and associated Velocyto loom file(s) into an AnnData object and generates visualization plots for RNA velocity analysis using scVelo. The AnnData object can be directly read from a file or accessed from memory to produce various styles of plots. This integrated approach facilitates the use of scVelo for trajectory analysis in Python's Scanpy library, allowing seamless transition between data processing in R and trajectory analysis in Python.
#' @param seu The Seurat object containing single-cell RNA sequencing data that needs to be analyzed using scVelo.
#' @param filename Path where the resulting AnnData object will be saved. This should be a path to an h5ad file.
#' @param velocyto.loompath Path(s) to the Velocyto-generated loom file which contains RNA velocity data.
#' @param cell.id.match.table An optional data frame for advanced users that maps cell IDs between the Seurat object and Velocyto loom file across multiple samples. It requires a strict format with three columns: cellid.seurat, cellid.velocyto, and velocyto.loompath, indicating the cell ID in the Seurat object, the corresponding cell ID in the Velocyto loom, and the loom file path for that sample, respectively. Default: NULL
#' @param prefix Prefix used to prepend to cell IDs in the Seurat object to match the corresponding IDs in the Velocyto loom file, reflecting sample or batch identifiers. Default: NULL
#' @param postfix Postfix appended to cell IDs in the Seurat object to match the corresponding IDs in the Velocyto loom file. Default: '-1'
#' @param remove_duplicates Logical flag indicating whether to remove duplicate cells in the AnnData object. If TRUE, duplicate cells are removed based on PCA and sum of gene expression values. Default: FALSE
#' @param conda_env Name of the Conda environment where the Python dependencies for scVelo and Scanpy are installed. This environment is used to run Python code from R. Default: 'seuratextend'
#' @return If remove_duplicates = TRUE, returns the filtered Seurat object with duplicate cells removed. Otherwise, does not return any object within R; instead, prepares and stores an AnnData object `adata` in the Python environment accessible via `reticulate`, and generates plots which can be viewed directly or saved to a file. The plots reflect the dynamics of RNA velocity in single-cell datasets.
#' @details This integrated functionality facilitates a seamless transition between converting Seurat objects to AnnData objects and plotting with scVelo. The primary metadata and dimension reduction data from the Seurat object are used to prepare the AnnData object, which is then utilized for generating plots. `SeuratExtend` enhances scVelo plotting capabilities in R, supporting a variety of customization options for visualizing single-cell RNA velocity data. Users can manipulate plot styles, color schemes, group highlights, and more, making it an essential tool for advanced single-cell analysis without the need for direct interaction with Python code.
#' @examples
#' library(Seurat)
#' library(SeuratExtend)
#'
#' # Download the example Seurat Object
#' mye_small <- readRDS(url("https://zenodo.org/records/10944066/files/pbmc10k_mye_small_velocyto.rds", "rb"))
#'
#' # Download the example velocyto loom file to tmp folder
#' loom_path <- file.path(tempdir(), "pbmc10k_mye_small.loom")
#' download.file("https://zenodo.org/records/10944066/files/pbmc10k_mye_small.loom", loom_path)
#'
#' # Set up the path for saving the AnnData object in the HDF5 (h5ad) format
#' adata_path <- file.path(tempdir(), "mye_small.h5ad")
#'
#' # Integrate Seurat Object and velocyto loom information into one AnnData object, which will be stored at the specified path.
#' scVelo.SeuratToAnndata(
#'   mye_small, # The downloaded example Seurat object
#'   filename = adata_path, # Path where the AnnData object will be saved
#'   velocyto.loompath = loom_path, # Path to the loom file
#'   prefix = "sample1_", # Prefix for cell IDs in the Seurat object
#'   postfix = "-1" # Postfix for cell IDs in the Seurat object
#' )
#'
#' # Generate a default UMAP plot colored by 'cluster' and save it as a PNG file
#' scVelo.Plot(color = "cluster", save = "umap1.png", figsize = c(5,4))
#'
#' # Generate a scatter style plot highlighting specific groups, using a custom color palette, with specified axis limits, and save it to a file
#' scVelo.Plot(
#'   style = "scatter",
#'   color = "cluster",
#'   groups = c("DC", "Mono CD14"),
#'   palette = color_pro(3, "light"),
#'   xlim = c(0, 10), ylim = c(0, 10),
#'   save = "umap2_specified_area.png"
#' )
#' @rdname scVelo.SeuratToAnndata
#' @export

scVelo.SeuratToAnndata <- function(
    seu,
    filename,
    velocyto.loompath,
    cell.id.match.table = NULL,
    prefix = NULL,
    postfix = "-1",
    remove_duplicates = FALSE,
    conda_env = "seuratextend"
) {
  library(Seurat)
  import("hdf5r")
  library(glue)
  library(reticulate)

  load_condaenv(conda_env = conda_env)

  if(!is.character(velocyto.loompath)) {
    stop("velocyto.loompath must be a character vector of file paths.")
  }
  if(length(velocyto.loompath) == 0) {
    stop("No valid loom file paths provided in velocyto.loompath.")
  }
  if(any(!file.exists(velocyto.loompath))) {
    stop("Loom file(s) not found: ",
         paste(velocyto.loompath[!file.exists(velocyto.loompath)], collapse = ", "))
  }

  # Initialize lists to store spliced and unspliced matrices
  spliced_list <- list()
  unspliced_list <- list()

  # Loop over each loom file
  for (file in velocyto.loompath) {
    # Open the loom file
    loom <- H5File$new(file, mode = "r")

    # Check if the required datasets and attributes exist in the loom file
    if (!loom$exists("layers/spliced") || !loom$exists("layers/unspliced") ||
        !loom$exists("col_attrs/CellID") || !loom$exists("row_attrs/Gene")) {
      warning(paste("Missing spliced/unspliced data or attributes in loom file:", file))
      next
    }

    # Read spliced and unspliced matrices
    spliced <- loom[["layers/spliced"]][,]
    unspliced <- loom[["layers/unspliced"]][,]

    # Assign row names and column names to the matrices
    cell_ids <- loom[["col_attrs/CellID"]][]
    if(is.null(cell.id.match.table)) {
      # Extract the barcode from the cell IDs
      cell_ids <- sub(".*:", "", cell_ids) # Remove everything before ":"
      cell_ids <- sub("x$", "", cell_ids) # Remove the trailing "x"
    }
    # Modify duplicated genes
    gene_ids <- make.unique(loom[["row_attrs/Gene"]][])

    rownames(spliced) <- cell_ids
    colnames(spliced) <- gene_ids
    rownames(unspliced) <- cell_ids
    colnames(unspliced) <- gene_ids

    spliced_list[[file]] <- spliced
    unspliced_list[[file]] <- unspliced

    # Close the loom file
    loom$close_all()
  }

  # Check if cell.id.match.table is not NULL
  if(!is.null(cell.id.match.table)) {
    # Validate the cell.id.match.table
    if(!all(c("cellid.seurat", "cellid.velocyto", "velocyto.loompath") %in% names(cell.id.match.table)))
      stop("cell.id.match.table must have columns: cellid.seurat, cellid.velocyto, and velocyto.loompath.")

    # Initialize merged spliced and unspliced list
    merged_spliced <- list()
    merged_unspliced <- list()

    # Loop over each loom file
    for(file in unique(cell.id.match.table$velocyto.loompath)) {
      spliced <- spliced_list[[file]]
      unspliced <- unspliced_list[[file]]

      # Get the subset of cell.id.match.table for this file
      subset_table <- subset(cell.id.match.table, velocyto.loompath == file)

      # Create a named vector for renaming: names are cellid.velocyto and values are cellid.seurat
      renaming_vector <- setNames(subset_table$cellid.seurat, subset_table$cellid.velocyto)

      # Filter rows of spliced and unspliced matrices that are in renaming_vector
      spliced <- spliced[rownames(spliced) %in% names(renaming_vector), , drop = FALSE]
      unspliced <- unspliced[rownames(unspliced) %in% names(renaming_vector), , drop = FALSE]

      # Rename the row names of spliced and unspliced matrices
      rownames(spliced) <- renaming_vector[rownames(spliced)]
      rownames(unspliced) <- renaming_vector[rownames(unspliced)]

      # Merge matrices
      merged_spliced[[file]] <- spliced
      merged_unspliced[[file]] <- unspliced
    }
  } else {
    # Validate and replicate prefix and postfix if needed
    if(length(prefix) != length(velocyto.loompath) && length(prefix) != 1)
      stop("Invalid length of prefix. It should be either the same length as velocyto.loompath or 1.")
    if(length(postfix) != length(velocyto.loompath) && length(postfix) != 1)
      stop("Invalid length of postfix. It should be either the same length as velocyto.loompath or 1.")

    if(length(prefix) == 1) prefix <- rep(prefix, length(velocyto.loompath))
    if(length(postfix) == 1) postfix <- rep(postfix, length(velocyto.loompath))

    # Initialize merged spliced and unspliced list
    merged_spliced <- list()
    merged_unspliced <- list()

    # Loop over each loom file
    for(i in seq_along(velocyto.loompath)) {
      spliced <- spliced_list[[i]]
      unspliced <- unspliced_list[[i]]

      # Modify row names using prefix and postfix
      rownames(spliced) <- paste0(prefix[i], rownames(spliced), postfix[i])
      rownames(unspliced) <- paste0(prefix[i], rownames(unspliced), postfix[i])

      # Filter rows based on colnames of Seurat object
      common_cells <- intersect(rownames(spliced), colnames(seu))

      # Check if common_cells is empty or significantly smaller
      if(length(common_cells) == 0)
        stop("No common cells found between the spliced/unspliced matrix and the Seurat object. Please double-check the prefix and postfix.")

      spliced <- spliced[common_cells, , drop = FALSE]
      unspliced <- unspliced[common_cells, , drop = FALSE]

      # Merge matrices
      merged_spliced[[i]] <- spliced
      merged_unspliced[[i]] <- unspliced
    }
  }

  # Identify common genes
  common_genes <- Reduce(intersect, lapply(merged_unspliced, colnames))

  # Subset each matrix to only include the common genes
  merged_unspliced <- lapply(merged_unspliced, function(mat) mat[, common_genes, drop = FALSE])
  merged_spliced <- lapply(merged_spliced, function(mat) mat[, common_genes, drop = FALSE])

  merged_spliced <- rlist::list.rbind(merged_spliced)
  merged_unspliced <- rlist::list.rbind(merged_unspliced)

  # Identify common cells and genes between merged matrices and Seurat object
  common_cells <- intersect(colnames(seu), rownames(merged_spliced))
  common_genes <- intersect(rownames(seu), colnames(merged_spliced))

  # Check if the number of common cells is less than the number of cells in Seurat object
  if(length(common_cells) < length(colnames(seu))) {
    warning(paste0("The number of common cells in the merged matrices (", length(common_cells),
                   ") is less than the number of cells in the Seurat object (", length(colnames(seu)),
                   "). Please double-check the input, including the order of the samples."))
    seu <- subset(seu, cells = common_cells)
  }

  # Check if the number of common genes is less than the number of genes in Seurat object
  if(length(common_genes) < length(rownames(seu))) {
    warning(paste0(length(rownames(seu)) - length(common_genes), " genes in the Seurat object not detected in the Velocyto matrix. Please double-check if the Seurat gene count matrix and Velocyto use the same reference genome."))
    seu <- subset(seu, features = common_genes)
  }

  # Reorder the rows and columns of merged_spliced and merged_unspliced
  merged_spliced <- merged_spliced[common_cells, common_genes]
  merged_unspliced <- merged_unspliced[common_cells, common_genes]

  # Create a temporary directory with a timestamp
  subfolder <- create_temp_dir()

  # Export loom to tmp folder
  tmp.loom.path <- file.path(subfolder, "scvelo.loom")
  Seu2Loom(
    seu, filename = tmp.loom.path,
    layers = list(spliced = merged_spliced, unspliced = merged_unspliced))

  python_code <- glue('
import scvelo as scv
import scanpy as sc
import numpy as np
import matplotlib.cm as cm

# Set scVelo settings
print(f"scVelo version: {{scv.__version__}}")
scv.settings.verbosity = 3
scv.settings.presenter_view = True
scv.set_figure_params("scvelo")

# Read the loom file
adata = sc.read("{tmp.loom.path}")
# Pre-process the data
scv.pp.filter_and_normalize(adata, min_shared_counts=20, n_top_genes=2000)
')

  py_run_string(python_code)

  # Get the list of dimension reductions
  dr_list <- Reductions(seu)

  # Check if dr_list is not empty
  if (length(dr_list) > 0) {
    adata.AddDR(seu, prefix = "X_", postfix = "", conda_env = conda_env)
    adata.AddDR(seu, prefix = "", postfix = "_cell_embeddings", conda_env = conda_env)
  }

  # Run Python code to check for duplicates
  python_code <- glue('
def get_duplicate_cells(data):
  """Check for duplicate cells in AnnData object."""
  from anndata import AnnData
  from scipy.sparse import issparse
  from collections import Counter
  import numpy as np
  import pandas as pd

  if isinstance(data, AnnData):
    X = data.X
    # Get initial size and PCA sum for comparison
    lst = list(np.sum(np.abs(data.obsm["X_pca"]), 1) + np.sum(X.toarray() if issparse(X) else X, 1))
  else:
    X = data
    lst = list(np.sum(X, 1).A1 if issparse(X) else np.sum(X, 1))

  idx_dup = []
  if len(set(lst)) < len(lst):
    vals = [val for val, count in Counter(lst).items() if count > 1]
    idx_dup = np.where(pd.Series(lst).isin(vals))[0]

    X_new = np.array(X[idx_dup].toarray() if issparse(X) else X[idx_dup])
    sorted_idx = np.lexsort(X_new.T)
    sorted_data = X_new[sorted_idx, :]

    row_mask = np.invert(np.append([True], np.any(np.diff(sorted_data, axis=0), 1)))
    idx = sorted_idx[row_mask]
    idx_dup = np.array(idx_dup)[idx]
  return len(idx_dup)

if "X_pca" not in adata.obsm.keys():
      sc.pp.pca(adata)

duplicates = get_duplicate_cells(adata)
')
  py_run_string(python_code)
  duplicates <- py_eval("duplicates")
  if (duplicates > 0 && remove_duplicates) {
    python_code <- "scv.pp.remove_duplicate_cells(adata)"
    py_run_string(python_code)
    message(paste0("Removed ", duplicates, " duplicate cells"))

    # Get remaining cells from adata
    remaining_cells <- py_eval("adata.obs_names.to_list()")
    # Subset Seurat object
    seu <- subset(seu, cells = remaining_cells)
  } else if (duplicates > 0) {
    message(paste0("Found ", duplicates, " duplicate cells, but skipping removal as remove_duplicates=FALSE"))
  }
  # Construct the Python code string to run scVelo and save the Anndata object
  python_code <- glue('
scv.pp.moments(adata, n_pcs=30, n_neighbors=30)
# Compute velocity and velocity graph
scv.tl.velocity(adata)
scv.tl.velocity_graph(adata)
adata.write("{filename}")
')

  # Run the Python code
  reticulate::py_run_string(python_code)

  if (remove_duplicates) {
    return(seu)
  } else {
    return(NULL)
  }
}

#' @title (Deprecated) Export Seurat Object and Velocyto Data to Loom for scVelo Analysis
#' @description This function converts a Seurat object and associated Velocyto loom files into a single loom file that integrates both spliced and unspliced RNA matrices. It is part of an older workflow and is soft deprecated in favor of `scVelo.SeuratToAnndata`. However, it is preserved for reproducibility of earlier work that cited its use.
#' @param seu Seurat object containing single-cell RNA sequencing data.
#' @param export.name Name of the resulting loom file to be created.
#' @param velocyto.loom.path Directory containing the Velocyto loom files; defaults to the current working directory. Default: getwd()
#' @param velocyto.loom.filenames Filenames of Velocyto output (loom files) that contain spliced and unspliced RNA data.
#' @param fxn.convert.loomtoseurat.cellname (Experimental) Function to convert cell names in Velocyto loom files to Seurat object cell names, if they differ.
#' @param prefix Optional prefix added to Seurat cell IDs, typically used to include sample names. Default: NULL
#' @param postfix Optional postfix added to Seurat cell IDs, typically used to denote unique identifiers such as '-1'. Default: NULL
#' @return Generates a loom file that integrates Seurat and Velocyto data, stored in the specified export path.
#' @details This set of functions supports the transition from Seurat objects to an integrated loom file via `scVelo.SeuratToLoom`, and from loom files to AnnData objects via `scVelo.RunBasic`, ready for basic to advanced scVelo analysis in Python. This older method is preserved for compatibility with previous studies and publications, ensuring reproducibility of earlier work. This older method is preserved for compatibility with previous studies and publications, but users are encouraged to transition to `scVelo.SeuratToAnndata` for a more streamlined and updated approach that directly integrates Seurat objects with scVelo analysis in Python's Scanpy library.
#' @seealso \code{\link[SeuratExtend:scVelo.SeuratToAnndata]{scVelo.SeuratToAnndata}} for the recommended method of preparing data for scVelo analysis.
#' @rdname scVelo-SeuratToLoom-and-RunBasic
#' @export

scVelo.SeuratToLoom <- function(
    seu,
    export.name,
    velocyto.loom.path = NULL,
    velocyto.loom.filenames,
    fxn.convert.loomtoseurat.cellname = NULL,
    cell.id.match.table = NULL,
    prefix = NULL,
    postfix = NULL
){
  library(Seurat)
  library(rlist)
  library(dplyr)
  library(rlang)
  import("loomR")
  import("SeuratDisk")

  check.metadata.na <- sapply(seu@meta.data, function(x) any(is.na(x)))
  if(any(check.metadata.na)) {
    message(paste0("NA(s) exist in meta data: ",
                   paste0(check.metadata.na %>% names(.)[.], collapse = ", "),
                   "\nColumn(s) removed"))
    seu@meta.data[,check.metadata.na] <- NULL
  }

  if(!is.null(velocyto.loom.path)) {
    loom.files <- file.path(velocyto.loom.path, velocyto.loom.filenames)
  } else {
    loom.files <- velocyto.loom.filenames
  }
  if(any(!file.exists(loom.files))) {
    stop(paste0("Loom file(s) not found: ",
                paste0(loom.files[!file.exists(loom.files)], collapse = ", ")
    )
    )
  }

  # orig.assay <- DefaultAssay(seu)
  DefaultAssay(seu) <- "RNA"
  genes <- rownames(seu)
  cells <- colnames(seu)
  seu.subset.genes <- F

  lfiles <- list()
  cellname <- list()
  # ambiguous <- list()
  spliced <- list()
  unspliced <- list()

  for (i in seq_along(velocyto.loom.filenames)) {
    lfiles[[i]] <- connect(filename = loom.files[i],
                           mode = "r+", skip.validate = T)

    if(is.null(fxn.convert.loomtoseurat.cellname)){
      cellname[[i]] <-
        strsplit(lfiles[[i]][["col_attrs"]][["CellID"]][], split = ":", fixed = T) %>%
        sapply(function(x) x[2])
      cellname[[i]] <- paste0(prefix[i],sub("x", "", cellname[[i]]),postfix[i])
    }else{
      cellname[[i]] <- fxn.convert.loomtoseurat.cellname(
        velocyto.loom.filenames = i,
        loom.cellname = lfiles[[i]][["col_attrs"]][["CellID"]][])
    }
    if(!any(cellname[[i]] %in% cells)) stop(paste0(velocyto.loom.filenames[i],": wrong prefix/postfix?"))

    g <- lfiles[[i]][["row_attrs"]][["Gene"]][] %>% make.unique()
    n <- length(setdiff(genes, g))
    if(n>0){
      message("There are ", n, " genes in the Seurat object which cannot be found in the Velocyto output. ",
              "Did you use different version of reference genome?")
      genes <- intersect(genes, g)
      seu.subset.genes <- T
    }
    # ambiguous[[i]] <-
    #   lfiles[[i]][["layers"]][["ambiguous"]][,] %>%
    #   `colnames<-`(g) %>%
    #   `rownames<-`(cellname[[i]]) %>%
    #   .[rownames(.) %in% cells, genes]
    spliced[[i]] <-
      lfiles[[i]][["layers"]][["spliced"]][,] %>%
      `colnames<-`(g) %>%
      `rownames<-`(cellname[[i]]) %>%
      .[!is.na(rownames(.)),] %>%
      .[rownames(.) %in% cells, genes]
    unspliced[[i]] <-
      lfiles[[i]][["layers"]][["unspliced"]][,] %>%
      `colnames<-`(g) %>%
      `rownames<-`(cellname[[i]]) %>%
      .[rownames(.) %in% cells, genes]
  }

  # ambiguous <- ambiguous %>%
  #   list.rbind() %>%
  #   .[cells,]
  # spliced <- spliced %>%
  #   list.rbind() %>%
  #   .[cells,]
  # unspliced <- unspliced %>%
  #   list.rbind() %>%
  #   .[cells,]
  spliced <- list.rbind(spliced)
  unspliced <- list.rbind(unspliced)
  cell.looms <- rownames(spliced)
  nc <- length(setdiff(cells, cell.looms))
  if(nc>0){
    message("There are ", nc, " cells in the Seurat object which cannot be found in the Velocyto output. ")
    cells <- intersect(cells, cell.looms)
  }
  spliced <- spliced[cells,]
  unspliced <- unspliced[cells,]

  seu.new <- subset(seu, cells = cells, features = genes)
  seu.new <- CreateSeuratObject(GetAssayData(seu.new, slot = "counts"), meta.data = seu.new@meta.data)
  # if(add.reduction) {
  #   seu.new <- NormalizeData(seu.new)
  #   seu.new <- FindVariableFeatures(seu.new)
  #   seu.new <- ScaleData(seu.new)
  #   for (i in names(seu@reductions)) {
  #     seu.new[[i]] <- CreateDimReducObject(embeddings = Embeddings(seu, reduction = i))
  #   }
  # }
  pfile <- as.loom(x = seu.new, filename = export.name, overwrite = T)

  pfile$add_layer(x = t(spliced), name = "spliced")
  pfile$add_layer(x = t(unspliced), name = "unspliced")

  pfile$close_all()
}

#' @param loom Path to the loom file created by `scVelo.SeuratToLoom`.
#' @param save.adata Name and path where the AnnData object should be saved, typically ending in '.h5ad'. Default: 'adata.obj'
#' @rdname scVelo-SeuratToLoom-and-RunBasic
#' @export

scVelo.RunBasic <- function(loom, save.adata = "adata.obj"){
  library(reticulate)
  code <- paste0("import scvelo as scv\n",
                 "import scanpy as sc\n",
                 "import numpy as np\n",
                 "import matplotlib.cm as cm\n",
                 'print(f"scVelo version: {scv.__version__}")\n',
                 "scv.settings.verbosity = 3\n",
                 "scv.settings.presenter_view = True\n",
                 "scv.set_figure_params('scvelo')\n",
                 "import loompy\n",
                 'ds = loompy.connect("', loom, '")\n',
                 'del ds.layers["norm_data"]\n',
                 'del ds.layers["scale_data"]\n',
                 'ds.close()\n',
                 'adata = sc.read("', loom, '")\n',
                 "adata\n",
                 "scv.pp.filter_and_normalize(adata, min_shared_counts=20, n_top_genes=2000)\n",
                 "scv.pp.moments(adata, n_pcs=30, n_neighbors=30)\n",
                 "scv.tl.velocity(adata)\n",
                 "scv.tl.velocity_graph(adata)\n",
                 "import pickle\n",
                 'pickle_out = open("', save.adata, '","wb")\n',
                 "pickle.dump(adata, pickle_out)\n",
                 "pickle_out.close()\n")
  py_run_string(code)
}

#' @param load.adata Path to a previously saved AnnData object (in h5ad format) which can be directly loaded to avoid re-running preprocessing. If NULL, reticulate will automatically use the existing AnnData object `adata` in the Python environment for plotting. Default: NULL.
#' @param style Style of the velocity plot, allowing for different visual representations such as 'stream', 'grid', or 'scatter'. Default: c("stream", "grid", "scatter").
#' @param basis The embedding to be used for plotting, typically 'umap_cell_embeddings' to represent UMAP reductions. Default: 'umap_cell_embeddings'.
#' @param color The variable by which to color the plot, usually a categorical variable like cluster identifiers or a continuous variable reflecting gene expression levels. Default: NULL.
#' @param groups Groups or clusters to highlight in the plot, useful for focusing on specific cell types or conditions within the dataset. Default: NULL.
#' @param palette Color palette to use for differentiating between groups or clusters within the plot. Allows customization of aesthetic presentation. Default: NULL.
#' @param alpha Opacity of the points in the plot, which can be adjusted to enhance visualization when dealing with densely packed points. Default: 0.15.
#' @param arrow_size Size of the arrows representing RNA velocity vectors in the plot, relevant only when `style` is set to 'scatter'. This can be adjusted to make the arrows more or less prominent based on visualization needs. Default: 3.
#' @param arrow_length Length of the arrows, which affects how far the arrows extend from their origin points. Relevant only when style is 'scatter', helping in interpreting the directionality and magnitude of cellular transitions. Default: 2.
#' @param dpi Resolution of the saved plot, useful when preparing figures for publication or presentations. Default: 300.
#' @param legend_fontsize Size of the font used in the plot legend, allowing for customization based on the figure's intended use or audience. Default: 9.
#' @param figsize Dimensions of the plot in inches, providing control over the size of the output figure to accommodate different analysis contexts. Default: c(7, 5).
#' @param xlim Limits for the x-axis, which can be set to focus on specific areas of the plot or to standardize across multiple plots. Default: NULL.
#' @param ylim Limits for the y-axis, similar in use to `xlim` for focusing or standardizing the y-axis view. Default: NULL.
#' @param save Path where the plot should be saved. If specified, the plot will be saved to the given location. Supports various file formats like PNG, PDF, SVG, etc. Default: NULL.
#' @rdname scVelo.SeuratToAnndata
#' @export

scVelo.Plot <- function(
    load.adata = NULL,
    style = c("stream", "grid", "scatter"),
    basis = "umap_cell_embeddings",
    color = NULL,
    groups = NULL,
    palette = NULL,
    alpha = 0.15,
    arrow_size = 3,
    arrow_length = 2,
    dpi = 300,
    legend_fontsize = 9,
    figsize = c(7, 5),
    xlim = NULL,
    ylim = NULL,
    save = NULL,
    conda_env = "seuratextend"
) {
  library(reticulate)
  library(glue)

  load_condaenv(conda_env = conda_env)

  # Load adata if load.adata is provided
  if (!is.null(load.adata)) {
    adata.Load(load.adata, conda_env = conda_env)
  }

  # Determine the scv plotting function
  plot_function <- switch(
    style[1],
    "stream" = "scv.pl.velocity_embedding_stream",
    "grid" = "scv.pl.velocity_embedding_grid",
    "scatter" = "scv.pl.velocity_embedding",
    stop("Invalid style provided.")
  )

  # Construct the function call
  args_list <- list(
    basis = shQuote(basis, type = "sh"),
    color = if (!is.null(color)) shQuote(color, type = "sh") else NULL,
    groups = r_vector_to_py(groups),
    palette = r_vector_to_py(palette),
    alpha = alpha,
    dpi = dpi,
    legend_fontsize = legend_fontsize,
    save = if (!is.null(save)) shQuote(save, type = "sh") else NULL,
    figsize = r_vector_to_py(figsize, type = "tuple"),
    xlim = r_vector_to_py(xlim, type = "tuple"),
    ylim = r_vector_to_py(ylim, type = "tuple")
  )

  if (style[1] == "scatter") {
    args_list$arrow_size <- arrow_size
    args_list$arrow_length <- arrow_length
  }

  # Filter out NULL arguments
  args_list <- args_list[!sapply(args_list, is.null)]

  # Convert to Python function call string
  args_string <- paste(paste(names(args_list), "=", args_list, sep = ""), collapse = ", ")
  py_command <- glue("{plot_function}(adata, {args_string})")

  # Execute the Python command
  reticulate::py_run_string("import scvelo as scv")
  reticulate::py_run_string(py_command)

  invisible(NULL)
}


