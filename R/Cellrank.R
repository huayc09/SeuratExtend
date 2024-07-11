#' @title Compute and Visualize Cell Trajectories Using CellRank
#' @description `Cellrank.Compute()` calculates cell trajectories using pre-existing pseudotime data in an AnnData object, providing an alternative to scVelo when it produces trajectories that may not align with established biological knowledge. `Cellrank.Plot()` visualizes these trajectories. While CellRank allows for trajectory modeling by choosing starting cells, it is crucial to base these decisions on validated biological insights to prevent misinterpretation of the data. For detailed examples, see the [CellRank documentation](https://cellrank.readthedocs.io/en/latest/notebooks/tutorials/kernels/300_pseudotime.html).
#' @param load.adata Path to a previously saved AnnData object (in h5ad format) which can be directly loaded to avoid re-running preprocessing. If NULL, reticulate will automatically use the existing AnnData object `adata` in the Python environment for plotting. Default: NULL.
#' @param time_key The key used to access the pseudotime data within the AnnData object, which is crucial for trajectory computation. This should match the column name in `adata` where the pseudotime data is stored.
#' @param conda_env Name of the Conda environment where the Python dependencies for cellrank and Scanpy are installed. This environment is used to run Python code from R, ensuring smooth integration and execution of the analysis. Default: 'seuratextend'
#' @examples
#' \dontrun{
#' library(Seurat)
#' library(SeuratExtend)
#'
#' # Load an example Seurat Object
#' mye_small <- readRDS(url("https://zenodo.org/records/10944066/files/pbmc10k_mye_small_velocyto.rds", "rb"))
#'
#' # Calculate diffusion map and pseudotime using Palantir
#' mye_small <- Palantir.RunDM(mye_small)
#' mye_small <- Palantir.Pseudotime(mye_small, start_cell = "sample1_GAGAGGTAGCAGTACG-1")
#'
#' # Retrieve pseudotime values and store them in the meta.data for easy access
#' ps <- mye_small@misc$Palantir$Pseudotime
#' mye_small$Pseudotime <- ps$Pseudotime
#'
#' # Convert the Seurat object to an AnnData object
#' Seu2Adata(mye_small)
#' # Compute cell trajectories using CellRank based on pseudotime
#' Cellrank.Compute(time_key = "Pseudotime")
#'
#' # Visualize cell trajectories using CellRank
#' # 'ms' dimension reduction is used for plotting
#' Cellrank.Plot(color = "cluster", basis = "ms")
#' }
#' @return These functions do not return any object directly. `Cellrank.Compute()` updates the AnnData object with the computed trajectories. `Cellrank.Plot()` generates visual representations of these trajectories in the AnnData object and can directly display or save the plots.
#' @details `Cellrank.Compute()` uses the provided pseudotime data to model the likelihood of cellular transitions, identifying potential pathways and fates within the developmental continuum. This step is critical for accurately capturing the dynamic nature of cell differentiation. `Cellrank.Plot()` then allows for a visual exploration of these pathways, highlighting differences and patterns that can guide further biological interpretation and analysis. Together, these functions provide a robust framework for trajectory analysis in single-cell studies.
#' @rdname Cellrank-Methods
#' @export

Cellrank.Compute <- function(
    load.adata = NULL,
    time_key,
    conda_env = "seuratextend"
) {
  library(reticulate)

  load_condaenv(conda_env = conda_env)

  # Initialize Python interface
  py <- reticulate::import_builtins()

  # Load adata if load.adata is provided
  if (!is.null(load.adata)) {
    adata.Load(load.adata)
  }

  # Pass the time_key parameter to Python
  py$time_key_from_r <- time_key

  # Execute the Python code
  reticulate::py_run_string("
import cellrank as cr
import scanpy as sc

if 'connectivities' not in adata.obsp:
    sc.pp.neighbors(adata, random_state=0)

pk = cr.kernels.PseudotimeKernel(adata, time_key=time_key_from_r)
pk.compute_transition_matrix()
")

  invisible(NULL)
}

#' @param load.adata Path to a previously saved AnnData object (in h5ad format) which can be directly loaded to avoid re-running preprocessing. If NULL, reticulate will automatically use the existing AnnData object `adata` in the Python environment for plotting. Default: NULL.
#' @param style Style of the cellrank plot, allowing for different visual representations such as 'stream', 'grid', or 'scatter'. Default: c("stream", "grid", "scatter").
#' @param basis The embedding to be used for plotting, typically "ms", "umap" or "umap_cell_embeddings".
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
#' @rdname Cellrank-Methods
#' @export

Cellrank.Plot <- function(
    load.adata = NULL,
    basis,
    color = NULL,
    groups = NULL,
    palette = NULL,
    alpha = 0.15,
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
    adata.Load(load.adata)
  }

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

  # Filter out NULL arguments
  args_list <- args_list[!sapply(args_list, is.null)]

  # Convert to Python function call string
  args_string <- paste(paste(names(args_list), "=", args_list, sep = ""), collapse = ", ")
  py_command <- glue("pk.plot_projection({args_string})")

  # Execute the Python command
  reticulate::py_run_string(py_command)

  invisible(NULL)
}
