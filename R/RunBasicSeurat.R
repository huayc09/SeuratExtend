#' @title Run Standard Seurat Pipeline
#' @description This function processes a Seurat object through various steps including normalization, PCA, and clustering based on specified parameters. It allows for conditional execution of each step based on prior executions and parameter changes.
#' @param seu A Seurat object that will be processed.
#' @param spe Either "human" or "mouse", used in calculating percent.mt. Default: NULL.
#' @param nFeature_RNA.min The minimum number of RNA features to include in the subset. Default: 500.
#' @param nFeature_RNA.max The maximum number of RNA features to include in the subset. Default: Inf.
#' @param percent.mt.max The maximum percentage of mitochondrial content allowed for cells in the subset. Default: 20.
#' @param dims The dimensions to be used in RunPCA and FindNeighbors. Default: 1:10.
#' @param resolution The resolution parameter for the clustering algorithm. Default: 0.5.
#' @param reduction The dimensional reduction technique to be used ('pca' or 'harmony'). If NULL, the method is determined based on available data in the Seurat object. Default: NULL.
#' @param harmony.by A string representing a metadata column in the Seurat object that categorizes data, which might need transformation for numeric-only issues before running Harmony. Default: NULL.
#' @param force.Normalize Forces the normalization step regardless of previous command history. Default: FALSE.
#' @param force.RunPCA Forces the PCA step regardless of previous command history. Default: FALSE.
#' @param vars.to.regress Variables to regress out during scaling. Default: NULL.
#' @param force.RunHarmony Forces the Harmony integration step regardless of previous command history or detection. Default: FALSE.
#' @param force.FindNeighbors Forces the neighbor finding step regardless of previous command history. Default: FALSE.
#' @param force.FindClusters Forces the clustering step regardless of previous command history. Default: FALSE.
#' @param force.RunUMAP Forces the UMAP step regardless of previous command history. Default: FALSE.
#' @return Returns the modified Seurat object after performing the requested operations and updates.
#' @details This function provides flexible and conditional execution of various data processing steps in a Seurat pipeline. Depending on the parameters set, the function intelligently determines whether to rerun specific steps based on changes in parameters or previous runs recorded in the Seurat object's command history. This allows for efficient reuse of existing results while ensuring updates are applied when necessary.
#' @examples
#' library(SeuratExtend)
#' pbmc <- RunBasicSeurat(pbmc, force.Normalize = TRUE)
#' DimPlot2(pbmc, group.by = "cluster")
#' @rdname RunBasicSeurat
#' @export

RunBasicSeurat <- function(
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
) {
  library(Seurat)
  if(is.null(seu@meta.data$percent.mt)) {
    if (is.null(spe) || length(spe) != 1 || !spe %in% c("human", "mouse")) {
      stop("Invalid or missing species. Please specify 'human' or 'mouse'.")
    }
    message("Calculating mitochondria percentage 'percent.mt' for ", spe, " genes.")
    mt.pattern <- switch (
      spe,
      mouse = "^mt-",
      human = "^MT-"
    )
    seu[["percent.mt"]] <- PercentageFeatureSet(seu, pattern = mt.pattern)
  }
  seu <- subset(
    seu,
    subset = nFeature_RNA > nFeature_RNA.min &
      nFeature_RNA < nFeature_RNA.max &
      percent.mt < percent.mt.max)

  # Normalization, RNA
  normalizeOrSCTCmdExists <- any(grepl(
    "^(NormalizeData\\.|SCTransform\\.)",
    names(seu@commands))
    )
  run_norm <- (!normalizeOrSCTCmdExists || force.Normalize)
  # Check the condition to run NormalizeData
  if (run_norm) {
    # Check and set default assay to 'RNA' if necessary
    if (DefaultAssay(seu) != "RNA") {
      DefaultAssay(seu) <- "RNA"
      message("Default assay set to 'RNA'.")
    }

    # Run NormalizeData
    seu <- NormalizeData(seu)
  }

  # Check if "RunPCA." command is missing, or if forced, or if normalization was just run
  pcaCmdExists <- any(grepl("^RunPCA\\.", names(seu@commands)))
  run_pca <- (!pcaCmdExists || force.RunPCA || run_norm)
  if (run_pca) {
    # Run the steps leading up to and including PCA
    seu <- FindVariableFeatures(seu)
    seu <- ScaleData(seu, vars.to.regress = vars.to.regress)
    seu <- RunPCA(seu)
  }

  run_harmony <- FALSE
  # Check if harmony.by is not NULL and whether to force running Harmony
  if (!is.null(harmony.by)) {
    # Check if Harmony should be run or if it has been run already
    if (force.RunHarmony || !"harmony" %in% names(seu@reductions) || run_pca) {
      import("harmony")
      # Determine if the harmony.by column is numeric or a factor of numeric values
      isnum_harmonyby <- (
        is.numeric(seu@meta.data[[harmony.by]]) || (
          is.factor(seu@meta.data[[harmony.by]]) &&
            all(is.numeric(levels(seu@meta.data[[harmony.by]])))
        )
      )
      if (isnum_harmonyby) {
        # Create a new column 'harmony_safe' with 'sample_' prefix to handle numeric issues
        seu@meta.data[["harmony_safe"]] <- paste0("sample_", seu@meta.data[[harmony.by]])
        harmony.by <- "harmony_safe"
      }
      seu <- RunHarmony(seu, group.by.vars = harmony.by)
      run_harmony <- TRUE
    } else {
      message("Harmony dimension reduction has already been detected in the dataset. If you wish to rerun Harmony regardless, please set 'force.RunHarmony = TRUE'.")
    }
  }

  # Determine the reduction method
  if (is.null(reduction)) {
    reduction <- if ("harmony" %in% names(seu@reductions)) "harmony" else "pca"
    message(sprintf("Using '%s' as the reduction method for FindNeighbors, FindClusters, and RunUMAP.", reduction))
  } else {
    message(sprintf("Reduction method set by user: Using '%s' for FindNeighbors, FindClusters, and RunUMAP.", reduction))
  }

  # Find neighbors, clusters, and run UMAP based on the chosen reduction
  defaultAssay <- DefaultAssay(seu)
  commandKey <- function(cmd) paste(cmd, defaultAssay, reduction, sep = ".")

  # Check if FindNeighbors needs to be run
  fnKey <- commandKey("FindNeighbors")
  graph.name <- paste0(defaultAssay, "_snn")
  run_findneighb <- (
    force.FindNeighbors ||
      run_pca || run_harmony ||
      !fnKey %in% names(seu@commands) ||
      !identical(seu@commands[[fnKey]]$dims, dims) ||
      !graph.name %in% names(seu@graphs)
  )
  if (run_findneighb) {
    seu <- FindNeighbors(seu, dims = dims, reduction = reduction)
  } else {
    message(sprintf("FindNeighbors with '%s' reduction and already up to date; not rerun.", reduction))
  }

  # Check if FindClusters needs to be run
  run_findclus <- (
    force.FindNeighbors ||
      run_findneighb ||
      !"FindClusters" %in% names(seu@commands) ||
      seu@commands[["FindClusters"]]$resolution != resolution
  )
  if (run_findclus) {
    seu <- FindClusters(seu, resolution = resolution)
  } else {
    message(sprintf("FindClusters with resolution %s already set; not rerun.", resolution))
  }

  # Check if RunUMAP needs to be run
  ruKey <- commandKey("RunUMAP")
  run_umap <- (
    force.RunUMAP ||
      run_pca || run_harmony ||
      !ruKey %in% names(seu@commands) ||
      !identical(seu@commands[[ruKey]]$dims, dims)
  )
  if (run_umap) {
    seu <- RunUMAP(seu, dims = dims, reduction = reduction)
  } else {
    message(sprintf("RunUMAP with '%s' reduction and dims already up to date; not rerun.", reduction))
  }
  return(seu)
}
