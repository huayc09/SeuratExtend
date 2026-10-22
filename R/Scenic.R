#' @title Import SCENIC Loom Files into Seurat
#' @description Imports SCENIC-generated loom files into Seurat objects for further analysis. This function allows the integration of gene regulatory network insights directly into the Seurat environment. If a Seurat object is specified, results are stored in `seu@misc$SCENIC` and a new 'TF' assay is created; if not, a list containing Regulons and RegulonsAUC is returned.
#' @param loom.path Path to the SCENIC-generated loom file.
#' @param seu Optional Seurat object. If specified, the SCENIC data (regulons and their activities) are imported directly into the object, facilitating further analyses. If not provided, the function returns a list containing the regulons and their activity matrices. Default: NULL.
#' @return If a Seurat object is provided, the function returns the modified Seurat object with SCENIC data integrated. If no Seurat object is provided, a list with two elements is returned: `Regulons` containing the regulons and their gene lists, and `RegulonsAUC` containing the activity matrix of these regulons.
#' @details This function is designed to ease the integration of SCENIC analysis into Seurat workflows, allowing users to explore and visualize transcription factor activities and their regulatory impacts on gene expression within their familiar Seurat analysis pipeline. The integration involves adding a new 'TF' assay for the regulon activity and storing detailed regulon information in the object's `misc` slot.
#' @examples
#' # Example of importing a SCENIC loom file with an existing Seurat object
#' library(SeuratExtend)
#' scenic_loom_path <- file.path(tempdir(), "pyscenic_integrated-output.loom")
#' \dontrun{
#' download.file("https://zenodo.org/records/10944066/files/pbmc3k_small_pyscenic_integrated-output.loom", scenic_loom_path)
#' pbmc <- ImportPyscenicLoom(scenic_loom_path, seu = pbmc)
#'
#' # Example without an existing Seurat object
#' scenic_output <- ImportPyscenicLoom(scenic_loom_path)
#'
#' # Accessing the SCENIC data in a Seurat object
#' tf_auc <- pbmc@misc$SCENIC$RegulonsAUC
#' head(tf_auc, 3:4)
#' tf_gene_list <- pbmc@misc$SCENIC$Regulons
#' head(tf_gene_list, 5)
#' }
#' @rdname ImportPyscenicLoom
#' @export

ImportPyscenicLoom <- function(loom.path, seu = NULL) {
  library(Seurat)

  lfile <- hdf5r::H5File$new(loom.path, mode = "r")
  on.exit(lfile$close_all(), add = TRUE)

  # pySCENIC stores the regulon activity and membership matrices as 2D
  # compound datasets under col_attrs/row_attrs. hdf5r reads them as
  # data frames with cells/genes as rows and regulons as columns, which is
  # the layout the rest of this function expects.
  required.datasets <- c("col_attrs/RegulonsAUC", "row_attrs/Regulons")
  missing.datasets <- required.datasets[!vapply(required.datasets, function(path) {
    parts <- strsplit(path, "/", fixed = TRUE)[[1]]
    lfile$exists(parts[1]) && lfile[[parts[1]]]$exists(parts[2])
  }, logical(1))]
  if(length(missing.datasets) > 0) {
    stop("Loom file is missing the pySCENIC dataset(s): ",
         paste(missing.datasets, collapse = ", "),
         " (is this a pySCENIC integrated-output loom file?)")
  }

  CellID <- lfile[["col_attrs"]][["CellID"]][]
  if(is.null(seu)) {} else if(length(CellID) != length(colnames(seu))) {
    stop("Loom file and Seurat object have different cell numbers")
  }else if(!identical(CellID, colnames(seu))){
    warning("Loom file and Seurat object have the same cell numbers but different cell names")
  }

  RegulonsAUC <- lfile[["col_attrs"]][["RegulonsAUC"]][]
  colnames(RegulonsAUC) <- sub("_(+)", "", colnames(RegulonsAUC), fixed = T)
  rownames(RegulonsAUC) <- CellID

  Regulons <- lfile[["row_attrs"]][["Regulons"]][]
  colnames(Regulons) <- sub("_(+)", "", colnames(Regulons), fixed = T)
  genes <- lfile[["row_attrs"]][["Gene"]][]
  Regulons <- apply(Regulons, 2, function(x) genes[which(x==1)])

  results <- list(RegulonsAUC = RegulonsAUC,
                  Regulons = Regulons)

  if(is.null(seu)) {
    return(results)
  } else {
    seu@misc[["SCENIC"]] <- results
    seu[["TF"]] <- CreateAssayObject(data = t(RegulonsAUC))
    return(seu)
  }
}
