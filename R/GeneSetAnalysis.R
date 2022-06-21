#' @title FUNCTION_TITLE
#' @description FUNCTION_DESCRIPTION
#' @param seu PARAM_DESCRIPTION, Default: NULL
#' @param genesets PARAM_DESCRIPTION
#' @param title PARAM_DESCRIPTION, Default: 'genesets'
#' @param ratio PARAM_DESCRIPTION, Default: 0.4
#' @param n.min PARAM_DESCRIPTION, Default: 1
#' @param n.max PARAM_DESCRIPTION, Default: Inf
#' @param slot PARAM_DESCRIPTION, Default: 'counts'
#' @param assay PARAM_DESCRIPTION, Default: 'RNA'
#' @param nCores PARAM_DESCRIPTION, Default: 1
#' @param aucMaxRank PARAM_DESCRIPTION, Default: NULL
#' @param export_to_matrix PARAM_DESCRIPTION, Default: F
#' @param verbose PARAM_DESCRIPTION, Default: TRUE
#' @param n.items.part PARAM_DESCRIPTION, Default: NULL
#' @return OUTPUT_DESCRIPTION
#' @details DETAILS
#' @examples
#' \dontrun{
#' if(interactive()){
#'  #EXAMPLE1
#'  }
#' }
#' @rdname GeneSetAnalysis
#' @export

GeneSetAnalysis <- function(
    seu = NULL,
    genesets,
    title = "genesets",
    ratio = 0.4, n.min = 1, n.max = Inf,
    slot = "counts", assay = "RNA", nCores = 1,
    aucMaxRank = NULL,
    export_to_matrix = F, verbose = TRUE,
    n.items.part = NULL) {

  library(SeuratObject)
  DefaultAssay(seu) <- assay
  genesets <- FilterGenesets(
    genes_in_data = rownames(seu),
    genesets = genesets,
    ratio = ratio,
    n.min = n.min,
    n.max = n.max,
    verbose = verbose)

  seu <- BuildAUCRank(
    seu,
    slot = slot,
    assay = assay,
    verbose = verbose)

  if(is.null(aucMaxRank)) aucMaxRank <- ceiling(0.05 * nrow(seu))
  AUC_matrix <- calcAUC_matrix(
    GenesetList = genesets,
    rankings = seu@misc$AUCell$cells_rankings,
    nCores = nCores,
    aucMaxRank = aucMaxRank,
    verbose = verbose,
    n.items.part = n.items.part)

  if(export_to_matrix) return(AUC_matrix)
  seu@misc[["AUCell"]][[title]] <- AUC_matrix
  return(seu)
}

# Internal ----------------------------------------------------------------

FilterGenesets <- function(
    genes_in_data,
    genesets,
    ratio = 0.4,
    n.min = 1,
    n.max = Inf,
    verbose = TRUE) {
  if(verbose) {
    message(paste0(
      Sys.time(),
      " Start filtering ",length(genesets)," gene set(s): "),
      "n(Genes) >= ", n.min, ", n(Genes) <= ", n.max,
      ", at least ", ratio * 100, "% of genes found in the datasets")
  }
  filter <- sapply(
    genesets,
    function(x){
      length(x) >= n.min &
        length(x) <= n.max &
        sum(genes_in_data %in% x)/length(x) > ratio
    })
  if(!any(filter)) stop("No gene set(s) pass the filter. Please check inputs.")
  genesets <- genesets[filter]
  if(verbose) {
    message(paste(
      Sys.time(),
      length(genesets),"gene set(s) passed the filter"))
  }
  return(genesets)
}
