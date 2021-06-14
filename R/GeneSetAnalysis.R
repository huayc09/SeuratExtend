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
#' @param nCores PARAM_DESCRIPTION, Default: getOption("nCores")
#' @param export_to_matrix PARAM_DESCRIPTION, Default: F
#' @param verbose PARAM_DESCRIPTION, Default: TRUE
#' @param n.items.part PARAM_DESCRIPTION, Default: 5e+05/ncol(seu) * parallel::detectCores()
#' @return OUTPUT_DESCRIPTION
#' @details DETAILS
#' @examples
#' \dontrun{
#' if(interactive()){
#'  #EXAMPLE1
#'  }
#' }
#' @seealso
#'  \code{\link[parallel]{detectCores}}
#' @rdname GeneSetAnalysis
#' @export
#' @importFrom parallel detectCores

GeneSetAnalysis <- function(seu = NULL, genesets, title = "genesets",
                            ratio = 0.4, n.min = 1, n.max = Inf,
                            slot = "counts", assay = "RNA", nCores = getOption("nCores"),
                            export_to_matrix = F, verbose = TRUE,
                            n.items.part = 5e5 / ncol(seu) * parallel::detectCores()){
  library(Seurat)
  library(dplyr)
  library(rlang)
  library(AUCell)
  library(rlist)

  message(paste(Sys.time(), "Start filtering gene sets"))
  DefaultAssay(seu) <- assay
  filter <- sapply(genesets,
                   function(x){
                     length(x) >= n.min &
                       length(x) <= n.max &
                       sum(rownames(seu) %in% x)/length(x) > (1 - ratio)
                   })
  GenesetList <- genesets[filter]
  nCores <- nCores %||% parallel::detectCores()
  if(is.null(seu@misc$AUCell[["cells_rankings"]])){
    if(verbose) message(Sys.time(), " Build AUC Rank")
    seu <- BuildAUCRank(seu, slot = slot, assay = assay, nCores = nCores)
  } else if(!identical(colnames(seu), colnames(seu@misc$AUCell$cells_rankings))) {
    if(verbose) message(Sys.time(), " Pre-existing cell ranking matrix has different cell IDs with current seurat object. ",
                        "Re-build AUC Rank")
    seu <- BuildAUCRank(seu, slot = slot, assay = assay, nCores = nCores)
  }
  AUC_matrix <- calcAUC_matrix(GenesetList, rankings = seu@misc$AUCell$cells_rankings,
                               nCores = nCores, n.items.part = n.items.part, verbose = verbose)
  if(export_to_matrix) return(AUC_matrix)
  seu@misc[["AUCell"]][[title]] <- AUC_matrix
  return(seu)
}

# seu <- readRDS("~/R documents/2020-2-10 EC PyMT and E0771/rds/PyMTEC_old.rds")
# ratio = 0.4
# n.min = 1
# n.max = Inf
# title = "genesets"
# slot = "counts"
# assay = "RNA"
# nCores = getOption("nCores")
# export_to_matrix = F
# genesets <- Genesets_data$mouse$GSEA_mouse_gene_transformed$`hallmark gene sets`
#
# seu <- GeneSetAnalysis(seu, genesets, nCores = 12)
# seu@misc$AUCell$genesets[1:5,1:5]
