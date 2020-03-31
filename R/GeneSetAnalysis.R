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
                            export_to_matrix = F){
  library(Seurat)
  library(dplyr)
  library(rlang)
  library(AUCell)
  library(rlist)

  message(paste(Sys.time(), "Start filtering gene sets"))
  filter <- sapply(genesets,
                   function(x){
                     length(x) >= n.min &
                       length(x) <= n.max &
                       sum(rownames(seu) %in% x)/length(x) > (1 - ratio)
                   })
  GenesetList <- genesets[filter]
  nCores <- nCores %||% parallel::detectCores()
  if(is.null(seu@misc$AUCell[["cells_rankings"]])){
    message(paste(Sys.time(), "Build AUC Rank"))
    seu <- BuildAUCRank(seu, slot = slot, assay = assay, nCores = nCores)
  }
  message(paste(Sys.time(), "Calculating", length(GenesetList), "gene set(s)"))
  n.items.part <- 1e6 / ncol(seu) * nCores
  splited_terms <- split(GenesetList, ceiling((1:length(GenesetList))/n.items.part))
  message(paste(Sys.time(), "Split gene set(s) into", length(splited_terms), "part(s)"))
  AUC_matrix <-
    splited_terms %>%
    lapply(function(x) AUCell_calcAUC(x, seu@misc$AUCell[["cells_rankings"]], nCores = nCores) %>% getAUC()) %>%
    list.rbind() %>%
    .[apply(., 1, sum)>0, ]
  message(paste0("\n", Sys.time(), " Done"))
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
# genesets <- Genesets$mouse$GSEA_mouse_gene_transformed$`hallmark gene sets`
#
# seu <- GeneSetAnalysis(seu, genesets)
# seu@misc$AUCell$genesets[1:5,1:5]
