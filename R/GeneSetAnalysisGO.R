#' @title FUNCTION_TITLE
#' @description FUNCTION_DESCRIPTION
#' @param seu PARAM_DESCRIPTION, Default: NULL
#' @param dataset PARAM_DESCRIPTION, Default: 'BP'
#' @param root PARAM_DESCRIPTION, Default: 'BP'
#' @param spe PARAM_DESCRIPTION, Default: getOption("spe")
#' @param ratio PARAM_DESCRIPTION, Default: 0.4
#' @param n.min PARAM_DESCRIPTION, Default: 1
#' @param n.max PARAM_DESCRIPTION, Default: Inf
#' @param only.end.terms PARAM_DESCRIPTION, Default: F
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
#' @rdname GeneSetAnalysisGO
#' @export
#' @importFrom parallel detectCores

GeneSetAnalysisGO<-function(seu = NULL, dataset = "BP", root = "BP", spe = getOption("spe"),
                            ratio = 0.4, n.min = 1, n.max = Inf, only.end.terms = F,
                            slot = "counts", assay = "RNA", nCores = getOption("nCores"),
                            export_to_matrix = F){
  check_spe(spe)
  DatabaseList<-list("BP"=c("BP" = "GO:0008150",
                            "immune_system_process" = "GO:0002376",
                            "response_to_stimulus" = "GO:0050896",
                            "signaling" = "GO:0023052",
                            "metabolic_process" = "GO:0008152",
                            "regulation_of_vasculature_development" = "GO:1901342",
                            "signal_transduction" = "GO:0007165"),
                     "MF"=c("MF" = "GO:0003674"),
                     "CC"=c("CC" = "GO:0005575"))
  library(Seurat)
  library(dplyr)
  library(rlang)
  library(AUCell)
  library(rlist)

  if(is.null(seu)){
    return(DatabaseList)
  }
  if(all(dataset %in% names(DatabaseList[[root]]))){
    GenesetNames <- GetAllChilrenGO(DatabaseList[[root]][dataset], spe = spe)
  }else{
    AllGeneset <- GetAllChilrenGO(DatabaseList[[root]][dataset], spe = spe)
    if(all(dataset %in% AllGeneset)){
      GenesetNames <- GetAllChilrenGO(dataset, spe = spe)
    }else{
      print(DatabaseList)
      return(seu)
    }
  }
  message(paste(Sys.time(), "Start filtering gene sets"))
  filter <- sapply(GO_Data[[spe]]$GO2Gene[GenesetNames],
                   function(x){
                     length(x) >= n.min &
                       length(x) <= n.max &
                       sum(rownames(seu) %in% x)/length(x) > (1 - ratio)
                   })
  GenesetNames <- GenesetNames[filter]
  if(only.end.terms) GenesetNames <- GetEndTermsGO(GenesetNames, spe)
  GenesetList <- GO_Data[[spe]]$GO2Gene[GenesetNames]
  nCores <- nCores %||% parallel::detectCores()
  if(is.null(seu@misc$AUCell[["cells_rankings"]])){
    message(paste(Sys.time(), "Build AUC Rank"))
    seu <- BuildAUCRank(seu, slot = slot, assay = assay, nCores = nCores)
  }
  message(paste(Sys.time(), "Calculating", length(GenesetList), "gene set(s)"))
  n.items.part <- 2e6 / ncol(seu) * nCores
  splited_terms <- split(GenesetList, ceiling((1:length(GenesetList))/n.items.part))
  message(paste(Sys.time(), "Split gene set(s) into", length(splited_terms), "part(s)"))
  AUC_matrix <-
    splited_terms %>%
    lapply(function(x) AUCell_calcAUC(x, seu@misc$AUCell[["cells_rankings"]], nCores = nCores) %>% getAUC()) %>%
    list.rbind() %>%
    .[apply(., 1, sum)>0, ]
  message(paste0("\n", Sys.time(), " Done"))
  if(export_to_matrix) return(AUC_matrix)
  seu@misc[["AUCell"]][["GO"]][[dataset]] <- AUC_matrix
  return(seu)
}

# dataset = "BP"
# root = "BP"
# spe = getOption("spe")
# ratio = 0.4
# n.min = 1
# n.max = Inf
# only.end.terms = F
# slot = "counts"
# assay = "RNA"
# nCores = getOption("nCores")
# export_to_matrix = F
# seu <- GeneSetAnalysisGO(seu)

