#' @title FUNCTION_TITLE
#' @description FUNCTION_DESCRIPTION
#' @param seu PARAM_DESCRIPTION, Default: NULL
#' @param dataset PARAM_DESCRIPTION, Default: 'BP'
#' @param root PARAM_DESCRIPTION, Default: 'BP'
#' @param spe PARAM_DESCRIPTION, Default: 'mouse'
#' @param ratio PARAM_DESCRIPTION, Default: 0.4
#' @param n.min PARAM_DESCRIPTION, Default: 1
#' @param n.max PARAM_DESCRIPTION, Default: Inf
#' @param only.end.terms PARAM_DESCRIPTION, Default: F
#' @param slot PARAM_DESCRIPTION, Default: 'counts'
#' @param assay PARAM_DESCRIPTION, Default: 'RNA'
#' @param nCores PARAM_DESCRIPTION, Default: 2
#' @return OUTPUT_DESCRIPTION
#' @details DETAILS
#' @examples
#' \dontrun{
#' if(interactive()){
#'  #EXAMPLE1
#'  }
#' }
#' @rdname GeneSetAnalysisGO
#' @export
GeneSetAnalysisGO<-function(seu = NULL, dataset = "BP", root = "BP", spe = "mouse",
                            ratio = 0.4, n.min = 1, n.max = Inf, only.end.terms = F,
                            slot = "counts", assay = "RNA", nCores = 2){
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

  if(is.null(seu)){
    return(DatabaseList)
  }
  if(all(dataset %in% names(DatabaseList[[root]]))){
    GenesetNames <- GetAllChilrenGO(DatabaseList[[root]][dataset], spe = spe)
  }else{
    AllGeneset <- GetAllChilrenGO(DatabaseList[[root]]["All"], spe = spe)
    if(all(dataset %in% AllGeneset)){
      GenesetNames <- GetAllChilrenGO(dataset, spe = spe)
    }else{
      print(DatabaseList)
      return(seu)
    }
  }
  Time1 <- Sys.time()
  filter <- lapply(GO_Data[[spe]]$GO2Gene[GenesetNames],
                   function(x){
                     length(x) >= n.min &
                       length(x) <= n.max &
                       sum(rownames(seu) %in% x)/length(x) > (1 - ratio)
                   })
  GenesetNames <- GenesetNames[unlist(filter)]
  if(only.end.terms) GenesetNames <- GetEndTermsGO(GenesetNames, spe)
  GenesetList <- GO_Data[[spe]]$GO2Gene[GenesetNames]
  if(is.null(seu@misc$AUCell[["cells_rankings"]])){
    seu <- BuildAUCRank(seu, slot = slot, assay = assay)
  }
  seu@misc[["AUCell"]][["GO"]][[dataset]] <-
    AUCell_calcAUC(GenesetList, seu@misc$AUCell[["cells_rankings"]], nCores = nCores) %>%
    getAUC() %>%
    .[apply(., 1, sum)>0, ]
  Time2 <- Sys.time()
  message <- paste0("\nTime usage: ", round(difftime(Time2, Time1, units='mins'), digits = 2), " mins")
  message(message)
  return(seu)
}

# library(magrittr)
# seu <- GeneSetAnalysisGO(seu)
# ScoreAndOrder(seu@misc[["AUCell"]][["GO"]][["All"]], f, n = 7) %>%
#   set_rownames(RenameGO(rownames(.))) %>%
#   Heatmap()
