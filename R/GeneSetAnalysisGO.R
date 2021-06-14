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
#' @rdname GeneSetAnalysisGO
#' @export
#' @importFrom parallel detectCores

GeneSetAnalysisGO<-function(seu = NULL, dataset = "BP", root = "BP", spe = getOption("spe"),
                            ratio = 0.4, n.min = 1, n.max = Inf, only.end.terms = F,
                            slot = "counts", assay = "RNA", nCores = getOption("nCores"),
                            export_to_matrix = F, verbose = TRUE,
                            n.items.part = 5e5 / ncol(seu) * parallel::detectCores()){
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
    AllGeneset <- GetAllChilrenGO(DatabaseList[[root]][root], spe = spe)
    if(all(dataset %in% AllGeneset)){
      GenesetNames <- GetAllChilrenGO(dataset, spe = spe)
    }else{
      print(DatabaseList)
      return(seu)
    }
  }
  message(paste(Sys.time(), "Start filtering gene sets"))
  DefaultAssay(seu) <- assay
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
    if(verbose) message(paste(Sys.time(), "Build AUC Rank"))
    seu <- BuildAUCRank(seu, slot = slot, assay = assay, nCores = nCores)
  } else if(!identical(colnames(seu), colnames(seu@misc$AUCell$cells_rankings))) {
    if(verbose) message(Sys.time(), " Pre-existing cell ranking matrix has different cell IDs with current seurat object. ",
            "Re-build AUC Rank")
    seu <- BuildAUCRank(seu, slot = slot, assay = assay, nCores = nCores)
  }
  AUC_matrix <- calcAUC_matrix(GenesetList, rankings = seu@misc$AUCell$cells_rankings,
                               nCores = nCores, n.items.part = n.items.part, verbose = verbose)
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
