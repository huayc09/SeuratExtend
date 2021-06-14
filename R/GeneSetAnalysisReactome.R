#' @title FUNCTION_TITLE
#' @description FUNCTION_DESCRIPTION
#' @param seu PARAM_DESCRIPTION, Default: NULL
#' @param parent PARAM_DESCRIPTION, Default: 'All'
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
#' @rdname GeneSetAnalysisReactome
#' @export
#' @importFrom parallel detectCores

GeneSetAnalysisReactome <-
  function(seu = NULL, parent = "All", spe = getOption("spe"),
           ratio = 0.4, n.min = 1, n.max = Inf, only.end.terms = F,
           slot = "counts", assay = "RNA", nCores = getOption("nCores"),
           export_to_matrix = F, verbose = TRUE,
           n.items.part = 5e5 / ncol(seu) * parallel::detectCores()){
    check_spe(spe)
    DatabaseList <- Reactome_Data[[spe]]$Roots
    library(Seurat)
    library(dplyr)
    library(rlang)
    library(AUCell)
    library(rlist)

    if(is.null(seu)){
      return(DatabaseList)
    }
    if(parent == "All"){
      GenesetNames <- names(Reactome_Data[[spe]]$Path2Gene)
    }else if(all(parent %in% DatabaseList)){
      GenesetNames <- GetAllChilrenReactome(names(DatabaseList)[DatabaseList %in% parent], spe = spe)
    }else if(all(parent %in% names(Reactome_Data[[spe]]$Path2Gene))){
      GenesetNames <- GetAllChilrenReactome(parent, spe = spe)
    }else{
      print(DatabaseList)
      return(seu)
    }
    message(paste(Sys.time(), "Start filtering gene sets"))
    DefaultAssay(seu) <- assay
    filter <- sapply(Reactome_Data[[spe]]$Path2Gene[GenesetNames],
                     function(x){
                       length(x) >= n.min &
                         length(x) <= n.max &
                         sum(rownames(seu) %in% x)/length(x) > (1 - ratio)
                     })
    GenesetNames <- GenesetNames[filter]
    if(only.end.terms) GenesetNames <- GetEndTermsReactome(GenesetNames, spe)
    GenesetList <- Reactome_Data[[spe]]$Path2Gene[GenesetNames]
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
    seu@misc[["AUCell"]][["Reactome"]][[parent]] <- AUC_matrix
    return(seu)
  }

# options(spe = "mouse", nCores = 12)
# parent = "All"
# ratio = 0.4
# n.min = 1
# n.max = Inf
# only.end.terms = F
# slot = "counts"
# assay = "RNA"
# nCores = 12
# seu <- readRDS("~/R documents/2020-2-10 EC PyMT and E0771/rds/PyMTEC_old.rds")
# seu <- GeneSetAnalysisReactome(seu)
