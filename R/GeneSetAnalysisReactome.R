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
           slot = "counts", assay = "RNA", nCores = getOption("nCores")){
    check_spe(spe)
    DatabaseList <- Reactome_Data[[spe]]$Roots
    library(Seurat)
    library(dplyr)
    library(rlang)

    if(is.null(seu)){
      return(DatabaseList)
    }
    if(parent == "All"){
      GenesetNames <- names(Reactome_Data[[spe]]$Path2Gene)
    }else if(all(parent %in% DatabaseList)){
      GenesetNames <- GetAllChilrenReactome(names(DatabaseList)[DatabaseList %in% parent], spe = spe)
    }else if(all(parent %in% names(DatabaseList))){
      GenesetNames <- GetAllChilrenReactome(parent, spe = spe)
    }else{
      print(DatabaseList)
      return(seu)
    }
    Time1 <- Sys.time()
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
      seu <- BuildAUCRank(seu, slot = slot, assay = assay, nCores = nCores)
    }
    seu@misc[["AUCell"]][["Reactome"]][[parent]] <-
      AUCell_calcAUC(GenesetList, seu@misc$AUCell[["cells_rankings"]], nCores = nCores) %>%
      getAUC() %>%
      .[apply(., 1, sum)>0, ]
    Time2 <- Sys.time()
    message <- paste0("\nTime usage: ", round(difftime(Time2, Time1, units='mins'), digits = 2), " mins")
    message(message)
    return(seu)
  }

# options(spe = "mouse", nCores = 12)
# spe = "mouse"
# parent = "All"
# ratio = 0.4
# n.min = 1
# n.max = Inf
# only.end.terms = F
# slot = "counts"
# assay = "RNA"
# nCores = 12
# seu <- readRDS("~/R documents/2020-2-10 EC PyMT and E0771/rds/PyMTEC_old.rds")
# GeneSetAnalysisReactome(seu)
