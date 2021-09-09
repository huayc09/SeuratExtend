#' @title FUNCTION_TITLE
#' @description FUNCTION_DESCRIPTION
#' @param seu PARAM_DESCRIPTION
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
#' @rdname BuildAUCRank
#' @export

BuildAUCRank <- function(seu, slot = "counts", assay = "RNA", nCores = getOption("nCores")){
  library(AUCell)
  library(Seurat)
  library(rlang)
  nCores <- nCores %||% parallel::detectCores()
  matr <- GetAssayData(seu, slot = slot, assay = assay)
  seu@misc$AUCell<-list()
  seu@misc$AUCell[["cells_rankings"]] <- AUCell_buildRankings(matr, plotStats=TRUE, nCores = nCores)
  return(seu)
}

#' @title FUNCTION_TITLE
#' @description FUNCTION_DESCRIPTION
#' @param GenesetList PARAM_DESCRIPTION
#' @param rankings PARAM_DESCRIPTION
#' @param nCores PARAM_DESCRIPTION, Default: NULL
#' @param verbose PARAM_DESCRIPTION, Default: TRUE
#' @param n.items.part PARAM_DESCRIPTION, Default: 5e+05/ncol(rankings) * parallel::detectCores()
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
#'  \code{\link[AUCell]{AUCell_calcAUC}}
#' @rdname calcAUC_matrix
#' @export
#' @importFrom parallel detectCores

calcAUC_matrix <-
  function(GenesetList, rankings, nCores = NULL, verbose = TRUE,
           n.items.part = 5e5 / ncol(rankings) * parallel::detectCores()) {
    import("AUCell")
    library(rlang)
    library(parallel)
    library(rlist)

    nCores <- nCores %||% detectCores()
    if(!is.list(GenesetList)) GenesetList <- list(geneset = GenesetList)
    sysname <- Sys.info()[["sysname"]]
    if(verbose) message(Sys.time(), " Calculating ", length(GenesetList), " gene set(s)")
    splited_terms <- split(GenesetList, ceiling(seq_along(GenesetList)/n.items.part))
    if(verbose) message(Sys.time(), " Split gene set(s) into ", length(splited_terms), " part(s)")
    if(sysname != "Windows" | nCores == 1) {
      AUC_matrix <-
        splited_terms %>%
        lapply(function(x){
          AUCell_calcAUC(x, rankings, nCores = nCores, verbose = verbose) %>%
            getAUC()
        }) %>%
        list.rbind() %>%
        .[apply(., 1, sum)>0, ]
    } else {
      AUCell_calcAUC2 <- function(x, rankings, nCores, verbose) {
        import("doParallel")
        import("doRNG")
        library(foreach)
        cl <- makeCluster(nCores)
        registerDoParallel(cl)
        AUC <- foreach(i = split(x, ceiling(seq_along(x) * nCores / length(x)))) %dopar%
          AUCell::AUCell_calcAUC(i, rankings = rankings, nCores = 1, verbose = verbose)
        stopCluster(cl)
        AUC <- lapply(AUC, getAUC) %>%
          list.rbind() %>%
          unique()
        return(AUC)
      }
      AUC <- list()
      for (i in seq_along(splited_terms)) {
        AUC[[i]] <- AUCell_calcAUC2(splited_terms[[i]], rankings = rankings, nCores = nCores, verbose = verbose)
      }
      AUC_matrix <- AUC %>%
        list.rbind() %>%
        .[apply(., 1, sum)>0, ]
    }
    message("\n", Sys.time(), " Done")
    return(AUC_matrix)
  }
