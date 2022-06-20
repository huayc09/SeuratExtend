# imported from AUCell package

AUCell_buildRankings <-
  new("MethodDefinition", .Data = function (
    exprMat, featureType = "genes",
    plotStats = TRUE, splitByBlocks = FALSE, BPPARAM = NULL,
    keepZeroesAsNA = FALSE, verbose = TRUE, nCores = NULL, mctype = NULL,
    ...) {
    .local <- function (
    exprMat, featureType = "genes", plotStats = TRUE,
    splitByBlocks = FALSE, BPPARAM = NULL, keepZeroesAsNA = FALSE,
    verbose = TRUE, nCores = NULL, mctype = NULL)
    {
      .AUCell_buildRankings(
        exprMat = exprMat, featureType = featureType,
        splitByBlocks = splitByBlocks, BPPARAM = BPPARAM,
        plotStats = plotStats, keepZeroesAsNA = keepZeroesAsNA,
        verbose = verbose, nCores = nCores, mctype = mctype)
    }
    .local(exprMat, featureType, plotStats, splitByBlocks, BPPARAM,
           keepZeroesAsNA, verbose, nCores, mctype, ...)
  }, target = new(
    "signature", .Data = "dgCMatrix", names = "exprMat",
    package = "Matrix"), defined = new(
      "signature", .Data = "dgCMatrix",
      names = "exprMat", package = "Matrix"), generic = "AUCell_buildRankings")

.AUCell_buildRankings <-
  function (exprMat, featureType = "genes", splitByBlocks = FALSE,
            keepZeroesAsNA = FALSE, BPPARAM = NULL, plotStats = TRUE,
            verbose = TRUE, nCores = NULL, mctype = NULL) {
    library(DelayedArray)
    if ((!splitByBlocks) && ("dgCMatrix" %in% class(exprMat)))
      stop("To use a dgCMatrix as input set 'splitByBlocks=TRUE'.")
    if (!is.null(nCores))
      warning("nCores is no longer used. It will be deprecated in the next AUCell version.")
    if (keepZeroesAsNA) {
      zeroesCoords <- which(exprMat == 0, arr.ind = TRUE)
    }
    nGenesDetected <- numeric(0)
    if (plotStats) {
      msg <- tryCatch(plotGeneCount(
        exprMat, plotStats = plotStats,
        verbose = verbose), error = function(e) {
          return(e)
        })
      if (methods::is(msg, "error")) {
        warning("There has been an error in plotGeneCount() [Message: ",
                msg$message, "]. Proceeding to calculate the rankings...",
                sep = "")
      }
      else {
        if (is.numeric(nGenesDetected))
          nGenesDetected <- msg
      }
    }
    rowNames <- rownames(exprMat)
    colNames <- colnames(exprMat)
    exprMat <- -exprMat
    if (splitByBlocks) {
      exprMat <- do.call(cbind, blockApply(
        DelayedArray(exprMat),
        FUN = DelayedMatrixStats::colRanks, ties.method = "random",
        preserveShape = TRUE, BPPARAM = BPPARAM, grid = colAutoGrid(exprMat)))
    }
    else {
      exprMat <- DelayedMatrixStats::colRanks(exprMat, ties.method = "random",
                                              preserveShape = TRUE, decreasing = TRUE)
    }
    rownames(exprMat) <- rowNames
    colnames(exprMat) <- colNames
    if (keepZeroesAsNA) {
      exprMat[which(zeroesCoords == 0, arr.ind = TRUE)] <- NA
    }
    names(dimnames(exprMat)) <- c(featureType, "cells")
    new("aucellResults", SummarizedExperiment::SummarizedExperiment(assays = list(ranking = exprMat)),
        nGenesDetected = nGenesDetected)
  }

#' @title Build AUCell Rank
#' @description Build AUCell ranking matrix using 'AUCell' package:
#'
#' Aibar et al. (2017) SCENIC: single-cell regulatory network inference and clustering.
#' Nature Methods. doi: 10.1038/nmeth.4463
#'
#' Aibar. et al. (2016) AUCell: Analysis of 'gene set' activity in single-cell RNA-seq data.
#' R/Bioconductor package.
#' @param seu Seurat object
#' @param slot Slot to pull feature data for, Default: 'counts'
#' @param assay Name of assay to use, Default: 'RNA'
#' @return Seurat object.
#' @details AUCell ranking matrix is saved in seu@misc$AUCell[["cells_rankings"]]
#' @examples
#' pbmc <- BuildAUCRank(pbmc)
#' pbmc@misc$AUCell[["cells_rankings"]]
#' @rdname BuildAUCRank
#' @export

BuildAUCRank <- function(seu, slot = "counts", assay = "RNA"){
  # import("AUCell")
  library(SeuratObject)

  # check_pkg_version("AUCell","1.18")
  seu@misc$AUCell<-list()
  seu@misc$AUCell[["cells_rankings"]] <-
    AUCell_buildRankings(
      exprMat = GetAssayData(seu, slot = slot, assay = assay),
      splitByBlocks = T, plotStats = F)
  return(seu)
}

#  Imported from AUCell package

getAUC <- function(object) {
  SummarizedExperiment::assays(object)[["AUC"]]
}

getRanking <- function(object) {
  SummarizedExperiment::assays(object)[["ranking"]]
}

AUCell_calcAUC <-
  new("MethodDefinition", .Data = function (
    geneSets, rankings,
    nCores = 1, normAUC = TRUE,
    aucMaxRank = ceiling(0.05 * nrow(rankings)), verbose = TRUE) {
    .AUCell_calcAUC(
      geneSets = geneSets, rankings = rankings,
      nCores = nCores, aucMaxRank = aucMaxRank, verbose = verbose)
  }, target = new(
    "signature", .Data = "list", names = "geneSets",
    package = "methods"), defined = new(
      "signature", .Data = "list",
      names = "geneSets", package = "methods"), generic = "AUCell_calcAUC")

.AUCell_calcAUC <- function (
    geneSets, rankings,
    nCores = 1, mctype = c("domc","snow")[1],
    normAUC = TRUE,
    aucMaxRank = ceiling(0.05 * nrow(rankings)), verbose = TRUE) {
  if (!is.list(geneSets))
    stop("geneSets should be a named list.")
  if (is.null(names(geneSets)))
    stop("geneSets should be a named list.")
  if (any(lengths(geneSets) <= 0)) {
    warning(
      "Ignoring the following empty sets: ",
      paste0(names(lengths(geneSets))[which(lengths(geneSets) <= 0)], collapse = ", "))
    geneSets <- geneSets[which(lengths(geneSets) > 0)]
  }
  if (length(geneSets) <= 0)
    stop("No geneSets provided or remaining.")
  if (nCores > length(geneSets))
    nCores <- length(geneSets)
  if ((aucMaxRank < 300) && verbose)
    warning("Using only the first ", aucMaxRank, " genes (aucMaxRank) to calculate the AUC.",
            immediate. = TRUE)
  if (aucMaxRank <= 0)
    stop("aucMaxRank should be a positive value.")
  if (!methods::is(rankings, "aucellResults") || SummarizedExperiment::assayNames(rankings) !=
      "ranking") {
    if (methods::is(rankings, "matrixWrapper")) {
      stop("These rankings were built with a previous AUCell version. ",
           "Please update them with updateAucellResults(..., objectType=\"ranking\")")
    }
    stop("Rankings should be a the object returned by AUCell_buildRankings()")
  }
  else {
    rankings <- getRanking(rankings)
  }
  if (!normAUC)
    .AUC.geneSet <- .AUC.geneSet_old()
  if (normAUC)
    .AUC.geneSet <- .AUC.geneSet_norm
  if (nCores == 1) {
    gSetName <- NULL
    aucMatrix <- sapply(names(geneSets), function(gSetName) .AUC.geneSet(
      geneSet = geneSets[[gSetName]],
      rankings = rankings, aucMaxRank = aucMaxRank, gSetName = gSetName))
    aucMatrix <- t(aucMatrix)
  }
  if (nCores > 1) {
    if (!mctype %in% c("domc"))
      stop("Valid 'mctype': 'doMC'")
    if (mctype == "snow") {
      cl <- parallel::makeCluster(nCores, type = "SOCK")
      doSNOW::registerDoSNOW(cl)
      if (verbose)
        message("Using ", length(cl), " cores with SNOW.")
      parallel::clusterExport(cl, c("geneSets", "rankings",
                                    "aucMaxRank", ".AUC.geneSet", ".auc", ".AUC.geneSet_norm",
                                    ".AUC.geneSet_old"), envir = environment())
      opts <- list(preschedule = TRUE)
      aucMatrix <- suppressWarnings(
        plyr::llply(.data = names(geneSets),
                    .fun = function(gSetName) setNames(list(
                      .AUC.geneSet(geneSet = geneSets[[gSetName]],
                                   rankings = rankings, aucMaxRank = aucMaxRank,
                                   gSetName = gSetName)), gSetName), .parallel = TRUE,
                    .paropts = list(.options.snow = opts), .inform = FALSE))
      aucMatrix <- do.call(rbind, unlist(aucMatrix, recursive = FALSE))
      parallel::stopCluster(cl)
    }
    if (mctype == "domc") {
      import("doMC")
      import("foreach")
      registerDoMC(nCores)
      if (verbose)
        message("Using ", foreach::getDoParWorkers(),
                " cores with doMC.")
      aucMatrix <- foreach::"%dopar%"(
        foreach(gSetName = names(geneSets)),
        {
          setNames(list(.AUC.geneSet(geneSet = geneSets[[gSetName]],
                                     rankings = rankings, aucMaxRank = aucMaxRank,
                                     gSetName = gSetName)), gSetName)
        })
      aucMatrix <- do.call(rbind, unlist(aucMatrix, recursive = FALSE))
    }
  }
  aucMatrix <- aucMatrix[intersect(names(geneSets), rownames(aucMatrix)),
                         , drop = FALSE]
  missingSets <- names(geneSets)[which(!names(geneSets) %in%
                                         rownames(aucMatrix))]
  if (length(missingSets) > 0)
    warning("The AUC for the following sets was not calculated: ",
            paste(missingSets, collapse = ", "))
  missingGenes <- as.matrix(aucMatrix[, c("missing", "nGenes"),
                                      drop = FALSE])
  missingPercent <- as.numeric(missingGenes[, "missing"])/as.numeric(missingGenes[,
                                                                                  "nGenes"])
  missingPercent <- setNames(missingPercent, rownames(missingGenes))
  if (all(missingPercent >= 0.8))
    stop("Fewer than 20% of the genes in the gene sets are included in the rankings.",
         "Check wether the gene IDs in the 'rankings' and 'geneSets' match.")
  if (any(missingPercent > 0.8)) {
    warning("The following gene sets will be excluded from the analysis",
            "(less than 20% of their genes are available):\n",
            paste(names(missingPercent)[which(missingPercent >=
                                                0.8)], collapse = ", "), sep = "", immediate. = TRUE)
    aucMatrix <- aucMatrix[which(missingPercent < 0.8),
                           , drop = FALSE]
  }
  missingGenes <- missingGenes[rownames(aucMatrix), , drop = FALSE]
  if (sum(missingGenes[, "missing"]) > 0) {
    msg1 <- "Genes in the gene sets NOT available in the dataset: \n"
    msg2 <- sapply(rownames(missingGenes)[
      which(missingGenes[, "missing"] > 0.01)],
      function(gSetName) {
        paste(
          "\t", gSetName, ": \t",
          missingGenes[gSetName,"missing"], " (",
          round(missingPercent[gSetName] * 100),
          "% of ", missingGenes[gSetName, "nGenes"],
          ")", sep = "")
      })
    if (verbose)
      message(msg1, paste(msg2, collapse = "\n"), sep = "")
  }
  aucMatrix <- aucMatrix[, 1:(ncol(aucMatrix) - 2), drop = FALSE]
  names(dimnames(aucMatrix)) <- c("gene sets", "cells")
  new("aucellResults", SummarizedExperiment::SummarizedExperiment(assays = list(AUC = aucMatrix)))
}

.AUC.geneSet_norm <-
  function (geneSet, rankings, aucMaxRank, gSetName = "") {
    geneSet <- unique(geneSet)
    nGenes <- length(geneSet)
    geneSet <- geneSet[which(geneSet %in% rownames(rankings))]
    missing <- nGenes - length(geneSet)
    gSetRanks <- rankings[which(rownames(rankings) %in% geneSet),
                          , drop = FALSE]
    rm(rankings)
    aucThreshold <- round(aucMaxRank)
    x_th <- 1:nrow(gSetRanks)
    x_th <- sort(x_th[x_th < aucThreshold])
    y_th <- seq_along(x_th)
    maxAUC <- sum(diff(c(x_th, aucThreshold)) * y_th)
    auc <- apply(gSetRanks, 2, .auc, aucThreshold, maxAUC)
    return(c(auc, missing = missing, nGenes = nGenes))
  }

.auc <-
  function (oneRanking, aucThreshold, maxAUC) {
    x <- unlist(oneRanking)
    x <- sort(x[x < aucThreshold])
    y <- seq_along(x)
    sum(diff(c(x, aucThreshold)) * y)/maxAUC
  }

#' @title FUNCTION_TITLE
#' @description FUNCTION_DESCRIPTION
#' @param GenesetList PARAM_DESCRIPTION
#' @param rankings PARAM_DESCRIPTION
#' @param nCores PARAM_DESCRIPTION, Default: 1
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
#' @seealso
#'  \code{\link[parallel]{detectCores}}
#'  \code{\link[rlist]{list.rbind}}
#' @rdname calcAUC_matrix
#' @export
#' @importFrom parallel detectCores
#' @importFrom rlist list.rbind

calcAUC_matrix <- function(
    GenesetList,
    rankings,
    nCores = 1,
    verbose = TRUE,
    n.items.part = NULL) {

  if(!is.list(GenesetList)) GenesetList <- list(geneset = GenesetList)
  nCores <- as.integer(nCores)
  if(!isTRUE(nCores >= 1 & nCores <= parallel::detectCores())) {
    message("nCore is set to 1")
    nCores <- 1
  }else if(nCores > length(GenesetList)) {
    nCores <- length(GenesetList)
    message("nCore is set to ", nCores)
  }
  if(!is.null(n.items.part)) {
    n.items.part <- as.integer(n.items.part)
    if(!isTRUE(n.items.part > 1 & n.items.part <= length(GenesetList))) {
      n.items.part <- NULL
    }
  }
  # sysname <- Sys.info()[["sysname"]]
  if(verbose) message(Sys.time(), " Calculating ", length(GenesetList), " gene set(s)")

  if(is.null(n.items.part) &  nCores == 1) {
    AUC_matrix <- AUCell_calcAUC(GenesetList, rankings, nCores = nCores, verbose = verbose)
    AUC_matrix <- getAUC(AUC_matrix)
  }else if(is.null(n.items.part)) {
    GenesetList2 <- split(GenesetList, ceiling(seq_along(GenesetList) * nCores / length(GenesetList)))
    import("doParallel")
    import("doRNG")
    import("foreach")
    cl <- makeCluster(nCores)
    registerDoParallel(cl)
    AUC <- foreach(i = GenesetList2) %dopar%
      getAUC(AUCell_calcAUC(i, rankings = rankings, nCores = 1, verbose = FALSE))
    stopCluster(cl)
    AUC_matrix <- rlist::list.rbind(AUC)
  }else{
    if(verbose) message(Sys.time(), " Split gene set(s) into ", n.items.part, " part(s)")
    splited_terms <- split(GenesetList, ceiling(seq_along(GenesetList) * n.items.part / length(GenesetList)))
    for (i in names(splited_terms)) {
      splited_terms[[i]] <- calcAUC_matrix(splited_terms[[i]], rankings = rankings, nCores = nCores, verbose = FALSE)
    }
    AUC_matrix <- rlist::list.rbind(splited_terms)
  }
  if(verbose) message(Sys.time(), " Done")
  return(AUC_matrix)
}


