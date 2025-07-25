# modified from AUCell package version 1.19

.AUCell_buildRankings <-
  function (exprMat, featureType = "genes",
            keepZeroesAsNA = FALSE, BPPARAM = NULL, plotStats = FALSE,
            verbose = TRUE) {
    import("DelayedArray")
    import("DelayedMatrixStats")
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

    # Check matrixStats version for useNames parameter compatibility
    matrixStats_version <- utils::packageVersion("matrixStats")
    if (matrixStats_version >= "0.60.0") {
      # For matrixStats >= 0.60.0, useNames parameter is available
      exprMat <- do.call(cbind, blockApply(
        DelayedArray(exprMat),
        FUN = colRanks, ties.method = "random", useNames = FALSE,
        preserveShape = TRUE, BPPARAM = BPPARAM, grid = colAutoGrid(exprMat)))
    } else {
      # For older versions, useNames parameter doesn't exist
      exprMat <- do.call(cbind, blockApply(
        DelayedArray(exprMat),
        FUN = colRanks, ties.method = "random",
        preserveShape = TRUE, BPPARAM = BPPARAM, grid = colAutoGrid(exprMat)))
    }

    rownames(exprMat) <- rowNames
    colnames(exprMat) <- colNames
    if (keepZeroesAsNA) {
      exprMat[which(zeroesCoords == 0, arr.ind = TRUE)] <- NA
    }
    names(dimnames(exprMat)) <- c(featureType, "cells")
    return(exprMat)
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
#' @param verbose Should the function show progress messages? Default: TRUE
#' @return Seurat object.
#' @details AUCell ranking matrix is saved in seu@misc$AUCell[["cells_rankings"]]
#' @examples
#' library(SeuratExtend)
#' options(max.print = 20)
#'
#' pbmc <- BuildAUCRank(pbmc)
#' pbmc@misc$AUCell[["cells_rankings"]][3,4]
#' @rdname BuildAUCRank
#' @export

BuildAUCRank <- function(seu, slot = "counts", assay = "RNA", verbose = TRUE){
  library(SeuratObject)
  cr_empty <- is.null(seu@misc$AUCell[["cells_rankings"]])
  gene_match <- identical(colnames(seu), colnames(seu@misc$AUCell$cells_rankings))
  if(!cr_empty & gene_match) {
    if(verbose) message(Sys.time(), " Using pre-exiting cell ranking matrix in 'misc' slot")
    return(seu)
  } else {
    if(!cr_empty & !gene_match) {
      if(verbose) message(
        Sys.time(),
        " Pre-existing cell ranking matrix has different cell IDs with current seurat object. ",
        "Re-build AUC Rank")
    }
    if(verbose) message(Sys.time(), " Build AUC Rank")
    seu@misc$AUCell<-list()
    seu@misc$AUCell[["cells_rankings"]] <-
      .AUCell_buildRankings(
        exprMat = .get_assay_data_compat(seu, slot = slot, assay = assay))
    return(seu)
  }
}

# modified from AUCell package version 1.19

# Method for class "aucellResults"
getAUC <- function(object) {
  SummarizedExperiment::assays(object)[["AUC"]]
}

getRanking <- function(object) {
  SummarizedExperiment::assays(object)[["ranking"]]
}

.AUCell_calcAUC <- function (
    geneSets, rankings,
    nCores = 1,
    normAUC = TRUE,
    aucMaxRank = ceiling(0.05 * nrow(rankings)),
    verbose = TRUE) {
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
  if (methods::is(rankings, "aucellResults")) {
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
    sysname <- Sys.info()[["sysname"]]
    if(sysname == "Windows") {
      import("doParallel")
      import("doRNG")
      import("foreach")
      cl <- makeCluster(nCores)
      registerDoParallel(cl)
      aucMatrix <- foreach(gSetName = names(geneSets)) %dopar% {
        .AUC.geneSet_norm(
          geneSet = geneSets[[gSetName]],
          rankings = rankings, aucMaxRank = ceiling(0.05 * nrow(rankings)),
          gSetName = gSetName)
      }
      stopCluster(cl)
      aucMatrix <- setNames(aucMatrix, names(geneSets))
      aucMatrix <- rlist::list.rbind(aucMatrix)
    } else {
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
  missingSets <- names(geneSets)[
    which(!names(geneSets) %in% rownames(aucMatrix))]
  if (length(missingSets) > 0)
    warning("The AUC for the following sets was not calculated: ",
            paste(missingSets, collapse = ", "))
  missingGenes <- as.matrix(aucMatrix[, c("missing", "nGenes"), drop = FALSE])
  missingPercent <- as.numeric(
    missingGenes[, "missing"])/as.numeric(missingGenes[,"nGenes"])
  missingPercent <- setNames(missingPercent, rownames(missingGenes))
  if (all(missingPercent >= 0.8))
    stop("Fewer than 20% of the genes in the gene sets are included in the rankings.",
         "Check wether the gene IDs in the 'rankings' and 'geneSets' match.")
  if (any(missingPercent > 0.8)) {
    warning("The following gene sets will be excluded from the analysis",
            "(less than 20% of their genes are available):\n",
            paste(names(missingPercent)[which(missingPercent >= 0.8)],
                  collapse = ", "), sep = "", immediate. = TRUE)
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
  return(aucMatrix)
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

calcAUC_matrix <- function(
    GenesetList,
    rankings,
    nCores = 1,
    aucMaxRank = ceiling(0.05 * nrow(rankings)),
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
  if(verbose) message(Sys.time(), " Calculating ", length(GenesetList), " gene set(s)")
  if(is.null(n.items.part)) {
    AUC_matrix <- .AUCell_calcAUC(
      GenesetList, rankings,
      nCores = nCores, aucMaxRank = aucMaxRank,
      verbose = verbose)
  }else{
    if(verbose) message(Sys.time(), " Split gene set(s) into ", n.items.part, " part(s)")
    splited_terms <- split(GenesetList, ceiling(seq_along(GenesetList) * n.items.part / length(GenesetList)))
    for (i in names(splited_terms)) {
      splited_terms[[i]] <- calcAUC_matrix(
        splited_terms[[i]], rankings = rankings, nCores = nCores, aucMaxRank = aucMaxRank, verbose = FALSE)
    }
    AUC_matrix <- rlist::list.rbind(splited_terms)
    names(dimnames(AUC_matrix)) <- c("gene sets", "cells")
  }
  if(verbose) message(Sys.time(), " Done")
  return(AUC_matrix)
}
