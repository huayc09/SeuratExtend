# Test file for R/AUCell_basic.R
#
# AUCell_basic.R is a vendored, modified copy of AUCell's ranking/AUC core.
# Only BuildAUCRank() is exported; the other seven functions are internal.
#
# SCOPE: the ranking matrix contract (shape, dimnames, value range), the
# misc-slot caching path, the nCores/aucMaxRank guards in calcAUC_matrix(),
# and the argument validation inside .AUCell_calcAUC().
#
# STYLE: characterisation tests - they pin the behaviour the shipped source
# actually has. Defects found by the audit are tagged `# KNOWN-ISSUE <id>`
# with their source line; when one of those is fixed, the tagged expectation
# is expected to flip, which is the point of pinning it.
#
# This file needs the DelayedArray / DelayedMatrixStats / matrixStats stack.
# (The AUCell *package* is NOT required: R/AUCell_basic.R is a vendored,
# modified copy of the AUCell v1.19 algorithm — no runtime call into AUCell.)
# Without it the ranking fixture cannot be built and everything that
# depends on it is skipped.

auc_deps_available <- function() {
  all(vapply(c("DelayedArray", "DelayedMatrixStats", "matrixStats"),
             requireNamespace, logical(1), quietly = TRUE))
}

# Built once per test run: building the 12627 x 500 ranking matrix is the
# expensive part (~3 s), and every assertion below reads it.
.rank_cache <- new.env(parent = emptyenv())

rank_fixture <- function() {
  if (is.null(.rank_cache$seu)) {
    data(pbmc, envir = environment())
    .rank_cache$seu <- suppressMessages(BuildAUCRank(pbmc, verbose = FALSE))
  }
  .rank_cache$seu
}

test_that("BuildAUCRank is the only export from this file", {
  exports <- getNamespaceExports("SeuratExtend")
  expect_true("BuildAUCRank" %in% exports)
  internals <- c(".AUCell_buildRankings", "getAUC", "getRanking",
                 ".AUCell_calcAUC", ".AUC.geneSet_norm", ".auc", "calcAUC_matrix")
  expect_length(intersect(internals, exports), 0)
})

# --- the ranking contract -------------------------------------------------

test_that("BuildAUCRank() stores a cell ranking matrix in misc$AUCell", {
  skip_if_not(auc_deps_available())
  seu <- rank_fixture()
  data(pbmc, envir = environment())

  expect_s4_class(seu, "Seurat")
  expect_true("AUCell" %in% names(seu@misc))
  rk <- seu@misc$AUCell[["cells_rankings"]]
  expect_true(is.matrix(rk))
  expect_identical(dim(rk), dim(pbmc))
  expect_identical(colnames(rk), colnames(pbmc))   # one column per cell
  expect_identical(rownames(rk), rownames(pbmc))   # one row per feature
  expect_identical(names(dimnames(rk)), c("genes", "cells"))
  # colRanks output: a permutation of 1..nFeatures in every column
  expect_identical(range(rk), c(1L, nrow(pbmc)))
  expect_identical(unname(sort(rk[, 1])), seq_len(nrow(pbmc)))
})

test_that("BuildAUCRank() takes slot and assay from its arguments", {
  skip_if_not(auc_deps_available())
  data(pbmc, envir = environment())
  out <- suppressMessages(BuildAUCRank(pbmc, slot = "data", assay = "RNA",
                                       verbose = FALSE))
  expect_identical(dim(out@misc$AUCell[["cells_rankings"]]), dim(pbmc))
})

# --- the caching path -----------------------------------------------------

test_that("a second call with the same cells reuses the cached ranking", {
  skip_if_not(auc_deps_available())
  seu <- rank_fixture()
  # KNOWN-ISSUE (L2-7 AUC2): the message says "pre-exiting" (sic) and the
  # guard it announces compares colnames only - gene_match is a misnomer.
  # Source: R/AUCell_basic.R:83-87.
  expect_message(
    again <- BuildAUCRank(seu, verbose = TRUE),
    "pre-exiting cell ranking matrix"
  )
  expect_identical(again@misc$AUCell[["cells_rankings"]],
                   seu@misc$AUCell[["cells_rankings"]])
})

test_that("a cell-number mismatch rebuilds the ranking and says so", {
  skip_if_not(auc_deps_available())
  seu <- rank_fixture()
  data(pbmc, envir = environment())
  sub <- pbmc[, 1:100]
  sub@misc$AUCell <- seu@misc$AUCell

  expect_message(
    rebuilt <- BuildAUCRank(sub, verbose = TRUE),
    "different cell IDs"
  )
  rk <- rebuilt@misc$AUCell[["cells_rankings"]]
  expect_identical(dim(rk), c(nrow(pbmc), 100L))
  expect_identical(colnames(rk), colnames(sub))
})

test_that("a rebuild resets the whole misc$AUCell slot, not just its entries", {
  skip_if_not(auc_deps_available())
  seu <- rank_fixture()
  data(pbmc, envir = environment())
  sub <- pbmc[, 1:100]
  sub@misc$AUCell <- seu@misc$AUCell
  sub@misc$AUCell$user_added_entry <- "keep me"

  # KNOWN-ISSUE (L2-7 AUC1): the rebuild path assigns `<- list()` before
  # writing cells_rankings, so any other entry the user put in
  # seu@misc$AUCell is destroyed silently. Source: R/AUCell_basic.R:96.
  rebuilt <- suppressMessages(BuildAUCRank(sub, verbose = TRUE))
  expect_identical(names(rebuilt@misc$AUCell), "cells_rankings")
  expect_false("user_added_entry" %in% names(rebuilt@misc$AUCell))
})

# --- calcAUC_matrix -------------------------------------------------------

test_that("calcAUC_matrix() returns one row per gene set, one column per cell", {
  skip_if_not(auc_deps_available())
  seu <- rank_fixture()
  data(pbmc, envir = environment())
  rk <- seu@misc$AUCell[["cells_rankings"]]
  sets <- list(setA = rownames(pbmc)[1:50], setB = rownames(pbmc)[60:120])

  auc <- suppressMessages(calcAUC_matrix(sets, rk, aucMaxRank = 100, verbose = FALSE))
  expect_true(is.matrix(auc))
  expect_identical(dim(auc), c(2L, ncol(pbmc)))
  expect_identical(rownames(auc), c("setA", "setB"))
  expect_identical(colnames(auc), colnames(pbmc))
  expect_identical(names(dimnames(auc)), c("gene sets", "cells"))
  expect_true(all(auc >= 0 & auc <= 1))
})

test_that("calcAUC_matrix() wraps a bare vector into a single gene set", {
  skip_if_not(auc_deps_available())
  seu <- rank_fixture()
  data(pbmc, envir = environment())
  rk <- seu@misc$AUCell[["cells_rankings"]]

  auc <- suppressMessages(calcAUC_matrix(rownames(pbmc)[1:50], rk,
                                         aucMaxRank = 100, verbose = FALSE))
  expect_identical(rownames(auc), "geneset")
  expect_identical(dim(auc), c(1L, ncol(pbmc)))
})

test_that("calcAUC_matrix() splits gene sets when n.items.part is given", {
  skip_if_not(auc_deps_available())
  seu <- rank_fixture()
  data(pbmc, envir = environment())
  rk <- seu@misc$AUCell[["cells_rankings"]]
  sets <- list(setA = rownames(pbmc)[1:50], setB = rownames(pbmc)[60:120])

  auc <- suppressMessages(calcAUC_matrix(sets, rk, aucMaxRank = 100,
                                         verbose = FALSE, n.items.part = 2))
  expect_identical(dim(auc), c(2L, ncol(pbmc)))
  expect_true(all(c("setA", "setB") %in% rownames(auc)))
})

test_that("calcAUC_matrix() clamps nCores to what the machine can do", {
  skip_if_not(auc_deps_available())
  seu <- rank_fixture()
  data(pbmc, envir = environment())
  rk <- seu@misc$AUCell[["cells_rankings"]]
  sets <- list(setA = rownames(pbmc)[1:50], setB = rownames(pbmc)[60:120])

  # above detectCores() -> falls back to 1
  expect_message(
    invisible(try(quietly(calcAUC_matrix(sets, rk, nCores = 1e6, aucMaxRank = 100,
                                         verbose = FALSE)), silent = TRUE)),
    "nCore is set to 1"
  )
  # above the number of gene sets -> clamped to the number of gene sets.
  # KNOWN-ISSUE (L2-7 AUC3): the clamp happens before the parallel branch, so
  # a request for N cores silently becomes "as many cores as gene sets" and
  # the parallel backend is then loaded. Source: R/AUCell_basic.R:267-273.
  expect_message(
    invisible(try(quietly(calcAUC_matrix(sets, rk, nCores = 3, aucMaxRank = 100,
                                         verbose = FALSE)), silent = TRUE)),
    "nCore is set to 2"
  )
})

test_that("calcAUC_matrix() rejects a non-positive aucMaxRank", {
  skip_if_not(auc_deps_available())
  seu <- rank_fixture()
  data(pbmc, envir = environment())
  rk <- seu@misc$AUCell[["cells_rankings"]]
  sets <- list(setA = rownames(pbmc)[1:50])

  expect_error(
    suppressMessages(calcAUC_matrix(sets, rk, aucMaxRank = 0, verbose = FALSE)),
    "positive value"
  )
})

# --- .AUCell_calcAUC argument validation ---------------------------------

test_that(".AUCell_calcAUC() requires a named list of gene sets", {
  skip_if_not(auc_deps_available())
  seu <- rank_fixture()
  data(pbmc, envir = environment())
  rk <- seu@misc$AUCell[["cells_rankings"]]
  calc <- SeuratExtend:::.AUCell_calcAUC

  expect_error(suppressMessages(calc(c("A", "B"), rk, aucMaxRank = 100, verbose = FALSE)),
               "named list")
  expect_error(suppressMessages(calc(list(rownames(pbmc)[1:50]), rk,
                                     aucMaxRank = 100, verbose = FALSE)),
               "named list")
})

test_that(".AUCell_calcAUC() drops empty gene sets with a warning", {
  skip_if_not(auc_deps_available())
  seu <- rank_fixture()
  data(pbmc, envir = environment())
  rk <- seu@misc$AUCell[["cells_rankings"]]

  expect_warning(
    suppressMessages(SeuratExtend:::.AUCell_calcAUC(
      list(emptySet = character(0), setB = rownames(pbmc)[1:50]), rk,
      aucMaxRank = 100, verbose = FALSE)),
    "Ignoring the following empty sets"
  )
})

test_that(".AUCell_calcAUC() refuses gene sets that are almost entirely absent", {
  skip_if_not(auc_deps_available())
  seu <- rank_fixture()
  rk <- seu@misc$AUCell[["cells_rankings"]]

  expect_error(
    suppressMessages(SeuratExtend:::.AUCell_calcAUC(
      list(x = "NOT_A_REAL_GENE_1", y = "NOT_A_REAL_GENE_2"), rk,
      aucMaxRank = 100, verbose = FALSE)),
    "Fewer than 20%"
  )
})

test_that(".AUCell_calcAUC() warns when aucMaxRank is small and verbose", {
  skip_if_not(auc_deps_available())
  seu <- rank_fixture()
  data(pbmc, envir = environment())
  rk <- seu@misc$AUCell[["cells_rankings"]]

  expect_warning(
    suppressMessages(SeuratExtend:::.AUCell_calcAUC(
      list(setA = rownames(pbmc)[1:50]), rk, aucMaxRank = 100, verbose = TRUE)),
    "aucMaxRank"
  )
})

test_that("normAUC = FALSE cannot run at all", {
  skip_if_not(auc_deps_available())
  seu <- rank_fixture()
  data(pbmc, envir = environment())
  rk <- seu@misc$AUCell[["cells_rankings"]]

  # KNOWN-ISSUE (L2-7 AUC4, P2): the non-normalised branch calls
  # .AUC.geneSet_old(), which does not exist anywhere in the package, so
  # normAUC = FALSE on a single core dies with "could not find function".
  # On >1 core the branch is unreachable anyway because the parallel body
  # hard-codes .AUC.geneSet_norm. Source: R/AUCell_basic.R:143-146,
  # 162-167. The argument is therefore dead in both directions.
  expect_error(
    suppressMessages(SeuratExtend:::.AUCell_calcAUC(
      list(setA = rownames(pbmc)[1:50]), rk,
      aucMaxRank = 100, normAUC = FALSE, verbose = FALSE)),
    "could not find function"
  )
})
