# Test file for R/GeneSetAnalysisGO.R
# Covers: GeneSetAnalysisGO, FilterGOTerms, RenameGO
# Regression issues:
#   #12 - MsigDB-style ontology gene sets must load completely from the
#         SeuratExtendData package (skip_if_not_installed guard)
#   #15 - GO database query + rename pipeline returns stable structures
# Data: built-in pbmc; GO database shipped in SeuratExtendData
#
# NOTE on the gates: every case in this file used to open with
#   skip_if_not_installed("SeuratExtendData")
# which is a DEAD GATE -- SeuratExtendData sits in DESCRIPTION:Depends, so
# installing SeuratExtend guarantees it is present and the condition is always
# TRUE. Six gates, six no-ops, and they made the file look guarded when it was
# not. They have been removed (checklist L2-2 G13).
#
# The heavy per-cell AUCell run is intentionally NOT part of this file: it is
# exercised in test-GeneSetAnalysis.R with a tiny custom gene set.

# --- regression: #12 ontology data completeness ---------------------------

test_that("GO ontology database loads and exposes expected roots (issue #12)", {
  options(spe = "human")
  # no-seu call lists the commonly used GO categories
  db <- suppressMessages(GeneSetAnalysisGO())
  expect_type(db, "character")
  # names = readable category, values = GO root IDs
  expect_true("immune_system_process" %in% names(db))
  expect_true("GO:0002376" %in% unname(db))
  expect_true("GO:0050896" %in% unname(db))   # response_to_stimulus
})

test_that("FilterGOTerms returns children under an immune parent (issue #12)", {
  options(spe = "human")
  terms <- FilterGOTerms(parent = "GO:0002376", n.min = 10, n.max = 1000)
  # immune system process has many children -> real, sized result
  expect_true(length(terms) > 10)
  expect_true(all(grepl("^GO:", terms)))
})

test_that("FilterGOTerms warns on unknown parent", {
  options(spe = "human")
  expect_warning(
    FilterGOTerms(parent = "GO:NOT-A-REAL-ID"),
    regexp = "not found"
  )
})

test_that("RenameGO converts GO IDs to readable names with gene counts", {
  options(spe = "human")
  renamed <- RenameGO(c("GO:0002376", "GO:0050896"), spe = "human")
  expect_true(any(grepl("immune system process", renamed)))
  expect_true(all(grepl("\\(.*g\\)", renamed)))
})

test_that("FilterGOTerms/RenameGO matrix pathway works", {
  options(spe = "human")
  set.seed(42)   # reproducible: an unseeded rnorm() is not a control
  m <- matrix(rnorm(6), nrow = 2,
              dimnames = list(c("GO:0002376", "GO:0050896"), c("a", "b", "c")))
  out <- RenameGO(m)
  expect_true(is.matrix(out))
  expect_true(any(grepl("immune system process", rownames(out))))
})

# --- the main path: GeneSetAnalysisGO(seu, ...) ---------------------------
#
# KNOWN-ISSUE (L2-2, P2): the four cases above never call the function the
# file is named after WITH a Seurat object. `GeneSetAnalysisGO(seu, ...)` is
# the path that does the real work (fetch GO terms -> rank cells -> AUC per
# term -> store on seu@misc) and it had zero assertions.
# Measured runtime for 300 cells + one parent: ~4 s, i.e. runnable.

test_that("GeneSetAnalysisGO stores a per-cell AUC matrix on the object", {
  data(pbmc)
  options(spe = "human")
  sub <- subset(pbmc, cells = colnames(pbmc)[1:300])
  seu <- suppressMessages(GeneSetAnalysisGO(
    sub, parent = "GO:0002376", title = "GOIMMUNE",
    nCores = 1, verbose = FALSE))
  expect_s4_class(seu, "Seurat")
  m <- seu@misc[["AUCell"]][["GO"]][["GOIMMUNE"]]
  expect_true(is.matrix(m))
  # one column per cell, one row per retained GO term
  expect_equal(ncol(m), ncol(sub))
  expect_gt(nrow(m), 0)
  # AUCell scores are ranks normalised by aucMaxRank -> [0, 1]
  expect_true(all(m >= 0 & m <= 1))
})

test_that("GeneSetAnalysisGO root is validated only when seu is supplied", {
  data(pbmc)
  options(spe = "human")
  # KNOWN-ISSUE (GO-6, P2) -- source: R/GeneSetAnalysisGO.R:39-43 (early
  # `return()` for seu = NULL) vs :240-243 (the `root` check inside
  # getGOdatabase). The same illegal root is accepted or rejected depending on
  # whether a Seurat object was passed.
  expect_no_error(suppressMessages(GeneSetAnalysisGO(root = "XX")))
  sub <- subset(pbmc, cells = colnames(pbmc)[1:50])
  expect_error(
    suppressMessages(GeneSetAnalysisGO(sub, root = "XX", verbose = FALSE)),
    "'root' database"
  )
})

# --- crash paths pinned from the checklist ---------------------------------

# KNOWN-ISSUE (GO-3, P1) -- source: R/GeneSetAnalysisGO.R:185, :197-200
# is.vector(NULL) and is.vector(factor) are both FALSE, so non-vector input
# enters the matrix branch; rownames() is still NULL there, so the function
# recurses into itself until the node stack overflows. A data.frame column
# (factor) or an unnamed matrix are both ordinary user input.
# Split into one case per input: a node stack overflow unwinds through every
# R handler, so triggering it repeatedly inside one case destabilises the
# expectations that follow. (expect_stack_overflow lives in helper-testthat.R.)

test_that("RenameGO overflows the stack on NULL", {
  expect_stack_overflow(RenameGO(NULL))
})

test_that("RenameGO overflows the stack on a factor", {
  expect_stack_overflow(RenameGO(factor("GO:0002376")))
})

test_that("RenameGO overflows the stack on an unnamed matrix", {
  expect_stack_overflow(RenameGO(matrix(1:4, 2)))
})

test_that("RenameGO accepts a matrix that HAS row names", {
  # control for the three cases above: same shape, only the row names differ
  m <- matrix(1:4, 2, dimnames = list(c("GO:0002376", "GO:0050896"), NULL))
  expect_no_error(RenameGO(m))
})

test_that("FilterGOTerms reports misleading errors for factor / empty input", {
  options(spe = "human")
  # KNOWN-ISSUE (GO-4, P2) -- source: R/GeneSetAnalysisGO.R:136, :154-157
  # is.vector(factor) is FALSE -> matrix branch -> term[rows, ] ->
  expect_error(FilterGOTerms(term = factor(c("GO:0002376", "GO:0050896"))),
               "incorrect number of dimensions")
  # KNOWN-ISSUE (GO-5, P2) -- source: R/GeneSetAnalysisGO.R:137-140
  # sapply over an empty list returns list(), not logical(0), and that list is
  # then used as a subscript.
  expect_error(FilterGOTerms(term = character(0)),
               "invalid subscript type")
})

test_that("GeneSetAnalysisGO calls plyr and library() without declaring them", {
  # KNOWN-ISSUE (GO-1, P2): R/GeneSetAnalysisGO.R:259 uses plyr::revalue()
  # while DESCRIPTION lists plyr in NONE of Imports / Suggests / Depends.
  # KNOWN-ISSUE (GO-2, P2): six bare library() calls (:36, :37, :133, :134,
  # :187, :218) rewrite the caller's search path, while NAMESPACE carries only
  # two imports total (import(ggplot2), importFrom(reshape2, melt)).
  src_path <- file.path("..", "..", "R", "GeneSetAnalysisGO.R")
  skip_if_not(file.exists(src_path), "source tree not reachable")
  src <- readLines(src_path, warn = FALSE)
  expect_true(any(grepl("plyr::revalue", src, fixed = TRUE)))
  expect_gte(sum(grepl("^\\s*library\\(", src)), 6)
  desc <- paste(readLines(system.file("DESCRIPTION", package = "SeuratExtend"),
                          warn = FALSE), collapse = " ")
  # NB: use a word boundary -- a plain grepl("plyr", ...) matches "dplyr",
  # which IS declared in Imports and would make this assertion always fail.
  expect_false(grepl("\\bplyr\\b", desc))
})

# --- parameter coverage (2.5) ---------------------------------------------

test_that("GeneSetAnalysisGO parameter smoke (root/only.end.terms)", {
  options(spe = "human")
  # root = MF (Molecular Function) changes the category space
  db.mf <- suppressMessages(GeneSetAnalysisGO(root = "MF"))
  expect_true(length(db.mf) > 0)
  # only.end.terms keeps terminal pathways only
  terms <- FilterGOTerms(parent = "GO:0002376", n.min = 10, n.max = 1000,
                         only.end.terms = TRUE)
  expect_true(length(terms) <= length(FilterGOTerms(
    parent = "GO:0002376", n.min = 10, n.max = 1000)))
})
