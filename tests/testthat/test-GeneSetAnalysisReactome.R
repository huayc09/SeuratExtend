# Test file for R/GeneSetAnalysisReactome.R
# Covers: GeneSetAnalysisReactome, FilterReactomeTerms, RenameReactome
# Regression issues:
#   #12 - Reactome gene sets load completely (skip_if_not_installed guard)
# Data: built-in pbmc; Reactome database shipped in SeuratExtendData
#
# NOTE on the gates: all four cases used to open with
# # which is a DEAD GATE -- SeuratExtendData is in DESCRIPTION:Depends, so it is
# always present and the condition is always TRUE. Four gates, four no-ops.
# They have been removed (checklist L2-3).
#
# NOTE on coverage: none of the original four cases ever called
# GeneSetAnalysisReactome(seu, ...), i.e. the function's main path -- filter
# terms -> rank cells -> AUC -> store on seu@misc -- had zero assertions. That
# path takes ~7 s on 147 pathways x 500 cells, so it is affordable; it is
# covered at the bottom of this file.

test_that("Reactome database lists Immune System root (issue #12)", {
  options(spe = "human")
  db <- suppressMessages(GeneSetAnalysisReactome())
  expect_type(db, "character")
  # names are Reactome IDs, values are readable pathway names
  expect_true("R-HSA-168256" %in% names(db))       # Immune System
  expect_true("Immune System" %in% unname(db))
})

test_that("FilterReactomeTerms returns children under Immune System", {
  options(spe = "human")
  terms <- FilterReactomeTerms(parent = "R-HSA-168256", n.min = 10,
                               n.max = 1000)
  expect_true(length(terms) > 10)
  expect_true(all(grepl("^R-HSA-", terms)))
})

test_that("FilterReactomeTerms warns on unknown parent", {
  options(spe = "human")
  expect_warning(
    FilterReactomeTerms(parent = "R-HSA-NOT-REAL"),
    regexp = "not found"
  )
})

test_that("RenameReactome converts IDs to readable names with gene counts", {
  options(spe = "human")
  renamed <- RenameReactome(c("R-HSA-109582", "R-HSA-112316"), spe = "human")
  expect_true(any(grepl("Hemostasis", renamed)))
  expect_true(all(grepl("\\(.*g\\)", renamed)))
})

# --- the main path: GeneSetAnalysisReactome(seu, ...) ---------------------
#
# KNOWN-ISSUE (L2-3, P2): none of the four cases above pass a Seurat object,
# so the path that does the real work (filter terms -> rank cells -> AUC ->
# store on seu@misc) had zero assertions. Measured runtime: ~7 s for
# 147 pathways x 500 cells, i.e. perfectly affordable.

test_that("GeneSetAnalysisReactome stores a per-cell AUC matrix on the object", {
  data(pbmc)
  options(spe = "human")
  seu <- suppressMessages(GeneSetAnalysisReactome(
    pbmc, parent = "R-HSA-168256", title = "REACTOME_IMMUNE",
    nCores = 1, verbose = FALSE))
  expect_s4_class(seu, "Seurat")
  m <- seu@misc[["AUCell"]][["Reactome"]][["REACTOME_IMMUNE"]]
  expect_true(is.matrix(m))
  expect_equal(ncol(m), ncol(pbmc))
  expect_gt(nrow(m), 0)
  expect_true(all(m >= 0 & m <= 1))
})

# --- crash / silence paths pinned from the checklist ----------------------

test_that("an illegal 'spe' yields three different behaviours", {
  # KNOWN-ISSUE (RCT-1/2/3, P2) -- the SAME bad species name, three outcomes
  # depending on which function receives it:
  #   GeneSetAnalysisReactome(spe = "bogus") -> NULL, silently, no message
  #   FilterReactomeTerms(spe = "bogus")     -> node stack overflow
  #   RenameReactome(id, spe = "bogus")      -> the input, unchanged, silently
  # Root cause: Internal.R check_spe() only tests is.null(spe), so "bogus"
  # passes through and Reactome_Data[["bogus"]] is NULL.
  expect_null(GeneSetAnalysisReactome(spe = "bogus"))
  expect_identical(RenameReactome("R-HSA-109582", spe = "bogus"), "R-HSA-109582")
})

test_that("FilterReactomeTerms overflows the stack on an illegal 'spe'", {
  expect_stack_overflow(FilterReactomeTerms(spe = "bogus"))
})

test_that("RenameReactome overflows the stack on NULL", {
  # KNOWN-ISSUE (RCT-4, P2) -- :121, :133-136: is.vector(NULL) is FALSE, so the
  # else branch runs and recurses on rownames(x), which is still NULL.
  expect_stack_overflow(RenameReactome(NULL))
})

test_that("RenameReactome overflows the stack on a factor", {
  expect_stack_overflow(RenameReactome(factor("R-HSA-109582")))
})

test_that("FilterReactomeTerms warns about an unknown parent, then ignores it", {
  options(spe = "human")
  # KNOWN-ISSUE (RCT-5, P2) -- :92-101. The guard warns but does not return or
  # stop, so execution falls through and the function returns the WHOLE
  # size-filtered library instead of the (empty) children of the unknown
  # parent. A caller who mistypes a parent id gets 1600+ plausible pathways
  # rather than an error. Contrast getReactomedatabase(), which stop()s.
  terms <- NULL
  expect_warning(
    terms <- FilterReactomeTerms(parent = "R-HSA-NOT-REAL",
                                 n.min = 10, n.max = 1000),
    "not found"
  )
  expect_gt(length(terms), 1000)
})

test_that("FilterReactomeTerms accepts ids only, the main entry accepts names", {
  options(spe = "human")
  # KNOWN-ISSUE (RCT-6, P2) -- :93 tests `all(parent %in% Ontology$id)`, i.e.
  # IDs only, while getReactomedatabase() (:182-187) also resolves readable
  # names. The same string therefore works through the main entry and is
  # silently ignored here.
  terms <- NULL
  expect_warning(
    terms <- FilterReactomeTerms(parent = "Immune System",
                                 n.min = 10, n.max = 1000),
    "not found"
  )
  expect_gt(length(terms), 1000)   # whole library, not the children
  # the identical name is accepted by the main entry point
  data(pbmc)
  expect_s4_class(
    suppressMessages(GeneSetAnalysisReactome(
      pbmc, parent = "Immune System", title = "BY_NAME",
      nCores = 1, verbose = FALSE)),
    "Seurat"
  )
})
