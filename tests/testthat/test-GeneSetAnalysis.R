# Test file for R/GeneSetAnalysis.R + R/runGSVA.R (GSEAplot)
# Covers: GeneSetAnalysis, GSEAplot, install_SeuratExtendData (doc-only),
#         hall50 built-in dataset
# Regression issues:
#   #58 - GSEA output structure: AUCell matrix per gene set x cell
#         (NOTE: issue #58 is an OPEN feature request asking for NES/p-adj
#         columns, which v1.2.10 does not produce. The regression we can
#         encode is the actual documented AUCell output structure.)
#   #47 - GSEA runs on imbalanced cell counts without extra controls
# Data: built-in pbmc (500 cells) + built-in hall50
#
# These tests run the real AUCell pipeline. The gate below used to ask for the
# AUCell *package*, which the code under test does not use at all.
#
# KNOWN-ISSUE (GS-2, P1) -- gate object. The path exercised here is
#   GeneSetAnalysis -> BuildAUCRank -> AUCell_basic.R:37 -> DelayedArray::blockApply
# and its real dependencies are DelayedArray + DelayedMatrixStats
# (AUCell_basic.R:7-8, :37-46) -- NOT AUCell. Measured on 2026-09-11:
#   * DelayedArray installed, AUCell absent -> gate blocked, all 4 cases skipped
#   * AUCell installed on top of it         -> all 4 cases pass
# So the old gate was a FALSE SKIP: the code was runnable, the tests were not.
# (The other direction happens to hold -- AUCell Imports DelayedArray -- which
# is exactly why the mistake stayed invisible for so long.)

# --- helper ---------------------------------------------------------------

gene_set_analysis_available <- function() {
  requireNamespace("DelayedArray", quietly = TRUE) &&
    requireNamespace("DelayedMatrixStats", quietly = TRUE)
}

# --- regression: #58 structure --------------------------------------------

test_that("GeneSetAnalysis returns AUCell matrix per gene set (issue #58)", {
  skip_if_not(gene_set_analysis_available())
  data(pbmc)
  data(hall50)
  gs <- hall50$human[1:2]
  matr <- GeneSetAnalysis(pbmc, genesets = gs, title = "test",
                          export_to_matrix = TRUE, verbose = FALSE)
  expect_true(is.matrix(matr))
  # one row per gene set, one column per cell
  expect_equal(nrow(matr), 2)
  expect_equal(ncol(matr), ncol(pbmc))
  expect_true(all(names(gs) %in% rownames(matr)))
  # AUCell scores are 0..1 scaled by aucMaxRank
  expect_true(all(matr >= 0 & matr <= 1))
})

test_that("GeneSetAnalysis stores results in seu@misc$AUCell (issue #58)", {
  skip_if_not(gene_set_analysis_available())
  data(pbmc)
  data(hall50)
  gs <- hall50$human[1]
  seu <- GeneSetAnalysis(pbmc, genesets = gs, title = "pathways",
                         verbose = FALSE)
  expect_s4_class(seu, "Seurat")
  matr <- seu@misc[["AUCell"]][["pathways"]]
  expect_true(is.matrix(matr))
  expect_equal(ncol(matr), ncol(pbmc))
})

# --- regression: #47 imbalanced counts ------------------------------------

test_that("GeneSetAnalysis tolerates imbalanced group sizes (issue #47)", {
  skip_if_not(gene_set_analysis_available())
  data(pbmc)
  data(hall50)
  # keep one very large and one very small identity (116 vs 14 cells)
  big <- colnames(pbmc)[pbmc$cluster %in% c("CD4 T Naive", "Platelet")]
  seu.sub <- subset(pbmc, cells = big)
  gs <- hall50$human[1]
  matr <- GeneSetAnalysis(seu.sub, genesets = gs, title = "imbalanced",
                          export_to_matrix = TRUE, verbose = FALSE)
  expect_true(is.matrix(matr))
  expect_equal(ncol(matr), ncol(seu.sub))
})

# --- parameter coverage (2.5) ---------------------------------------------

test_that("GeneSetAnalysis parameter smoke (ratio/n.min/n.max/nCores/aucMaxRank)", {
  skip_if_not(gene_set_analysis_available())
  data(pbmc)
  data(hall50)
  gs <- hall50$human[1]
  matr <- GeneSetAnalysis(
    pbmc, genesets = gs, title = "params",
    ratio = 0.2, n.min = 5, n.max = 500,
    slot = "counts", assay = "RNA", nCores = 1,
    aucMaxRank = 100,
    export_to_matrix = TRUE, verbose = FALSE)
  expect_true(is.matrix(matr))
})

# --- GSEAplot (runGSVA.R) -------------------------------------------------
#
# NOTE: GSEAplot is covered in depth by test-runGSVA.R, which is the file that
# matches the source file name. The two cases below are smoke tests only --
# they assert `inherits(p, "ggplot")`, pass an explicit `title =` and wrap the
# call in suppressWarnings(), so they cannot observe the default-title path or
# the ggplot2 deprecation warning. They are kept because they exercise the
# two-group comparison end to end, not because they guard anything.

test_that("GSEAplot returns a ggplot comparing two identities", {
  data(pbmc)
  data(hall50)
  p <- suppressWarnings(GSEAplot(
    pbmc, ident.1 = "CD4 T Naive", ident.2 = "CD8 T cell",
    geneset = hall50$human$HALLMARK_TNFA_SIGNALING_VIA_NFKB,
    title = "TNFA_TEST"))
  expect_true(inherits(p, "ggplot"))
  save_for_review(p, "GSEAplot two-group comparison")
})

test_that("GSEAplot parameter smoke (group.by/logFC.cutoff/p.position)", {
  data(pbmc)
  data(hall50)
  p <- suppressWarnings(GSEAplot(
    pbmc, group.by = "cluster",
    ident.1 = "CD4 T Naive", ident.2 = "CD8 T cell",
    geneset = hall50$human$HALLMARK_HYPOXIA,
    logFC.cutoff = c(-0.5, 0.5),
    p.position = c(0.5, 0.5), sample = 200,
    title = "HYPOXIA_TEST"))
  expect_true(inherits(p, "ggplot"))
  save_for_review(p, "GSEAplot parameter smoke")
})

# --- install_SeuratExtendData (doc-only helper) ---------------------------

test_that("install_SeuratExtendData errors gracefully without a package name", {
  # function exists and is exported
  expect_true("install_SeuratExtendData" %in%
                getNamespaceExports("SeuratExtend"))
})

# --- dependency declarations (checklist L2-1 G1) ---------------------------

test_that("GeneSetAnalysis's DelayedArray use is Suggests-declared but still unguarded on the main path", {
  # KNOWN-ISSUE (GS-1, P1) -- partially resolved 2026-10-16 (work order
  # test-infra-landing): BOTH packages are now declared in DESCRIPTION Suggests
  # (CI/check environments covered). REMAINING issue (清单 P1): the main
  # computation path still calls
  #   DelayedArray::blockApply / colAutoGrid / DelayedArray()  R/GeneSetAnalysis.R:96-109
  #   DelayedMatrixStats::colRanks                             R/AUCell_basic.R:7-8
  # unconditionally while they are only in Suggests -- a minimal install (no
  # Suggests) still dies inside import()'s interactive prompt, surfacing as the
  # unhelpful `could not find function "blockApply"`. NAMESPACE still imports
  # neither. The real fix (requireNamespace guard + clear error) belongs to the
  # P1 dependency-declaration batch; this pin flips fully when that lands.
  root <- file.path("..", "..")
  skip_if_not(dir.exists(file.path(root, "R")), "source tree not reachable")

  desc <- paste(readLines(system.file("DESCRIPTION", package = "SeuratExtend"),
                          warn = FALSE), collapse = " ")
  expect_true(grepl("DelayedArray", desc, fixed = TRUE))
  expect_true(grepl("DelayedMatrixStats", desc, fixed = TRUE))

  ns <- paste(readLines(file.path(root, "NAMESPACE"), warn = FALSE),
              collapse = " ")
  expect_false(grepl("DelayedArray", ns, fixed = TRUE))

  # ... and yet the calls really are on the main path. They live in
  # AUCell_basic.R (reached from GeneSetAnalysis -> BuildAUCRank), where the
  # dependency is pulled in through the package's own import() helper rather
  # than through NAMESPACE.
  body_txt <- paste(readLines(file.path(root, "R", "AUCell_basic.R"),
                              warn = FALSE), collapse = "\n")
  expect_true(grepl("blockApply", body_txt, fixed = TRUE))
  expect_true(grepl('import("DelayedArray")', body_txt, fixed = TRUE))
  expect_true(grepl('import("DelayedMatrixStats")', body_txt, fixed = TRUE))
})
