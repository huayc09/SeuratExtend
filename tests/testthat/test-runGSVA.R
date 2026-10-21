# Test file for R/runGSVA.R
# Covers: GSEAplot -- the only exported function in this file. Note the file is
#         named runGSVA.R but contains no function of that name (GSV-8).
# Regression issues: L2-4 checklist (10 items: 2x P1 / 1x P2 / 6x P3 / 1x pending)
# Data: built-in pbmc + hall50
#
# GSEAplot is slow (~3-8 s per call: it ranks every gene in the object), so the
# file sticks to the calls it actually needs and pins the rest at source level.
#
# Two warnings come out of every successful call and are swallowed here so they
# neither pollute the log nor make the counts order-dependent:
#   * "Using `size` aesthetic for lines was deprecated in ggplot2 3.4.0" --
#     raised once per session, pinned statically at GSV-4 below;
#   * "Removed N rows containing missing values ... (`geom_tile()`/`geom_bar()`)" --
#     the rank strip and the hit bar both draw outside the scale range.
# This matches how test-GeneSetAnalysis.R already wraps its GSEAplot calls.

gs_human <- function() hall50$human$HALLMARK_INTERFERON_GAMMA_RESPONSE

# --- happy path -----------------------------------------------------------

test_that("GSEAplot() returns the ggarrange object it documents", {
  # Source: R/runGSVA.R:150 (ggarrange) / :152 (return).
  skip_if_not_installed("ggpubr")
  data(pbmc)
  p <- suppressWarnings(
    GSEAplot(pbmc, ident.1 = "CD4 T Naive", title = "MYTITLE",
             geneset = gs_human()))
  expect_s3_class(p, "ggarrange")
  # it keeps ggarrange's class vector, so it still looks like a ggplot
  expect_s3_class(p, "ggplot")
})

test_that("GSEAplot() rejects a geneset with no overlap", {
  # Reference (checklist exclusions e4/e5): the guard at :64-65 is clear and
  # also catches NULL / character(0).
  skip_if_not_installed("ggpubr")
  data(pbmc)
  expect_error(GSEAplot(pbmc, ident.1 = "CD4 T Naive", title = "T",
                        geneset = c("NOPE1", "NOPE2")),
               "No genes in 'geneset' found in expression matrix")
  expect_error(GSEAplot(pbmc, ident.1 = "CD4 T Naive", title = "T",
                        geneset = character(0)),
               "No genes in 'geneset' found in expression matrix")
})

# --- regression: the width crash (GSV-1, fixed) ----------------------------

test_that("GSEAplot() works when 256 or fewer genes survive filtering", {
  # GSV-1 (P1) was FIXED in work order
  # runs/2026-10-21-cellweave-loom-gseaplot (reported by CellWeave): R/runGSVA.R
  # now initializes `width <- 0` before the `if (nrow(scores) > 256)` branch.
  # Previously `width` was assigned only inside that branch and then used
  # unconditionally as `width = width + 1`, so any dataset yielding <= 256
  # ranked genes died ("object 'width' not found", or "non-numeric argument to
  # binary operator" when BiocGenerics' `width` generic was attached). The local
  # initialization shadows any attached generic, so both paths are gone.
  # Source: R/runGSVA.R: width <- 0 immediately before the branch.
  skip_if_not_installed("ggpubr")
  data(pbmc)
  small <- subset(pbmc, features = rownames(pbmc)[1:200])
  p <- expect_no_error(
    suppressWarnings(GSEAplot(small, ident.1 = "CD4 T Naive", title = "T",
                              geneset = rownames(small)[1:50])))
  expect_true(inherits(p, "ggplot"))
})

test_that("GSEAplot() is fine once the gene count clears 256", {
  # Positive control for GSV-1: the very same call succeeds one step above the
  # threshold, so the failure really is the branch, not the subsetting.
  # Source: R/runGSVA.R:116 / :121.
  skip_if_not_installed("ggpubr")
  data(pbmc)
  mid <- subset(pbmc, features = rownames(pbmc)[1:400])
  expect_no_error(suppressWarnings(
    GSEAplot(mid, ident.1 = "CD4 T Naive", title = "T",
             geneset = rownames(mid)[1:50])))
})

# --- KNOWN-ISSUE: parameters that do not do what they say -----------------

test_that("GSEAplot()'s slot argument is dead -- the call hardcodes 'data'", {
  # KNOWN-ISSUE (GSV-3, P2): the formal exists (:34) but :51 writes
  # slot = "data" straight into Seu2Matr, so the user's value is never read.
  # Source: R/runGSVA.R:34 (formal) / :51 (hardcoded).
  skip_if_not_installed("ggpubr")
  data(pbmc)
  p_data <- suppressWarnings(
    GSEAplot(pbmc, ident.1 = "CD4 T Naive", title = "T",
             geneset = gs_human(), slot = "data"))
  p_cnt <- suppressWarnings(
    GSEAplot(pbmc, ident.1 = "CD4 T Naive", title = "T",
             geneset = gs_human(), slot = "counts"))
  expect_identical(class(p_data), class(p_cnt))
  body_txt <- paste(deparse(body(SeuratExtend::GSEAplot)), collapse = "\n")
  expect_true(grepl('slot = "data"', body_txt, fixed = TRUE))
})

test_that("GSEAplot()'s default title is the deparsed geneset call", {
  # KNOWN-ISSUE (GSV-2, P1): title defaults to deparse(substitute(geneset))
  # (:39), so a call written as GSEAplot(pbmc, geneset = hall50$human$X) gets
  # the whole expression pasted into the plot title, prefixed at :67.
  # Source: R/runGSVA.R:39 (default) / :67 (first use, where it is evaluated).
  fm <- formals(SeuratExtend::GSEAplot)
  expect_identical(paste(deparse(fm$title), collapse = ""),
                   "deparse(substitute(geneset))")
  body_txt <- paste(deparse(body(SeuratExtend::GSEAplot)), collapse = "\n")
  expect_true(grepl('paste0("Enrichment Plot: ", title)', body_txt, fixed = TRUE))
})

# --- KNOWN-ISSUE: pinned at the source level ------------------------------

test_that("GSEAplot() attaches three packages with a bare library()", {
  # KNOWN-ISSUE (GSV-6, P3): same bare-library family as the rest of the
  # package. Source: R/runGSVA.R:42-44.
  body_txt <- paste(deparse(body(SeuratExtend::GSEAplot)), collapse = "\n")
  for (pkg in c("ggplot2", "ggpubr", "scales")) {
    expect_true(grepl(sprintf("library(%s)", pkg), body_txt, fixed = TRUE))
  }
})

test_that("GSEAplot() uses grid:: while grid is not declared in DESCRIPTION", {
  # KNOWN-ISSUE (GSV-7, P3): grid::textGrob at :95 (same shape at
  # R/Heatmap.R:291). grid ships with R itself, so this is a packaging-hygiene
  # gap rather than a break.
  # Source: R/runGSVA.R:95.
  body_txt <- paste(deparse(body(SeuratExtend::GSEAplot)), collapse = "\n")
  expect_true(grepl("grid::textGrob", body_txt, fixed = TRUE))
  desc <- read.dcf(system.file("DESCRIPTION", package = "SeuratExtend"))
  expect_false(grepl("grid", paste(unlist(desc), collapse = " ")))
})

test_that("GSEAplot() is the only export in a file named runGSVA.R", {
  # KNOWN-ISSUE (GSV-8, P3): the module name and its only function disagree,
  # which is why this test file had to be named after the file and not the
  # export.
  ns <- asNamespace("SeuratExtend")
  expect_true(exists("GSEAplot", envir = ns, inherits = FALSE))
  expect_false(exists("runGSVA", envir = ns, inherits = FALSE))
  expect_true("GSEAplot" %in% getNamespaceExports("SeuratExtend"))
})

test_that("GSEAplot() draws the connector with a deprecated size aesthetic", {
  # KNOWN-ISSUE (GSV-4, P3): geom_hline(..., size = 1) at :105 -- ggplot2 has
  # preferred linewidth since 3.4.0. Verified to really warn: the first draw in
  # a fresh session emits "Using `size` aesthetic for lines was deprecated in
  # ggplot2 3.4.0". It is a once-per-session lifecycle warning, so whether it
  # surfaces here or in an earlier test depends on run order -- a behavioural
  # assertion would be a false green, hence the source-level pin.
  # Source: R/runGSVA.R:105.
  body_txt <- paste(deparse(body(SeuratExtend::GSEAplot)), collapse = "\n")
  expect_true(grepl("size = 1", body_txt, fixed = TRUE))
  expect_false(grepl("linewidth", body_txt, fixed = TRUE))
})
