# Test file for R/CalcStats.R
# Covers: CalcStats (S3: Seurat + default), CalcScoreGeneral, CalcScoreGeneral_v3,
#         ScoreAndOrder (aliases)
# Regression: #15 (CalcStats -> Heatmap pipeline structure), #58 (GSEA output)
# Data: built-in pbmc (500 cells x 12627 genes, cluster identity 9 levels)
#
# NOTE on CalcStats.default matrix convention (from man page): matr has
# FEATURES in rows and CELLS in columns; f is per-cell grouping (length =
# ncol(matr)). The Seurat method internally passes t = TRUE after fetching
# cells-in-rows data.
#
# NOTE on how this file must be RUN: it calls unexported S3 methods by bare
# name (CalcStats.default). That is correct under the real entry point
# (tests/testthat.R -> test_check(), whose test env is a clone of the package
# namespace) and under pkgload::load_all(), but it FAILS under a plain
# library(SeuratExtend) runner, which only exposes exported names. Recheck
# scripts must therefore use load_all()/test_check(), otherwise these cases
# report phantom "could not find function" errors. See
# deploy/_tests-gap-0916/run-single.R.

# --- structure & values ---------------------------------------------------

test_that("CalcStats returns correct structure and values", {
  data(pbmc)
  res <- CalcStats(pbmc, features = c("CD3D", "MS4A1"), method = "zscore")
  expect_s3_class(res, "data.frame")
  expect_equal(nrow(res), 2)
  # one row per feature; columns are identity levels (cluster, 9 levels)
  expect_equal(ncol(res), length(levels(Idents(pbmc))))
  expect_true(all(c("CD3D", "MS4A1") %in% rownames(res)))
})

test_that("CalcStats zscore row is centered (manual-verifiable small case)", {
  data(pbmc)
  # zscore is computed per feature row across identity columns (9 clusters).
  # Standardization subtracts the row mean, so each row must sum to zero.
  res <- CalcStats(pbmc, features = "CD3D", method = "zscore")
  row.z <- as.numeric(res["CD3D", ])
  expect_equal(sum(row.z), 0, tolerance = 1e-10)
  # with 9 groups the values cannot be all identical
  expect_gt(max(row.z), min(row.z))
})

test_that("CalcStats method='mean' reflects biology of the marker", {
  data(pbmc)
  # CD3D is a T-cell marker: higher mean in CD4 T Naive than in B cell
  res <- CalcStats(pbmc, features = c("CD3D", "MS4A1"), method = "mean")
  expect_gt(res["CD3D", "CD4 T Naive"], res["CD3D", "B cell"])
})

# --- parameter coverage (2.5: every parameter with non-default value) -----

test_that("CalcStats accepts the five non-logfc method values", {
  data(pbmc)
  for (m in c("mean", "median", "zscore", "tscore", "p")) {
    res <- CalcStats(pbmc, features = c("CD3D", "MS4A1"), method = m)
    expect_s3_class(res, "data.frame")
    expect_equal(nrow(res), 2)
  }
})

test_that("CalcStats method='logfc' warns to enable exp.transform", {
  data(pbmc)
  # log2FC on log-normalised data without expm1() back-transform is
  # misleading -> the package warns (documented behaviour)
  expect_warning(
    CalcStats(pbmc, features = "CD3D", method = "logfc"),
    regexp = "exp.transform"
  )
  # with exp.transform = TRUE the warning is gone
  expect_no_warning(
    CalcStats(pbmc, features = "CD3D", method = "logfc",
              exp.transform = TRUE)
  )
})

test_that("CalcStats honors exp.transform, order and n", {
  data(pbmc)
  # exp.transform = TRUE (non-default)
  res <- CalcStats(pbmc, features = "CD3D", method = "mean",
                   exp.transform = TRUE)
  expect_s3_class(res, "data.frame")

  # order + n: top-1 per cluster returns at most n rows per cluster
  res.top <- CalcStats(pbmc, features = c("CD3D", "MS4A1", "GAPDH"),
                       method = "zscore", order = "value", n = 1)
  expect_s3_class(res.top, "data.frame")
  expect_lte(nrow(res.top), length(levels(Idents(pbmc))) * 1)

  # order = "p" exercises the p.threshold branch
  res.p <- CalcStats(pbmc, features = c("CD3D", "MS4A1"),
                     method = "tscore", order = "p", p.threshold = 0.05)
  expect_s3_class(res.p, "data.frame")
})

test_that("CalcStats honors group.by, cells, slot, assay, priority", {
  data(pbmc)
  # group.by: non-default meta.data column -> 2 columns (condition1/2)
  res.g <- CalcStats(pbmc, features = "CD3D", group.by = "condition")
  expect_equal(ncol(res.g), 2)

  # cells: subset + explicit group.by (single group; mean stays finite)
  b.cells <- colnames(pbmc)[pbmc$cluster == "B cell"]
  res.c <- CalcStats(pbmc, features = "CD3D", cells = b.cells,
                     group.by = "cluster", method = "mean")
  expect_s3_class(res.c, "data.frame")
  expect_equal(ncol(res.c), 1)

  # slot + assay: read raw counts from the RNA assay (mean method)
  res.s <- CalcStats(pbmc, features = "CD3D", slot = "counts",
                     assay = "RNA", method = "mean")
  expect_s3_class(res.s, "data.frame")

  # priority: "none" skips the expr-over-meta.data shortcut
  res.p <- CalcStats(pbmc, features = "CD3D", priority = "none")
  expect_s3_class(res.p, "data.frame")
})

# --- error handling -------------------------------------------------------

test_that("CalcStats.default errors when f length mismatches cells", {
  # features in rows (4), cells in columns (3); f length 2 is wrong
  matr <- matrix(rnorm(12), nrow = 4, ncol = 3)
  colnames(matr) <- paste0("cell", 1:3)
  expect_error(CalcStats.default(matr, f = c("a", "b")),
               regexp = "same length")
})

test_that("CalcStats.default accepts matrix input and returns df", {
  set.seed(42)
  matr <- matrix(abs(rnorm(20)), nrow = 4, ncol = 5,
                 dimnames = list(paste0("g", 1:4), paste0("cell", 1:5)))
  f <- factor(c("a", "a", "b", "b", "b"))
  res <- CalcStats.default(matr, f = f, method = "mean")
  expect_s3_class(res, "data.frame")
  expect_equal(nrow(res), 4)                 # one row per feature
  expect_equal(colnames(res), c("a", "b"))   # one column per group
})

test_that("CalcStats rejects invalid method", {
  data(pbmc)
  # KNOWN-ISSUE (CalcStats-1, P2) -- R/CalcStats.R:80-99. CalcStats.default has
  # no match.arg() anywhere in the file (pinned statically in the case below),
  # so an unknown `method` falls through switch() at :90-98 -- which has no
  # unnamed default branch -- and the value stays NULL. The error the user
  # finally sees is raised by a downstream tidyr call and blames an internal
  # helper they never called:
  #   no applicable method for 'drop_na' applied to an object of class "NULL"
  # Neither the argument name nor the offending value appears, so it cannot be
  # acted on. Contrast VlnPlot2(sort = "bogus"), which does validate:
  #   "sort must be FALSE, TRUE, 'median', or 'mean'"
  #
  # This case used to assert only `expect_error(...)`, which passes for any
  # error whatsoever -- including one that has nothing to do with `method`.
  # (A first rewrite of it asserted `expect_false(grepl("method", err))`, which
  # is itself wrong: the message contains the word "method".)
  err <- tryCatch({
    CalcStats(pbmc, features = "CD3D", method = "bogus")
    "no condition raised"
  }, error = function(e) conditionMessage(e))
  expect_match(err, "no applicable method")
  # the message points at an internal helper, not at the argument under test
  expect_match(err, "drop_na")
  # and it never echoes the offending value, so the caller cannot tell which
  # argument is wrong
  expect_false(grepl("bogus", err, fixed = TRUE))
  # no match.arg()-style validation text either
  expect_false(grepl("should be one of", err, fixed = TRUE))
})

test_that("CalcStats.default has no match.arg() to validate method", {
  # Root cause of CalcStats-1, pinned statically so the check is independent of
  # test execution order. Source: R/CalcStats.R.
  src_path <- file.path("..", "..", "R", "CalcStats.R")
  skip_if_not(file.exists(src_path), "source tree not reachable")
  src <- readLines(src_path, warn = FALSE)

  # nothing validates `method` anywhere in the file
  expect_false(any(grepl("match.arg", src, fixed = TRUE)))
  # ... and the value is dispatched by a bare switch() instead
  expect_true(any(grepl("scores <- switch(", src, fixed = TRUE)))
  # the signature still advertises the documented default
  expect_true(any(grepl('method = "zscore"', src, fixed = TRUE)))
})

# --- aliases --------------------------------------------------------------

test_that("CalcScoreGeneral / ScoreAndOrder alias CalcStats.default", {
  set.seed(42)
  matr <- matrix(abs(rnorm(20)), nrow = 4, ncol = 5,
                 dimnames = list(paste0("g", 1:4), paste0("cell", 1:5)))
  f <- factor(c("a", "a", "b", "b", "b"))
  expect_equal(CalcScoreGeneral(matr, f, method = "mean"),
               CalcStats.default(matr, f, method = "mean"))
  expect_equal(ScoreAndOrder(matr, f, method = "mean"),
               CalcStats.default(matr, f, method = "mean"))
})

test_that("CalcScoreGeneral_v3 works on Seurat objects", {
  data(pbmc)
  expect_s3_class(CalcScoreGeneral_v3(pbmc, features = "CD3D"), "data.frame")
})

# --- regression: #15 pipeline structure -----------------------------------

test_that("CalcStats output feeds Heatmap (issue #15)", {
  data(pbmc)
  res <- CalcStats(pbmc, features = c("CD3D", "MS4A1"), method = "zscore")
  # structure contract for the downstream Heatmap:
  #   features as rows, groups as columns, no NA after zscore on 9 groups
  expect_true(all(!is.na(res)))
  expect_false(any(duplicated(rownames(res))))
})
