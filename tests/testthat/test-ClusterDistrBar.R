# Test file for R/ClusterDistrBar.R
# Covers: ClusterDistrBar, BarOfCluster (deprecated wrapper), ClusterDistrPlot
#         (boxplot-by-condition extension)
# Regression issues:
#   #11 - factor level order MUST be respected: reversing levels reverses the
#         stacking order (core factor-handling case)
# Data: built-in pbmc (500 cells; cluster = 9 identities, orig.ident = 2 samples)

# --- basic structure & matrix mode ----------------------------------------

test_that("ClusterDistrBar returns a ggplot with expected data", {
  data(pbmc)
  p <- ClusterDistrBar(origin = pbmc$orig.ident, cluster = pbmc$cluster)
  expect_s3_class(p, "ggplot")
  expect_true(length(p$layers) >= 1)
  expect_true(all(c("Var1", "Var2", "value") %in% colnames(p$data)))
  save_for_review(p, "ClusterDistrBar basic bar plot")
})

test_that("ClusterDistrBar percent mode sums to 100 per sample", {
  data(pbmc)
  m <- ClusterDistrBar(origin = pbmc$orig.ident, cluster = pbmc$cluster,
                       plot = FALSE)
  expect_true(is.matrix(m))
  # default: percent = TRUE, no rev/normalize -> columns sum to 100
  expect_equal(as.numeric(round(colSums(m))), c(100, 100))
  # rownames = cluster levels, colnames = sample levels
  expect_equal(rownames(m), levels(factor(pbmc$cluster)))
  expect_equal(colnames(m), levels(factor(pbmc$orig.ident)))
})

test_that("ClusterDistrBar count mode returns integer totals", {
  data(pbmc)
  m <- ClusterDistrBar(origin = pbmc$orig.ident, cluster = pbmc$cluster,
                       percent = FALSE, plot = FALSE)
  expect_true(all(m == round(m)))
  # total cell number is preserved (500 cells)
  expect_equal(sum(m), ncol(pbmc))
})

test_that("ClusterDistrBar rev+normalize transposes the matrix", {
  data(pbmc)
  m <- ClusterDistrBar(origin = pbmc$orig.ident, cluster = pbmc$cluster,
                       rev = TRUE, normalize = TRUE, plot = FALSE)
  # after rev: rows are samples (2), columns are clusters (9)
  expect_equal(dim(m), c(2, length(levels(factor(pbmc$cluster)))))
})

# --- regression: #11 factor order -----------------------------------------

test_that("ClusterDistrBar respects reversed factor levels (issue #11)", {
  data(pbmc)
  cl_orig <- factor(pbmc$cluster)
  cl_rev  <- factor(pbmc$cluster, levels = rev(levels(cl_orig)))

  p_orig <- ClusterDistrBar(origin = pbmc$orig.ident, cluster = cl_orig,
                            flip = FALSE)
  p_rev  <- ClusterDistrBar(origin = pbmc$orig.ident, cluster = cl_rev,
                            flip = FALSE)

  # first fill level (Var1) follows the input factor levels
  expect_equal(as.character(p_orig$data$Var1[1]), as.character(cl_orig[1]))
  expect_equal(levels(p_orig$data$Var1), levels(cl_orig))
  expect_equal(levels(p_rev$data$Var1), levels(cl_rev))
  # reversal actually flips the stacking order
  expect_false(identical(as.character(levels(p_orig$data$Var1)),
                         as.character(levels(p_rev$data$Var1))))
  save_for_review(p_rev, "ClusterDistrBar reversed factor levels (issue #11)")
})

# --- parameter coverage (2.5) ---------------------------------------------

test_that("ClusterDistrBar accepts width, stack, cols, border, flip", {
  data(pbmc)
  # stack = FALSE -> dodged bars; flip = FALSE -> vertical bars
  p <- ClusterDistrBar(origin = pbmc$orig.ident, cluster = pbmc$cluster,
                       flip = FALSE, stack = FALSE, width = 0.5,
                       cols = "light", border = "black")
  expect_s3_class(p, "ggplot")
  save_for_review(p, "ClusterDistrBar dodge no-flip parameters")
})

test_that("ClusterDistrBar rejects length-mismatched inputs", {
  data(pbmc)
  expect_error(
    ClusterDistrBar(origin = pbmc$orig.ident[1:10], cluster = pbmc$cluster),
    regexp = "same length"
  )
})

# --- deprecated wrapper ----------------------------------------------------

test_that("BarOfCluster warns deprecated and delegates (issue #11 aliases)", {
  data(pbmc)
  expect_warning(
    p <- BarOfCluster(meta = pbmc@meta.data, cluster = "cluster",
                      origin = "orig.ident"),
    regexp = "deprecated"
  )
  expect_s3_class(p, "ggplot")
  save_for_review(p, "BarOfCluster deprecated wrapper")
})

# --- boxplot-by-condition extension ---------------------------------------

test_that("ClusterDistrPlot with condition returns boxplot ggplot", {
  data(pbmc)
  p <- ClusterDistrPlot(origin = pbmc$sample_id, cluster = pbmc$cluster,
                        condition = pbmc$condition)
  expect_s3_class(p, "ggplot")
  expect_true(length(p$layers) >= 1)
  save_for_review(p, "ClusterDistrPlot boxplot by condition")
})

test_that("ClusterDistrPlot rejects rev with condition", {
  data(pbmc)
  expect_error(
    ClusterDistrPlot(origin = pbmc$sample_id, cluster = pbmc$cluster,
                     condition = pbmc$condition, rev = TRUE),
    regexp = "rev=TRUE"
  )
})

test_that("ClusterDistrPlot supports violin/style/box parameters", {
  data(pbmc)
  p <- ClusterDistrPlot(origin = pbmc$sample_id, cluster = pbmc$cluster,
                        condition = pbmc$condition,
                        violin = TRUE, box = FALSE, style = "outline",
                        pt = FALSE)
  expect_s3_class(p, "ggplot")
  save_for_review(p, "ClusterDistrPlot violin style outline")
})

# --- KNOWN-ISSUE pins (2026-09-16 recheck against 问题清单0915) ------------

test_that("ClusterDistrBar returns a bare vector for a single cluster", {
  data(pbmc)
  # KNOWN-ISSUE (CDB-1, P1) -- R/ClusterDistrBar.R:71 (normalize) / :73
  # (percent), failing at :84. `apply(X, 2, f)` collapses to a vector when the
  # result would have a single row, so a one-cluster input does NOT return the
  # matrix the @return contract at :15 promises. Every case above uses the
  # 9-cluster pbmc data, so the degeneracy was never exercised.
  m <- ClusterDistrBar(origin = pbmc$orig.ident,
                       cluster = rep("only", ncol(pbmc)), plot = FALSE)
  expect_null(dim(m))          # i.e. NOT a 1 x 2 matrix
  expect_true(is.numeric(m))
  expect_length(m, length(levels(factor(pbmc$orig.ident))))
})

test_that("ClusterDistrBar crashes when asked to PLOT a single cluster", {
  data(pbmc)
  # the plot=TRUE half of the same defect: the 1-row vector result reaches the
  # plotting code, which expects a data.frame, and never returns a ggplot.
  err <- tryCatch({
    ClusterDistrBar(origin = pbmc$orig.ident,
                    cluster = rep("only", ncol(pbmc)))
    "no condition raised"
  }, error = function(e) conditionMessage(e))
  expect_match(err, "replacement has|data has|rows")
})
