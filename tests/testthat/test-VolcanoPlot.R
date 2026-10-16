# Test file for R/VolcanoPlot.R (+ generics.R S3 dispatch)
# Covers: VolcanoPlot (Seurat + default methods)
# Data: built-in pbmc

# --- Seurat method ---------------------------------------------------------

test_that("VolcanoPlot returns a ggplot for two identities", {
  data(pbmc)
  p <- VolcanoPlot(pbmc, ident.1 = "CD4 T Naive", ident.2 = "CD8 T cell")
  expect_true(inherits(p, "ggplot"))
  save_for_review(p, "VolcanoPlot two identities")
})

# --- default (matrix) method ----------------------------------------------

test_that("VolcanoPlot.default accepts matrix + factor", {
  set.seed(42)
  # 4 features x 60 cells, two groups
  matr <- matrix(abs(rnorm(240)), nrow = 4,
                 dimnames = list(c("g1", "g2", "g3", "g4"), NULL))
  f <- factor(rep(c("a", "b"), each = 30))
  p <- VolcanoPlot.default(matr, f = f)
  expect_true(inherits(p, "ggplot"))
})

test_that("VolcanoPlot x/y/threshold parameters are accepted", {
  data(pbmc)
  p <- VolcanoPlot(pbmc, ident.1 = "CD4 T Naive", ident.2 = "CD8 T cell",
                   x = "logFC", y = "tscore", top.n = 20,
                   x.threshold = 0.5, y.threshold = 2)
  expect_true(inherits(p, "ggplot"))
})

# --- errors ----------------------------------------------------------------

test_that("VolcanoPlot generic dispatches correctly", {
  data(pbmc)
  expect_true(inherits(VolcanoPlot(pbmc), "ggplot"))
})

# --- KNOWN-ISSUE pins (2026-09-16 recheck against 问题清单0915) ------------

test_that("VolcanoPlot renders an empty plot for a non-existent identity", {
  data(pbmc)
  # KNOWN-ISSUE (VP-4, P2) -- the identity is never validated anywhere in the
  # call chain. Measured with BOTH identities wrong: no error, 12627 rows,
  # every `significant` flag NA and `point_color` uniformly grey50 -> a
  # well-formed plot carrying ZERO visible points.
  # Contrast WaterfallPlot, which errors on the same input ("missing value
  # where TRUE/FALSE needed"). A silent empty figure is more dangerous than an
  # error, not evidence of robustness.
  p <- VolcanoPlot(pbmc, ident.1 = "NOPE1", ident.2 = "NOPE2")
  expect_s3_class(p, "ggplot")
  expect_true(all(is.na(p$data$significant)))
  expect_equal(unique(p$data$point_color), "grey50")
})

test_that("VolcanoPlot validates y but not x", {
  data(pbmc)
  # KNOWN-ISSUE (VP-2, P2) -- R/VolcanoPlot.R:117 wraps `y` in match.arg while
  # `x` is never validated anywhere in the file, and the file contains no
  # warning() or stop() of its own. Same function, two parameters, two
  # different contracts.
  expect_error(VolcanoPlot(pbmc, ident.1 = "B cell", y = "bogus"),
               "should be one of")
  expect_error(VolcanoPlot(pbmc, ident.1 = "B cell", x = "bogus"),
               "non-numeric argument")
})

test_that("VolcanoPlot reaches ggrepel through a bare library() call", {
  # KNOWN-ISSUE (VP-1, P1) -- partially resolved 2026-10-16 (work order
  # test-infra-landing): ggrepel is now declared in DESCRIPTION Suggests, so
  # check/CI environments see it. REMAINING issue (清单 P1): R/VolcanoPlot.R:183
  # still attaches ggrepel with a bare library() instead of
  # requireNamespace()/`::` -- a minimal install stays fragile. The full fix
  # belongs to the P1 dependency-declaration batch; this pin flips fully then.
  root <- file.path("..", "..")
  skip_if_not(dir.exists(file.path(root, "R")), "source tree not reachable")
  src <- paste(readLines(file.path(root, "R", "VolcanoPlot.R"), warn = FALSE),
               collapse = "\n")
  expect_true(grepl("library(ggrepel)", src, fixed = TRUE))
  expect_true(grepl("geom_text_repel", src, fixed = TRUE))
  desc <- paste(readLines(system.file("DESCRIPTION", package = "SeuratExtend"),
                          warn = FALSE), collapse = " ")
  expect_true(grepl("ggrepel", desc, fixed = TRUE))
})
