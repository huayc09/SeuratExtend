# Test file for R/FeaturePlot.R
# Covers: FeaturePlot3 (3-colour 3D-style co-expression), FeaturePlot3.grid
# Data: built-in pbmc
#
# Spec 2.7 notes FeaturePlot3 depends on rgl (3D). On a headless CI the 3D
# scene itself may not render, but the ggarrange legend pipeline that the
# default single-feature path exercises is device-independent. We keep the
# calls minimal and skip only the interactive rgl scene.

test_that("FeaturePlot3 returns a plot object", {
  data(pbmc)
  # ggpubr::get_legend warns when the legend is assembled from components
  p <- suppressWarnings(FeaturePlot3(pbmc, feature.1 = "CD3D"))
  expect_true(inherits(p, "ggplot"))
  save_for_review(p, "FeaturePlot3 single feature")
})

test_that("FeaturePlot3 accepts custom colors", {
  data(pbmc)
  p <- suppressWarnings(FeaturePlot3(pbmc, feature.1 = "CD3D", color = "rgb"))
  expect_true(inherits(p, "ggplot"))
})

test_that("FeaturePlot3 stops when no feature given", {
  data(pbmc)
  expect_error(FeaturePlot3(pbmc, feature.1 = NA, feature.2 = NA,
                            feature.3 = NA),
               regexp = "No feature to plot")
})

test_that("FeaturePlot3.grid renders a feature grid", {
  skip_if_not_installed("ggtext")
  data(pbmc)
  # FIXED by runs/2026-10-22-deps-guard-silentfix (same family as FUN-1): the
  # grid pipeline's `tp.c <- list.rbind(tp.c)` had no rlist import in NAMESPACE
  # and only worked because R/ClusterDistrBar.R happened to attach rlist first.
  # NAMESPACE now declares importFrom(rlist, list.rbind) (added with the
  # show_col2 roxygen exports), so the call resolves inside the package
  # namespace, order-independently and without a manual library() shim.
  root <- file.path("..", "..")
  skip_if_not(file.exists(file.path(root, "NAMESPACE")),
              "dev-layout source inspection only (not available under R CMD check)")
  ns_lines <- readLines(file.path(root, "NAMESPACE"), warn = FALSE)
  expect_true(any(grepl("importFrom(rlist,list.rbind)", ns_lines, fixed = TRUE)))
  fn_lines <- readLines(file.path(root, "R", "FeaturePlot.R"), warn = FALSE)
  expect_true(any(grepl("tp.c <- list.rbind(tp.c)", fn_lines, fixed = TRUE)))

  p <- FeaturePlot3.grid(pbmc, features = c("CD3D", "MS4A1", "LYZ"))
  expect_true(inherits(p, "ggplot"))
  save_for_review(p, "FeaturePlot3.grid multiple features")
})
