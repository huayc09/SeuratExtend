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
  data(pbmc)
  # KNOWN-ISSUE (R/FeaturePlot.R:305, not yet in any checklist): the grid
  # pipeline ends in `tp.c <- list.rbind(tp.c)` WITHOUT the rlist:: prefix,
  # while NAMESPACE only declares import(ggplot2) + importFrom(reshape2, melt).
  # rlist is therefore reachable only because R/ClusterDistrBar.R:55 / :205 do a
  # bare library(rlist) at run time -- this test used to pass purely because
  # test-ClusterDistrBar.R happens to run first. Satisfy the prerequisite
  # explicitly so the test is order-independent, and pin the missing import
  # statically. (Same family as the show_col2()/Functions.R:213 case pinned in
  # test-Colors.R.)
  root <- file.path("..", "..")
  ns_lines <- readLines(file.path(root, "NAMESPACE"), warn = FALSE)
  expect_false(any(grepl("rlist", ns_lines, fixed = TRUE)))
  fn_lines <- readLines(file.path(root, "R", "FeaturePlot.R"), warn = FALSE)
  expect_true(any(grepl("tp.c <- list.rbind(tp.c)", fn_lines, fixed = TRUE)))

  had_rlist <- "package:rlist" %in% search()
  if (!had_rlist) {
    skip_if_not(requireNamespace("rlist", quietly = TRUE), "rlist absent")
    library(rlist)
    on.exit(detach("package:rlist", unload = FALSE), add = TRUE)
  }

  p <- FeaturePlot3.grid(pbmc, features = c("CD3D", "MS4A1", "LYZ"))
  expect_true(inherits(p, "ggplot"))
  save_for_review(p, "FeaturePlot3.grid multiple features")
})
