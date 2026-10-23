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

test_that("FeaturePlot3.grid orders by expression with a >2-dim reduction", {
  # FIXED (problems ledger, FeaturePlot P1): the ordering used hard-coded
  # columns 3:5, which only hold the expression columns when Embeddings()
  # returns two dimensions. With reduction = "pca" (50 dims) it ordered by
  # PC_3..PC_5 instead of CD3D/MS4A1/LYZ. The expression columns are now
  # located as ncol(emb) + 1:3 in both FeaturePlot3 and FeaturePlot3.grid.
  skip_if_not_installed("ggtext")
  data(pbmc)
  feats <- c("CD3D", "MS4A1", "LYZ")
  p <- suppressWarnings(FeaturePlot3.grid(pbmc, features = feats,
                                          reduction = "pca"))
  expect_true(inherits(p, "ggplot"))

  # Row order in the plotting data must follow the mean scaled expression
  # (same transformation as sca() with the default color.range).
  tp <- as.data.frame(lapply(feats, function(f) Seurat::FetchData(pbmc, f)[, 1]))
  rownames(tp) <- colnames(pbmc)
  scaled <- sapply(tp, function(x) {
    if (max(x) == min(x)) rep(0.1, length(x)) else 0.1 + 0.85 * (x - min(x)) / (max(x) - min(x))
  })
  ordered_expr <- rownames(tp)[order(rowMeans(scaled))]
  plotted <- rownames(p$data)[p$data$index == 1]
  expect_identical(plotted, ordered_expr)
  # Monotone in mean expression.
  expect_false(is.unsorted(rowMeans(scaled)[match(plotted, rownames(tp))]))

  # Positive control: the old hard-coded 3:5 columns (PC_3..PC_5) give a
  # different order, so the assertion above really exercises the fix.
  emb <- Seurat::Embeddings(pbmc, reduction = "pca")
  expect_false(identical(plotted, rownames(tp)[order(rowMeans(emb[, 3:5]))]))
})
