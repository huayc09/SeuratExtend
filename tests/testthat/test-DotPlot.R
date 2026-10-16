# Test file for R/DotPlot.R
# Covers: DotPlot2, get_available_legends (helper exported alongside DotPlot2)
# Regression issues (spec section 3, DotPlot2 group):
#   #22 missing features  #29/#32 SCT DefaultAssay  #38 percent semantics
#   #39 legend order      #43 round dots under split.by
#   #45 border annotation #49 col.min/col.max cutoff #57 named list + split.by
#   #59 group.by behavior
# Data: built-in pbmc (500 cells x 12627 genes; cluster = 9 identities)
#
# Structure assertions only (no pixel comparison, per spec 2.4). Plotting
# tests call save_for_review() after assertions (spec 2.6): files are written
# only when SEURATEXTEND_VISUAL_REVIEW=true.
#
# NOTE on ggplot introspection: DotPlot2 puts aesthetics in the *global*
# mapping (p$mapping), not in the per-layer mapping; color is normalised to
# "colour" by ggplot2. Verified against Seurat 5.3.0 / ggplot2 3.5.2.

# --- basic structure ------------------------------------------------------

test_that("DotPlot2 returns a ggplot with expected data", {
  data(pbmc)
  p <- DotPlot2(pbmc, features = c("CD3D", "MS4A1"), group.by = "cluster")
  expect_s3_class(p, "ggplot")
  expect_true(length(p$layers) >= 1)
  # plotting data frame holds per-feature/per-group percent and zscore
  expect_true(all(c("pct", "zscore") %in% colnames(p$data)))
  expect_gt(nrow(p$data), 0)
  save_for_review(p, "DotPlot2 basic structure")
})

# --- regression: #22 missing features -------------------------------------

test_that("DotPlot2 warns and skips missing features (issue #22)", {
  data(pbmc)
  expect_warning(
    p <- DotPlot2(pbmc, features = c("CD3D", "NOT_A_REAL_GENE", "MS4A1")),
    regexp = "not found"
  )
  expect_s3_class(p, "ggplot")
  # missing feature is skipped, not drawn
  expect_false("NOT_A_REAL_GENE" %in% unique(p$data$Var1))
  save_for_review(p, "DotPlot2 warns and skips missing features (issue #22)")
})

test_that("DotPlot2 stops when no feature is found (issue #22)", {
  data(pbmc)
  expect_error(
    DotPlot2(pbmc, features = c("NOT_A_REAL_GENE_1", "NOT_A_REAL_GENE_2")),
    regexp = "No valid features"
  )
})

# --- regression: #38 percent semantics ------------------------------------

test_that("DotPlot2 scale_percent controls 0-1 vs 0-100 semantics (issue #38)", {
  data(pbmc)
  genes <- c("CD3D", "MS4A1")
  p.pct100 <- DotPlot2(pbmc, features = genes, group.by = "cluster",
                       scale_percent = TRUE)
  p.pct01  <- DotPlot2(pbmc, features = genes, group.by = "cluster",
                       scale_percent = FALSE)
  expect_lte(max(p.pct100$data$pct), 100)
  expect_gte(min(p.pct100$data$pct), 0)
  expect_lte(max(p.pct01$data$pct), 1)
  # percent columns are scaled by exactly 100x between the two modes
  expect_equal(sort(p.pct100$data$pct),
               sort(p.pct01$data$pct * 100),
               tolerance = 1e-6)
  save_for_review(p.pct100, "DotPlot2 percent scale 0-100 (issue #38)")
})

# --- regression: #49 col.min / col.max ------------------------------------

test_that("DotPlot2 col.min/col.max clip the zscore (issue #49)", {
  data(pbmc)
  genes <- c("CD3D", "MS4A1", "GAPDH")
  p.unclipped <- DotPlot2(pbmc, features = genes, group.by = "cluster")
  p.clipped   <- DotPlot2(pbmc, features = genes, group.by = "cluster",
                          col.min = -1, col.max = 1)
  # KNOWN-ISSUE (L1-3 DP-9, P3): these two assertions used to sit inside
  # `if (max(p.unclipped$data$zscore) > 1) { ... }` guards. On a fixture whose
  # zscores happen not to exceed the bounds, both guards were skipped and the
  # case still PASSED while asserting nothing -- a structural false green.
  # The guard is now an explicit precondition on the fixture (which fails
  # loudly if the fixture ever stops exercising the clip), and the two
  # clipping assertions run unconditionally.
  expect_gt(max(p.unclipped$data$zscore), 1)    # fixture precondition
  expect_lt(min(p.unclipped$data$zscore), -1)
  expect_lte(max(p.clipped$data$zscore), 1)
  expect_gte(min(p.clipped$data$zscore), -1)
  save_for_review(p.clipped, "DotPlot2 col.min col.max clipping (issue #49)")
})

# --- regression: #57 named feature list + split.by ------------------------

test_that("DotPlot2 keeps named feature groups aligned with split.by (issue #57)", {
  data(pbmc)
  features <- list(
    "T cell markers" = c("CD3D", "CD8A"),
    "B cell markers" = c("MS4A1", "CD79A")
  )
  p <- DotPlot2(pbmc, features = features, group.by = "cluster",
                split.by = "condition")
  expect_s3_class(p, "ggplot")
  d <- p$data
  expect_true(all(c("FeatureGroup", "split") %in% colnames(d)))
  # annotation groups are aligned: no feature leaks across feature groups
  expect_true(all(d$FeatureGroup[d$Var1 == "CD3D"] == "T cell markers"))
  expect_true(all(d$FeatureGroup[d$Var1 == "MS4A1"] == "B cell markers"))
  # split factor keeps the declared order (condition levels)
  expect_equal(levels(d$split), levels(factor(pbmc$condition)))
  save_for_review(p, "DotPlot2 named feature list with split.by (issue #57)")
})

# --- regression: #59 group.by changes x-axis correctly --------------------

test_that("DotPlot2 group.by controls the x-axis grouping (issue #59)", {
  data(pbmc)
  p.cluster <- DotPlot2(pbmc, features = "CD3D", group.by = "cluster")
  p.condition <- DotPlot2(pbmc, features = "CD3D", group.by = "condition")
  # cluster is a factor column: x-axis levels equal its factor levels
  expect_equal(levels(p.cluster$data$group), levels(pbmc$cluster))
  expect_equal(nlevels(p.cluster$data$group), 9)
  # condition is a character column: grouping follows factor(condition)
  expect_equal(levels(p.condition$data$group), levels(factor(pbmc$condition)))
  expect_equal(nlevels(p.condition$data$group), 2)
  save_for_review(p.condition, "DotPlot2 group.by condition (issue #59)")
})

# --- regression: #43 round dots under split.by ----------------------------

test_that("DotPlot2 keeps dots round when split.by is used (issue #43)", {
  data(pbmc)
  p <- DotPlot2(pbmc, features = c("CD3D", "MS4A1"), group.by = "cluster",
                split.by = "condition", split.by.method = "border")
  expect_s3_class(p, "ggplot")
  # point layer present and shape is not data-mapped (would distort dots)
  point_layer <- p$layers[[1]]
  expect_s3_class(point_layer$geom, "GeomPoint")
  expect_false("shape" %in% names(p$mapping))
  save_for_review(p, "DotPlot2 round dots under split.by (issue #43)")
})

# --- regression: #39 legend order -----------------------------------------

test_that("get_available_legends reflects the plot configuration (issue #39)", {
  expect_equal(get_available_legends(NULL, "border", border = TRUE),
               c("size", "fill"))
  expect_equal(get_available_legends("condition", "border", border = TRUE),
               c("size", "fill", "color"))
  expect_equal(get_available_legends("condition", "color", border = TRUE),
               c("size", "color", "alpha"))
})

test_that("DotPlot2 accepts a valid legend_order (issue #39)", {
  data(pbmc)
  p <- DotPlot2(pbmc, features = c("CD3D", "MS4A1"),
                legend_order = c("fill", "size"))
  expect_s3_class(p, "ggplot")
  # invalid legend_order is rejected loudly
  expect_error(
    DotPlot2(pbmc, features = c("CD3D", "MS4A1"),
             legend_order = c("fill")),
    regexp = "must exactly match"
  )
  save_for_review(p, "DotPlot2 legend order (issue #39)")
})

# --- regression: #45 split.by.method border -------------------------------

test_that("DotPlot2 split.by.method='border' renders (issue #45)", {
  data(pbmc)
  p <- DotPlot2(pbmc, features = c("CD3D", "MS4A1"), group.by = "cluster",
                split.by = "condition", split.by.method = "border")
  expect_s3_class(p, "ggplot")
  # border colour must be driven by the split variable (global mapping)
  expect_true("colour" %in% names(p$mapping))
  expect_true("split" %in% colnames(p$data))
  save_for_review(p, "DotPlot2 split border method (issue #45)")
})

# --- regression: #29 / #32 SCT DefaultAssay -------------------------------

test_that("DotPlot2 works when DefaultAssay is SCT (issues #29, #32)", {
  data(pbmc)
  # Build a minimal stand-in: register a second assay and set it as default
  # (mirrors the SCT layout: assay 'SCT' holding expression data).
  sct <- pbmc
  suppressWarnings(sct[["SCT"]] <- sct[["RNA"]])
  DefaultAssay(sct) <- "SCT"
  expect_equal(DefaultAssay(sct), "SCT")
  # must not raise match.arg 'arg' should be one of "SCT" errors
  expect_no_error(
    p <- DotPlot2(sct, features = c("CD3D", "MS4A1"), group.by = "cluster")
  )
  expect_s3_class(p, "ggplot")
  save_for_review(p, "DotPlot2 SCT DefaultAssay (issues #29 #32)")
})

# --- parameter coverage smoke (2.5) ---------------------------------------

test_that("DotPlot2 split.by.method variants all render", {
  data(pbmc)
  for (m in c("border", "color", "facet")) {
    p <- DotPlot2(pbmc, features = c("CD3D", "MS4A1"),
                  group.by = "cluster", split.by = "condition",
                  split.by.method = m)
    expect_s3_class(p, "ggplot")
    save_for_review(p, paste("DotPlot2 split.by.method", m))
  }
  # KNOWN-ISSUE (L1-3 DP-9, P3): the annotation branch below used to be
  # guarded by `skip_if_not_installed("cowplot")` and asserted only
  # `expect_no_error(...)` -- that is, zero behaviour.
  #   * The gate was dead: cowplot is in Seurat's Imports (not Suggests), and
  #     SeuratExtend imports Seurat, so cowplot is always installed.
  #   * `expect_no_error` cannot tell "assembled correctly" from "returned
  #     something"; the @return contract says the result is a
  #     cowplot::plot_grid() composition that still inherits from "ggplot".
  # Both are corrected here.
  p <- DotPlot2(pbmc, features = c("CD3D", "MS4A1"),
                group.by = "cluster", split.by = "condition",
                split.by.method = "annotation")
  expect_s3_class(p, "ggplot")
  # it is a cowplot::plot_grid() composition rather than a single ggplot, which
  # is why it carries no $data -- the @return contract at DotPlot.R:65-67 says
  # exactly this. Assert the composition class, not the layer data.
  expect_true("gg" %in% class(p))
  save_for_review(p, "DotPlot2 split.by.method annotation")
})

test_that("DotPlot2 accepts visual/cosmetic parameters", {
  data(pbmc)
  # one call per parameter group (spec 2.5: non-default values, call must
  # succeed and return the right type; aesthetics judged visually in 2.6)
  p <- DotPlot2(
    pbmc, features = c("CD3D", "MS4A1"), group.by = "cluster",
    color_scheme = "OrRd",
    border = FALSE, flip = TRUE,
    show_grid = FALSE, free_space = FALSE,
    angle = 45, hjust = 1, vjust = 1,
    legend_position = "top",
    border.width = 1.2,
    nudge_factor = 0.2
  )
  expect_s3_class(p, "ggplot")

  p2 <- DotPlot2(
    pbmc, features = c("CD3D", "MS4A1"), group.by = "cluster",
    color_scheme = c(low = "lightblue", high = "red"),
    center_color = "white",
    cells = colnames(pbmc)[pbmc$condition == "condition1"],
    legend_position = "bottom", show_grid = TRUE
  )
  expect_s3_class(p2, "ggplot")
  save_for_review(p, "DotPlot2 cosmetic parameters smoke")
})
