# Test file for R/Vlnplots.R
# Covers: VlnPlot2 (S3: Seurat + default), StackedViolin (alias of default),
#         StackedViolin_v3 (alias of Seurat method)
# Regression issues:
#   #10 - stats.method="wilcox.test": statistics must actually be computed and
#         appear as a layer (GeomBracket), not merely "no error"
#   #51 - statistical comparison must render: annotation aesthetic complete
#   #21 - sort must reorder groups by median expression
# Data: built-in pbmc (500 cells; cluster = 9 identities)

# --- basic structure ------------------------------------------------------

test_that("VlnPlot2 returns a ggplot with violin + box + point data", {
  data(pbmc)
  p <- VlnPlot2(pbmc, features = c("CD3D", "MS4A1"), group.by = "cluster")
  expect_s3_class(p, "ggplot")
  # long-format data: group (f), feature, expression value
  expect_true(all(c("f", "feature", "value") %in% colnames(p$data)))
  expect_true(all(c("CD3D", "MS4A1") %in% levels(p$data$feature)))
  # default: violin + boxplot + points (at least 3 layers)
  expect_gte(length(p$layers), 3)
  save_for_review(p, "VlnPlot2 basic violin box point")
})

test_that("VlnPlot2 respects group.by levels order", {
  data(pbmc)
  p <- VlnPlot2(pbmc, features = "CD3D", group.by = "cluster")
  expect_equal(levels(p$data$f), levels(pbmc$cluster))
})

# --- regression: #21 sort -------------------------------------------------

test_that("VlnPlot2 sort reorders groups by descending median (issue #21)", {
  data(pbmc)
  p.unsorted <- VlnPlot2(pbmc, features = "CD3D", group.by = "cluster")
  p.sorted   <- VlnPlot2(pbmc, features = "CD3D", group.by = "cluster",
                         sort = TRUE)
  # expected order: levels by descending median expression of CD3D
  sub <- p.unsorted$data[p.unsorted$data$feature == "CD3D", ]
  med <- tapply(sub$value, sub$f, median, na.rm = TRUE)
  expected <- names(sort(med, decreasing = TRUE))
  expect_equal(levels(p.sorted$data$f), expected)
  # sorting must differ from the original order for this data
  expect_false(identical(levels(p.sorted$data$f), levels(pbmc$cluster)))
  save_for_review(p.sorted, "VlnPlot2 sort by median expression (issue #21)")
})

test_that("VlnPlot2 sort accepts 'mean' and rejects invalid sort values", {
  data(pbmc)
  p <- VlnPlot2(pbmc, features = "CD3D", group.by = "cluster", sort = "mean")
  expect_s3_class(p, "ggplot")
  expect_error(
    VlnPlot2(pbmc, features = "CD3D", group.by = "cluster", sort = "bogus"),
    regexp = "sort must be"
  )
})

# --- regression: #10 wilcox stats layer -----------------------------------

test_that("VlnPlot2 stat.method='wilcox.test' adds bracket layer (issue #10)", {
  data(pbmc)
  p <- VlnPlot2(pbmc, features = "CD3D", group.by = "cluster",
                stat.method = "wilcox.test")
  expect_s3_class(p, "ggplot")
  # statistics must be visible as a bracket geom layer
  geoms <- sapply(p$layers, function(l) class(l$geom)[1])
  expect_true("GeomBracket" %in% geoms)
  save_for_review(p, "VlnPlot2 wilcox stats brackets (issue #10)")
})

test_that("VlnPlot2 stat.method with comparisons is accepted (issue #51)", {
  data(pbmc)
  p <- VlnPlot2(pbmc, features = "CD3D", group.by = "cluster",
                stat.method = "wilcox.test",
                comparisons = list(c("B cell", "CD8 T cell")))
  expect_s3_class(p, "ggplot")
  geoms <- sapply(p$layers, function(l) class(l$geom)[1])
  expect_true("GeomBracket" %in% geoms)
  save_for_review(p, "VlnPlot2 wilcox with comparisons (issue #51)")
})

test_that("VlnPlot2 defaults to no stats when stat.method='none'", {
  data(pbmc)
  p <- VlnPlot2(pbmc, features = "CD3D", group.by = "cluster",
                stat.method = "none")
  geoms <- sapply(p$layers, function(l) class(l$geom)[1])
  expect_false("GeomBracket" %in% geoms)
})

# --- split.by / f2 --------------------------------------------------------

test_that("VlnPlot2 split.by creates f2 column", {
  data(pbmc)
  p <- VlnPlot2(pbmc, features = "CD3D", group.by = "cluster",
                split.by = "condition")
  expect_s3_class(p, "ggplot")
  expect_true("f2" %in% colnames(p$data))
  expect_equal(levels(p$data$f2), levels(factor(pbmc$condition)))
  save_for_review(p, "VlnPlot2 split.by f2")
})

# --- parameter coverage (2.5) ---------------------------------------------

test_that("VlnPlot2 visual parameters are accepted", {
  data(pbmc)
  p <- VlnPlot2(
    pbmc, features = "CD3D", group.by = "cluster",
    cols = "light",
    style = "fill",
    violin = TRUE, box = FALSE,
    width = 0.7,
    show.mean = TRUE,
    pt = FALSE,
    hide.outlier = TRUE,
    strip.position = "bottom",
    ncol = 1
  )
  expect_s3_class(p, "ggplot")
  save_for_review(p, "VlnPlot2 fill violin show.mean parameters")
})

test_that("VlnPlot2 outline style + points is accepted", {
  data(pbmc)
  p <- VlnPlot2(pbmc, features = "CD3D", group.by = "cluster",
                style = "outline", pt = TRUE, pt.style = "jitter",
                pt.size = 0.5, pt.alpha = 0.5)
  expect_s3_class(p, "ggplot")
  save_for_review(p, "VlnPlot2 outline style jitter")
})

test_that("VlnPlot2 rejects invalid pt.style", {
  data(pbmc)
  expect_error(VlnPlot2(pbmc, features = "CD3D", group.by = "cluster",
                        pt.style = "bogus"),
               regexp = "pt.style")
})

test_that("VlnPlot2 default method works on a matrix", {
  set.seed(42)
  matr <- matrix(abs(rnorm(500)), nrow = 5, ncol = 100,
                 dimnames = list(paste0("g", 1:5), paste0("c", 1:100)))
  f <- factor(rep(c("a", "b"), each = 50))
  p <- VlnPlot2.default(matr, f = f)
  expect_s3_class(p, "ggplot")
})

# --- aliases --------------------------------------------------------------

test_that("StackedViolin aliases VlnPlot2.default; _v3 aliases Seurat method", {
  set.seed(42)
  matr <- matrix(abs(rnorm(500)), nrow = 5, ncol = 100,
                 dimnames = list(paste0("g", 1:5), paste0("c", 1:100)))
  f <- factor(rep(c("a", "b"), each = 50))
  expect_equal(StackedViolin(matr, f = f), VlnPlot2.default(matr, f = f))

  data(pbmc)
  expect_equal(StackedViolin_v3(pbmc, features = "CD3D", group.by = "cluster"),
               VlnPlot2(pbmc, features = "CD3D", group.by = "cluster"))
})

# --- KNOWN-ISSUE pins (2026-09-16 recheck against 问题清单0915) ------------

test_that("VlnPlot2 draws the mean line when box = FALSE", {
  data(pbmc)
  # FIXED (VLN-5) in work order runs/2026-10-22-deps-guard-silentfix: the
  # mean/median crossbars used to live inside the `box` branches of
  # R/Vlnplots.R (VlnPlot2_Plot), so `show.mean = TRUE` was silently ignored
  # unless box = TRUE. The crossbars are now added whenever show.mean = TRUE,
  # independently of `box`. Measured before the fix:
  #   show.mean = TRUE, box = TRUE  -> GeomViolin, GeomPoint, GeomBoxplot,
  #                                    GeomCrossbar, GeomCrossbar
  #   show.mean = TRUE, box = FALSE -> GeomViolin, GeomPoint   (silently gone)
  geoms <- function(p) sapply(p$layers, function(l) class(l$geom)[1])
  crossbars <- function(p) Filter(function(l) inherits(l$geom, "GeomCrossbar"),
                                  p$layers)
  p_box  <- VlnPlot2(pbmc, features = "CD3D", group.by = "cluster",
                     show.mean = TRUE, box = TRUE)
  p_none <- VlnPlot2(pbmc, features = "CD3D", group.by = "cluster",
                     show.mean = TRUE, box = FALSE)
  expect_true("GeomCrossbar" %in% geoms(p_box))
  expect_true("GeomCrossbar" %in% geoms(p_none))
  expect_length(p_none$layers, 4)   # violin + points + mean + median

  # box = TRUE still carries exactly one Mean and one Median layer, and a
  # single colour scale for them (no duplicated layers / scales).
  cb <- crossbars(p_box)
  expect_length(cb, 2)
  expect_setequal(vapply(cb, function(l) l$stat_params$fun, character(1)),
                  c("mean", "median"))
  colour_scales <- Filter(function(s) "colour" %in% s$aesthetics,
                          p_box$scales$scales)
  expect_length(colour_scales, 1)
})

test_that("VlnPlot2 reports dropped comparisons and an empty result", {
  data(pbmc)
  # FIXED (VLN-3) in work order runs/2026-10-22-deps-guard-silentfix:
  # comparisons referencing group names that do not exist used to be dropped
  # silently (no error, warning or message), leaving the caller with a plot
  # that carries no statistics (R/Vlnplots.R, vlnplot2_Stat). The filter now
  # reports dropped comparisons and states when no valid comparison remains.
  expect_message(
    p <- VlnPlot2(pbmc, features = "CD3D", group.by = "cluster",
                  stat.method = "wilcox.test",
                  comparisons = list(c("B cell", "NO_SUCH_GROUP"))),
    "NO_SUCH_GROUP"
  )
  expect_message(
    VlnPlot2(pbmc, features = "CD3D", group.by = "cluster",
             stat.method = "wilcox.test",
             comparisons = list(c("B cell", "NO_SUCH_GROUP"))),
    "No valid comparisons"
  )
  expect_s3_class(p, "ggplot")
  geoms <- sapply(p$layers, function(l) class(l$geom)[1])
  expect_false("GeomBracket" %in% geoms)   # nothing to draw
})

test_that("VlnPlot2 reports a missing feature differently on its two paths", {
  data(pbmc)
  # KNOWN-ISSUE (VLN-2, P2) -- measured 2026-09-16 with a call-stack capture.
  # VlnPlot2 has two entry points and they do NOT agree on how a feature that
  # does not exist is reported:
  #
  #   Seurat path : VlnPlot2.Seurat -> Seu2Matr -> Seurat::FetchData, so the
  #                 caller gets SEURAT's warning
  #                 "The following requested variables were not found: ...".
  #                 R/Vlnplots.R:321 is on a different branch entirely and
  #                 never runs here.
  #   matrix path : R/Vlnplots.R:321 emits message("<feature> not found") and
  #                 raises no warning.
  #
  # So the checklist wording ("VlnPlot2 uses message()") holds only for the
  # matrix entry point. An earlier version of this case asserted the message on
  # the Seurat path, which is why it failed: the hand-written message is
  # unreachable there and the wording belongs to Seurat.
  expect_warning(
    VlnPlot2(pbmc, features = c("CD3D", "NOT_A_REAL_GENE"), group.by = "cluster"),
    "requested variables were not found"
  )

  set.seed(42)
  matr <- matrix(abs(rnorm(500)), nrow = 5, ncol = 100,
                 dimnames = list(paste0("g", 1:5), paste0("c", 1:100)))
  f <- factor(rep(c("a", "b"), each = 50))
  expect_message(VlnPlot2.default(matr, f = f, features = c("g1", "NOT_A_GENE")),
                 "NOT_A_GENE not found")
})
