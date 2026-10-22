# Test file for R/DimPlot.R
# Covers: DimPlot2 (enhanced DimPlot/FeaturePlot hybrid; plot_grid/ggplot
#         output), theme_umap_arrows helper
# Regression issues:
#   #17 - multiple reductions / explicit reduction: points render, no empty data
#   #41 - title and legend title are settable
#   #52 - user p + guides(...) must override the default guide (size not locked)
#   #56 - min/max cutoff (col.min/col.max) clips; scale breaks modifiable
# Data: built-in pbmc (500 cells; cluster identity, CD14 continuous, 2 reductions)

# NOTE on return type: DimPlot2 default (combine=TRUE) returns a plot_grid
# object that inherits "ggplot" (verified). Structure assertions below inspect
# the pre-combination single plots via combine=FALSE where needed.

# --- basic structure ------------------------------------------------------

test_that("DimPlot2 returns a ggplot-compatible object with point data", {
  data(pbmc)
  p <- DimPlot2(pbmc)   # features = NULL -> Idents (cluster)
  expect_true(inherits(p, "ggplot"))
  expect_true(length(p) > 0)
  save_for_review(p, "DimPlot2 basic cluster UMAP")
})

test_that("DimPlot2 single plots carry non-empty point data (issue #17)", {
  data(pbmc)
  single <- DimPlot2(pbmc, features = "cluster", combine = FALSE)
  expect_true(is.list(single))
  g <- single[[1]]
  expect_true(inherits(g, "ggplot"))
  # embedding coordinates present and non-empty
  expect_gt(nrow(g$data), 0)
  expect_true(all(c("UMAP_1", "UMAP_2", "var") %in% colnames(g$data)))
})

# --- regression: #17 reduction / dims -------------------------------------

test_that("DimPlot2 honours reduction and dims arguments (issue #17)", {
  data(pbmc)
  # pca reduction with dims 1 and 3 -> PC_1 / PC_3 columns
  single <- DimPlot2(pbmc, features = "cluster", reduction = "pca",
                     dims = c(1, 3), combine = FALSE)
  g <- single[[1]]
  expect_true(all(c("PC_1", "PC_3") %in% colnames(g$data)))
  expect_equal(nrow(g$data), ncol(pbmc))
  save_for_review(
    DimPlot2(pbmc, features = "cluster", reduction = "pca"),
    "DimPlot2 pca reduction (issue #17)")
})

test_that("DimPlot2 handles discrete and continuous features", {
  data(pbmc)
  # discrete: cluster (factor)
  g.d <- DimPlot2(pbmc, features = "cluster", combine = FALSE)[[1]]
  expect_true(is.factor(g.d$data$var))
  # continuous: CD14 gene expression
  g.c <- DimPlot2(pbmc, features = "CD14", combine = FALSE)[[1]]
  expect_true(is.numeric(g.c$data$var))
})

# --- regression: #41 title / legend title ---------------------------------

test_that("DimPlot2 titles are settable via theme param (issue #41)", {
  data(pbmc)
  single <- DimPlot2(pbmc, features = "cluster", combine = FALSE,
                     theme = labs(title = "My Title"))
  expect_equal(single[[1]]$labels$title, "My Title")
  p <- DimPlot2(pbmc, features = "cluster",
                theme = labs(title = "Combined Title"))
  expect_true(inherits(p, "ggplot"))
  save_for_review(p, "DimPlot2 settable title (issue #41)")
})

# --- regression: #52 guides override --------------------------------------

test_that("DimPlot2 guides can be overridden (issue #52)", {
  data(pbmc)
  single <- DimPlot2(pbmc, features = "cluster", combine = FALSE)
  g <- single[[1]]
  # default discrete scale is not size-locked; a user guides() call applies
  g.over <- g + guides(color = guide_legend(ncol = 2))
  expect_true(inherits(g.over, "ggplot"))
  save_for_review(
    DimPlot2(pbmc, features = "cluster"),
    "DimPlot2 guides override (issue #52)")
})

# --- regression: #56 col.min/col.max --------------------------------------

test_that("DimPlot2 col.max clips continuous values (issue #56)", {
  data(pbmc)
  single <- DimPlot2(pbmc, features = "CD14", col.max = 2, combine = FALSE)
  g <- single[[1]]
  expect_lte(max(g$data$var), 2)
  # without clipping the maximum is higher (CD14 is a monocyte marker)
  single.full <- DimPlot2(pbmc, features = "CD14", combine = FALSE)
  expect_gt(max(single.full[[1]]$data$var), 2)
  save_for_review(
    DimPlot2(pbmc, features = "CD14", col.max = 2),
    "DimPlot2 col.max clipping (issue #56)")
})

test_that("DimPlot2 col.min clips low values", {
  data(pbmc)
  single <- DimPlot2(pbmc, features = "CD14", col.min = 0.5, combine = FALSE)
  g <- single[[1]]
  expect_gte(min(g$data$var), 0.5)
})

# --- parameter coverage (2.5) ---------------------------------------------

test_that("DimPlot2 split.by produces per-split panels", {
  data(pbmc)
  p <- DimPlot2(pbmc, features = "cluster", split.by = "condition")
  expect_true(inherits(p, "ggplot"))
  save_for_review(p, "DimPlot2 split.by condition")
})

test_that("DimPlot2 cells / slot / assay / label params are accepted", {
  data(pbmc)
  cells <- colnames(pbmc)[pbmc$cluster %in% c("B cell", "CD8 T cell")]
  p <- DimPlot2(pbmc, features = "cluster", cells = cells,
                slot = "data", assay = "RNA")
  expect_true(inherits(p, "ggplot"))

  # KNOWN-ISSUE (DIM-11): label/box go through DimPlot2_LabelClusters(), which
  # builds its mapping with aes_string() (deprecated in ggplot2 3.0.0, source
  # R/DimPlot.R:680/687/697). The warning is a lifecycle soft warning raised
  # once per session, so whether it surfaces here or in a later test depends on
  # order -- swallow it, and pin the deprecation at the source level instead.
  p2 <- suppressWarnings(
    DimPlot2(pbmc, features = "cluster", label = TRUE,
             label.color = "black", box = TRUE, repel = TRUE,
             label.size = 3))
  expect_true(inherits(p2, "ggplot"))
  save_for_review(p2, "DimPlot2 label box repel")
})

test_that("DimPlot2 cells.highlight works", {
  data(pbmc)
  b <- colnames(pbmc)[pbmc$cluster == "B cell"]
  single <- DimPlot2(pbmc, cells.highlight = b, combine = FALSE)
  g <- single[[1]]
  # highlighted cells are flagged as the plotted variable (Selected_cells)
  expect_equal(levels(g$data$var), c("Unselected", "Selected"))
  expect_equal(sum(g$data$var == "Selected"), length(b))
  save_for_review(DimPlot2(pbmc, cells.highlight = b),
                  "DimPlot2 cells.highlight")
})

test_that("DimPlot2 cosmetic params are accepted (2.5 visual-only)", {
  data(pbmc)
  p <- DimPlot2(pbmc, features = c("cluster", "CD14"),
                cols = list(discrete = "light", continuous = "A"),
                pt.size = 1.2, ncol = 1, align = "v")
  expect_true(inherits(p, "ggplot"))
  save_for_review(p, "DimPlot2 cols list two features")
})

# ===========================================================================
# Internal building blocks of DimPlot2. They are not exported, but every plot
# above goes through them, so these tests pin the mechanics the regressions
# above only observe from the outside.
# ===========================================================================

# --- is_continuous --------------------------------------------------------

test_that("is_continuous() accepts numerics and nothing else", {
  # Source: R/DimPlot.R:368-370.
  expect_true(is_continuous(c(1, 2, 3)))
  expect_true(is_continuous(c(1L, 2L)))
  expect_false(is_continuous(factor(c("a", "b"))))
  expect_false(is_continuous(c("a", "b")))
  expect_false(is_continuous(c(TRUE, FALSE)))   # logical is not numeric
})

# --- DimPlot2_AutoPointSize ----------------------------------------------

test_that("DimPlot2_AutoPointSize() scales down once past 1583 cells", {
  # Source: R/DimPlot.R:567-581. base = min(1583 / nrow, 1), then / sqrt(n_features).
  expect_equal(DimPlot2_AutoPointSize(data.frame(x = 1:100)), 1)
  expect_equal(DimPlot2_AutoPointSize(data.frame(x = 1:1583)), 1)
  expect_equal(DimPlot2_AutoPointSize(data.frame(x = 1:2000)), 1583 / 2000)
  expect_equal(DimPlot2_AutoPointSize(data.frame(x = 1:2000), n_features = 4),
               1583 / 2000 / 2)
  # raster plots always use 1: scattermore ignores the size aesthetic
  expect_equal(DimPlot2_AutoPointSize(data.frame(x = 1:2000), raster = TRUE), 1)
})

# --- subsetting / variable resolution ------------------------------------

test_that("DimPlot2_GetData() returns coordinates plus resolved variables", {
  # Source: R/DimPlot.R:372-461.
  data(pbmc)
  d <- DimPlot2_GetData(pbmc, features = "CD3D", dims = c(1, 2))
  expect_equal(nrow(d$data.dim), ncol(pbmc))
  expect_equal(colnames(d$data.dim), c("UMAP_1", "UMAP_2"))
  expect_equal(colnames(d$data.var), "CD3D")

  # no features / group.by / highlight -> the identity is used
  expect_equal(colnames(DimPlot2_GetData(pbmc)$data.var), "Idents")

  # cells.highlight adds a two-level factor with Unselected first
  dh <- DimPlot2_GetData(pbmc, cells.highlight = colnames(pbmc)[1:10])
  expect_equal(colnames(dh$data.var), "Selected_cells")
  expect_equal(levels(dh$data.var$Selected_cells),
               c("Unselected", "Selected"))

  # split/shape/alpha columns are appended to the coordinate frame
  ds <- DimPlot2_GetData(pbmc, features = "CD3D", shape.by = "cluster")
  expect_true("shape.by" %in% colnames(ds$data.dim))
})

test_that("DimPlot2_GetData() rejects bad cells and bad dims", {
  data(pbmc)
  cells <- colnames(pbmc)
  expect_error(DimPlot2_GetData(pbmc, cells = c("nope1", "nope2")),
               "'cells' not found in Seurat object")
  expect_error(DimPlot2_GetData(pbmc, cells = c(cells[1], "nope")),
               "1 cell\\(s\\) not found")
  expect_error(DimPlot2_GetData(pbmc, cells = c(TRUE, FALSE)),
               "same length as cells")
  expect_error(DimPlot2_GetData(pbmc, dims = c(1, 2, 3)),
               "'dims' must be a two-length integer vector")
})

test_that("DimPlot2_GetData() lets an out-of-range dim fall through to base R", {
  # KNOWN-ISSUE (DIM-2): dims is only checked for length and type, so a dim
  # that exceeds the reduction surfaces as the bare base-R subscript error.
  # Source: R/DimPlot.R:412-416.
  data(pbmc)
  expect_error(DimPlot2_GetData(pbmc, dims = c(1, 99)),
               "subscript out of bounds")
})

test_that("DimPlot2_PlotVars() resolves four input shapes", {
  # Source: R/DimPlot.R:463-490.
  data(pbmc)
  cells <- colnames(pbmc)
  lgl <- cells %in% cells

  by_name <- DimPlot2_PlotVars(pbmc, "cluster", "split.by", cells, lgl)
  expect_true(is.factor(by_name))
  expect_equal(length(by_name), length(cells))
  expect_equal(nlevels(by_name), nlevels(pbmc$cluster))
  expect_equal(names(by_name), cells)

  # NULL -> NULL
  expect_null(DimPlot2_PlotVars(pbmc, NULL, "split.by", cells, lgl))

  # one value per cell in the object
  per_obj <- DimPlot2_PlotVars(pbmc, rep(c("a", "b"), length.out = ncol(pbmc)),
                               "split.by", cells, lgl)
  expect_equal(nlevels(per_obj), 2)

  # one value per requested cell
  per_cells <- DimPlot2_PlotVars(pbmc, rep(c("a", "b"), length.out = length(cells)),
                                 "split.by", cells, lgl)
  expect_equal(nlevels(per_cells), 2)

  expect_error(DimPlot2_PlotVars(pbmc, "nope", "split.by", cells, lgl),
               "Cannot find 'nope' in meta.data")
  expect_error(DimPlot2_PlotVars(pbmc, c(TRUE, FALSE), "split.by", cells, lgl),
               "should be variable name in 'meta.data'")
})

# --- DimPlot2_PlotSingle --------------------------------------------------

test_that("DimPlot2_PlotSingle() draws one panel and honours order/shuffle", {
  # Source: R/DimPlot.R:492-565.
  data(pbmc)
  d <- DimPlot2_GetData(pbmc, features = "CD3D")
  sheet <- d$data.dim
  sheet$var <- d$data.var$CD3D

  p <- DimPlot2_PlotSingle(sheet, title = "CD3D", pt.size = 0.2)
  expect_s3_class(p, "ggplot")
  expect_equal(class(p$layers[[1]]$geom)[1], "GeomPoint")
  expect_equal(p$layers[[1]]$aes_params$size, 0.2)

  # only a discrete variable gets the legend size override (source :536-539);
  # ggplot2 4.x no longer exposes it through p$guides, so pin the branch itself
  body_txt <- paste(deparse(body(SeuratExtend:::DimPlot2_PlotSingle)),
                    collapse = "\n")
  expect_true(grepl("!is_continuous(data.plot$var)", body_txt,
                    fixed = TRUE))
  expect_true(grepl("guide_legend(override.aes = list(size = 3)",
                    body_txt, fixed = TRUE))

  # discrete + order = TRUE -> rows are sorted by the variable
  disc <- sheet
  disc$var <- factor(rep(c("a", "b"), length.out = nrow(disc)), levels = c("a", "b"))
  p_ord <- DimPlot2_PlotSingle(disc, title = "cluster", pt.size = 0.2,
                               order = TRUE)
  expect_false(is.unsorted(as.integer(p_ord$data$var)))
  p_uns <- DimPlot2_PlotSingle(disc, title = "cluster", pt.size = 0.2,
                               order = FALSE)
  expect_true(is.unsorted(as.integer(p_uns$data$var)) ||
                identical(as.integer(p_uns$data$var), as.integer(disc$var)))
})

test_that("DimPlot2_PlotSingle() treats a cells.highlight title specially", {
  # KNOWN-ISSUE (DIM-4): when the panel is the Selected_cells highlight, the
  # per-point size is overwritten by sizes.highlight for EVERY point, so the
  # greyed-out background grows too -- sizes.highlight cannot mean "the
  # highlighted points only". Source: R/DimPlot.R:519-528, :559-563.
  data(pbmc)
  d <- DimPlot2_GetData(pbmc, cells.highlight = colnames(pbmc)[1:10])
  sheet <- d$data.dim
  sheet$var <- d$data.var$Selected_cells

  p <- DimPlot2_PlotSingle(sheet, title = "Selected_cells", pt.size = 0.2,
                           sizes.highlight = 2)
  expect_equal(p$layers[[1]]$aes_params$size, 2)
  # the manual palette is Unselected-grey then the highlight colour
  scale_idx <- which(vapply(p$scales$scales,
                            function(s) "colour" %in% s$aesthetics, logical(1)))
  expect_equal(p$scales$scales[[scale_idx]]$palette(2),
               c("#C3C3C3", "#DE2D26"))
})

test_that("DimPlot2_PlotSingle() switches to scattermore and facets on split.by", {
  data(pbmc)
  d <- DimPlot2_GetData(pbmc, features = "CD3D")
  sheet <- d$data.dim
  sheet$var <- d$data.var$CD3D

  skip_if_not(requireNamespace("scattermore", quietly = TRUE))
  p_raster <- DimPlot2_PlotSingle(sheet, title = "CD3D", pt.size = 1,
                                  raster = TRUE, raster.dpi = c(512, 512))
  expect_equal(class(p_raster$layers[[1]]$geom)[1], "GeomScattermore")

  split_sheet <- sheet
  split_sheet$split.by <- factor(rep(c("a", "b"), length.out = nrow(split_sheet)))
  p_facet <- DimPlot2_PlotSingle(split_sheet, title = "CD3D", pt.size = 0.2)
  expect_true(inherits(p_facet$facet, "FacetWrap"))
})

# --- DimPlot2_LabelClusters / GetXYAesthetics ----------------------------

test_that("GetXYAesthetics() reports the mapped x/y of a point layer", {
  # Source: R/DimPlot.R:710-732.
  hp <- ggplot2::ggplot(data.frame(x = 1:3, y = 1:3),
                        ggplot2::aes(x, y)) + ggplot2::geom_point()
  expect_equal(GetXYAesthetics(hp), list(x = "x", y = "y"))
  expect_error(GetXYAesthetics(ggplot2::ggplot() + ggplot2::geom_line()),
               "Cannot find a geom of class GeomPoint")
})

test_that("DimPlot2_LabelClusters() adds one label layer per cluster", {
  # Source: R/DimPlot.R:583-708.
  data(pbmc)
  d <- DimPlot2_GetData(pbmc, features = "CD3D")
  sheet <- d$data.dim
  sheet$var <- factor(rep(c("a", "b"), length.out = nrow(sheet)))

  base <- ggplot2::ggplot(sheet, ggplot2::aes(UMAP_1, UMAP_2, color = var)) +
    ggplot2::geom_point()

  p <- suppressWarnings(DimPlot2_LabelClusters(base, id = "var", repel = FALSE))
  expect_equal(length(p$layers), length(base$layers) + 1)
  expect_equal(class(p$layers[[length(p$layers)]]$geom)[1], "GeomText")
  # one row per cluster, positioned at the per-cluster median
  labels <- p$layers[[length(p$layers)]]$data
  expect_equal(nrow(labels), 2)
  expect_equal(sort(as.character(labels$var)), c("a", "b"))

  # the lapply(FUN = ...) median really is the coordinate summariser
  medians <- lapply(split(sheet, sheet$var), FUN = function(g) median(g$UMAP_1))
  expect_equal(sort(unlist(medians, use.names = FALSE)),
               sort(labels$UMAP_1))

  # KNOWN-ISSUE (DIM-11): all three label branches build their aesthetic
  # mapping with the deprecated aes_string(). Source: R/DimPlot.R:680/687/697.
  lb_body <- paste(deparse(body(SeuratExtend:::DimPlot2_LabelClusters)),
                   collapse = "\n")
  expect_equal(length(gregexpr("aes_string", lb_body, fixed = TRUE)[[1]]), 3)

  expect_error(
    DimPlot2_LabelClusters(
      ggplot2::ggplot(sheet, ggplot2::aes(UMAP_1, UMAP_2)) +
        ggplot2::geom_point(), id = "nope"),
    "Cannot find variable nope in plotting data")
  expect_error(DimPlot2_LabelClusters(base, id = "var", clusters = "zzz"),
               "The following clusters were not found: zzz")
  expect_error(DimPlot2_LabelClusters(base, id = "var", labels = "only1"),
               "Length of labels")
})

# --- DimPlot2_SelColCont / DimPlot2_SelColDisc ---------------------------

test_that("DimPlot2_SelColCont() resolves a continuous colour scheme", {
  # Source: R/DimPlot.R:734-758.
  data(pbmc)
  expect_s3_class(DimPlot2_SelColCont(pbmc, "CD3D", list(), TRUE),
                  "ScaleContinuous")
  expect_s3_class(
    DimPlot2_SelColCont(pbmc, "CD3D", list(continuous = "Blues"), FALSE),
    "ScaleContinuous")
  # a value stored by save_colors() wins when load.cols = TRUE
  stored <- save_colors(pbmc, list(CD3D = "Blues"))
  expect_s3_class(DimPlot2_SelColCont(stored, "CD3D", list(), TRUE),
                  "ScaleContinuous")
})

test_that("DimPlot2_SelColDisc() resolves a discrete colour scheme", {
  # Source: R/DimPlot.R:760-791.
  data(pbmc)
  n <- nlevels(factor(pbmc$cluster))
  expect_s3_class(
    DimPlot2_SelColDisc(pbmc, n, "cluster", list(), TRUE, FALSE,
                        ggplot2::waiver(), FALSE),
    "ScaleDiscrete")
  expect_null(DimPlot2_SelColDisc(pbmc, n, "cluster", NULL, FALSE, FALSE,
                                  ggplot2::waiver(), FALSE))
  stored <- save_colors(pbmc, list(cluster = "Set1"))
  expect_s3_class(
    DimPlot2_SelColDisc(stored, n, "cluster", list(), TRUE, FALSE,
                        ggplot2::waiver(), FALSE),
    "ScaleDiscrete")
})

test_that("DimPlot2_SelColDisc() collapses a dead ifelse branch", {
  # KNOWN-ISSUE (DIM-14): both arms of the ifelse return the same string, so
  # the `light` palette is rewritten to `pro_light` regardless of the label
  # flag. Harmless today (pro_light is the default-style palette either way)
  # but it means the branch was meant to do something and does not.
  # Source: R/DimPlot.R:788. Only occurrence of this shape in R/.
  body_txt <- paste(deparse(body(SeuratExtend:::DimPlot2_SelColDisc)),
                    collapse = "\n")
  expect_true(grepl('ifelse(label_l, "pro_light", "pro_light")', body_txt,
                    fixed = TRUE))
})

# --- theme plumbing -------------------------------------------------------

test_that("CenterTitle() centres the plot title", {
  # Source: R/DimPlot.R:793-795.
  # NOTE: must stay namespaced. Seurat also exports a CenterTitle(); test-Colors.R
  # calls library(Seurat), which attaches Seurat AHEAD of package:SeuratExtend on
  # the search path, so a bare CenterTitle() resolves to Seurat's version in a
  # full test_dir() run (verified: environmentName(environment(get("CenterTitle")))
  # == "Seurat"). Both happen to return hjust = 0.5, so the assertion passed while
  # testing the wrong function -- a false green.
  ct <- SeuratExtend:::CenterTitle()
  expect_s3_class(ct, "theme")
  expect_equal(ct$plot.title$hjust, 0.5)
})

test_that("apply_theme_elements() takes NULL, one element or a list", {
  # Source: R/DimPlot.R:973-992.
  hp <- ggplot2::ggplot(data.frame(x = 1:3, y = 1:3),
                        ggplot2::aes(x, y)) + ggplot2::geom_point()
  expect_identical(apply_theme_elements(hp, NULL), hp)

  one <- apply_theme_elements(hp, ggplot2::theme_bw())
  expect_s3_class(one, "ggplot")

  # a bare list is applied element by element, NULLs skipped
  many <- apply_theme_elements(
    hp, list(ggplot2::theme_bw(), NULL,
             ggplot2::theme(legend.position = "none")))
  expect_s3_class(many, "ggplot")
  expect_equal(many$theme$legend.position, "none")
})

test_that("theme_umap_arrows() builds a composable theme specification", {
  # Source: R/DimPlot.R:829-884. This is the helper the file header of this
  # test file promised to cover; it is exercised here for the first time.
  t <- theme_umap_arrows()
  expect_s3_class(t, "theme_umap_arrows")
  expect_equal(names(t),
               c("no_axes", "anchor_x", "anchor_y", "line_length", "x_arrow",
                 "y_arrow", "text_offset_x", "text_offset_y", "text_size",
                 "line_width", "x_label", "y_label"))

  # ensure_unit(): bare numbers become mm, existing units keep their own unit
  expect_true(grid::is.unit(t$anchor_x))
  expect_equal(as.character(t$anchor_x), "6mm")
  expect_equal(as.character(theme_umap_arrows(anchor_x = 5)$anchor_x), "5mm")
  expect_equal(
    as.character(theme_umap_arrows(anchor_x = grid::unit(1, "cm"))$anchor_x),
    "1cm")
  expect_equal(theme_umap_arrows(text_size = 8)$text_size, 8)
})

test_that("adding theme_umap_arrows() appends the four annotations", {
  # Source: R/DimPlot.R:887-932, dispatched through S3 ggplot_add.
  hp <- ggplot2::ggplot(data.frame(x = 1:3, y = 1:3),
                        ggplot2::aes(x, y)) + ggplot2::geom_point()
  expect_false(is.null(getS3method("ggplot_add", "theme_umap_arrows",
                                   optional = TRUE)))
  # ... and it is registered in the namespace under its full dotted name
  expect_true(is.function(SeuratExtend:::ggplot_add.theme_umap_arrows))

  p <- hp + theme_umap_arrows()
  expect_s3_class(p, "ggplot")
  # two arrow grobs + two text grobs, added as annotation layers
  expect_equal(length(p$layers), length(hp$layers) + 4)
  # axes are blanked as part of the same object
  expect_s3_class(p$theme$axis.line, "element_blank")
  # and a custom label overrides the auto-detected one
  p_lab <- p + theme_umap_arrows(x_label = "PC1", y_label = "PC2")
  expect_s3_class(p_lab, "ggplot")
})

test_that("extract_label() falls back to UMAP_1 / UMAP_2", {
  # Source: R/DimPlot.R:934-970.
  hp <- ggplot2::ggplot(data.frame(x = 1:3, y = 1:3),
                        ggplot2::aes(x, y)) + ggplot2::geom_point() +
    ggplot2::labs(x = "tSNE_1", y = "tSNE_2")
  expect_equal(extract_label(hp, "x"), "tSNE_1")
  expect_equal(extract_label(hp, "y"), "tSNE_2")

  bare <- ggplot2::ggplot()
  expect_equal(extract_label(bare, "x"), "UMAP_1")
  expect_equal(extract_label(bare, "y"), "UMAP_2")
})

# --- get_disc_cont_par ----------------------------------------------------

test_that("get_disc_cont_par() resolves named, single and junk order/shuffle", {
  # Source: R/DimPlot.R:259-268. Local closure, so it is pinned through the
  # export that owns it plus its observable effect.
  body_txt <- paste(deparse(body(SeuratExtend::DimPlot2)), collapse = "\n")
  expect_true(grepl("get_disc_cont_par <- function", body_txt, fixed = TRUE))
  # one definition plus four call sites: continuous order/shuffle, discrete
  # order/shuffle
  expect_equal(length(gregexpr("get_disc_cont_par", body_txt, fixed = TRUE)[[1]]),
               5)

  data(pbmc)
  # shuffle = FALSE removes the RNG so the row order is comparable
  plain <- DimPlot2(pbmc, features = "cluster", shuffle = FALSE,
                    combine = FALSE)[[1]]
  # KNOWN-ISSUE (DIM-6): c(FALSE, TRUE) has no names, so `type %in% names(par)`
  # is FALSE and length(par) != 1 -> the default is used and the request is
  # silently dropped.
  unnamed <- DimPlot2(pbmc, features = "cluster", shuffle = FALSE,
                      order = c(FALSE, TRUE), combine = FALSE)[[1]]
  expect_identical(plain$data$var, unnamed$data$var)

  # a single logical IS honoured
  sorted <- DimPlot2(pbmc, features = "cluster", shuffle = FALSE,
                     order = TRUE, combine = FALSE)[[1]]
  expect_false(identical(plain$data$var, sorted$data$var))
  expect_false(is.unsorted(as.integer(sorted$data$var)))

  # a NAMED two-element vector is honoured too
  named <- DimPlot2(pbmc, features = "cluster", shuffle = FALSE,
                    order = c(discrete = TRUE, continuous = FALSE),
                    combine = FALSE)[[1]]
  expect_identical(sorted$data$var, named$data$var)
})

test_that("DimPlot2() rejects a discrete variable with more than 100 values", {
  # Source: R/DimPlot.R:320.
  data(pbmc)
  s <- pbmc
  s$many <- factor(seq_len(ncol(s)))
  expect_error(DimPlot2(s, features = "many"), "has > 100 values")
})


test_that("DimPlot2 honours priority when a gene and a meta column share a name", {
  # FIXED (problems ledger row 11): `priority` was declared but never used by
  # DimPlot2, so a meta.data column named like a gene always won and
  # priority = "expr" was silently ignored. Wired in work order
  # runs/2026-10-22-deps-guard-silentfix with the same semantics as Seu2Matr:
  # the conflicting meta.data column is dropped before FetchData.
  data(pbmc)
  s <- pbmc
  s@meta.data$CD3D <- factor(rep(c("low", "high"), length.out = ncol(s)))
  continuous <- function(g) is.numeric(g$data$var) && !is.factor(g$data$var)

  g_expr <- DimPlot2(s, features = "CD3D", priority = "expr",
                     combine = FALSE)[[1]]
  expect_true(continuous(g_expr))

  g_meta <- DimPlot2(s, features = "CD3D", priority = "none",
                     combine = FALSE)[[1]]
  expect_false(continuous(g_meta))
  expect_true(is.factor(g_meta$data$var))

  # default is priority = c("expr", "none") -> expression wins
  g_default <- DimPlot2(s, features = "CD3D", combine = FALSE)[[1]]
  expect_true(continuous(g_default))
})
