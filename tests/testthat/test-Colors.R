# Test file for R/Colors.R
# Covers: color_pro, color_iwh, ryb2rgb, save_colors, show_col2 helpers,
#         plus the colour-maths helpers (seg, seg.x, sca, ang) and the
#         auto-scale builders (handle_color_scheme, scale_cont_auto,
#         scale_disc_auto and their wrappers).
# Data: built-in pbmc

# --- color_pro / color_iwh -------------------------------------------------

test_that("color_pro returns requested number of hex colors", {
  cols <- color_pro(5, col.space = "default")
  expect_type(cols, "character")
  expect_length(cols, 5)
  expect_true(all(grepl("^#", cols)))
  # non-default col.space accepted
  cols2 <- color_pro(5, col.space = "light")
  expect_length(cols2, 5)
})

test_that("color_iwh returns requested number of colors", {
  cols <- color_iwh(5, 1)
  expect_type(cols, "character")
  expect_length(cols, 5)
})

# --- ryb2rgb ---------------------------------------------------------------

test_that("ryb2rgb converts RYB values to hex", {
  out <- ryb2rgb(c(0.1, 0.9, 0.5))
  expect_type(out, "character")
  expect_true(grepl("^#", out))
})

# --- save_colors (issue #8 adjacency: stored colors drive plots) ----------

test_that("save_colors stores var_colors in misc and validates input", {
  data(pbmc)
  seu <- save_colors(pbmc, list("cluster" = "pro_blue", "CD14" = "D",
                                "CD3D" = c("#EEEEEE", "black")))
  expect_s4_class(seu, "Seurat")
  expect_true(all(c("cluster", "CD14", "CD3D") %in%
                    names(seu@misc[["var_colors"]])))
  # non-list input errors
  expect_error(save_colors(pbmc, "not_a_list"), regexp = "must be a list")
})

test_that("loaded var_colors drive DimPlot2 (save_colors pipeline)", {
  data(pbmc)
  seu <- save_colors(pbmc, list("cluster" = "pro_blue"))
  # load.cols=TRUE reads the stored colors without error
  single <- DimPlot2(seu, features = "cluster", combine = FALSE,
                     load.cols = TRUE)
  expect_true(inherits(single[[1]], "ggplot"))
})

# --- show_col2 helper (documented export, #8) ------------------------------

test_that("show_col2 returns a ggplot grid of palettes", {
  # FIXED (FUN-1, P1, same defect as the test-Functions.R pin): show_col2()
  # called list.rbind()/`%>%`/NoAxes() without the package namespace importing
  # them, so this case used to require a manual library(rlist) shim and was
  # order-dependent on whichever neighbour test happened to attach the
  # namespace. The roxygen block now imports all three, so no attach is needed.
  root <- file.path("..", "..")
  skip_if_not(file.exists(file.path(root, "NAMESPACE")),
              "dev-layout source inspection only")
  ns_lines <- readLines(file.path(root, "NAMESPACE"), warn = FALSE)
  expect_true(any(grepl('importFrom(magrittr,"%>%")', ns_lines, fixed = TRUE)))
  expect_true(any(grepl("importFrom(rlist,list.rbind)", ns_lines, fixed = TRUE)))
  expect_true(any(grepl("importFrom(Seurat,NoAxes)", ns_lines, fixed = TRUE)))
  fn_lines <- readLines(file.path(root, "R", "Functions.R"), warn = FALSE)
  expect_true(any(grepl("pal_mtx <- list.rbind(", fn_lines, fixed = TRUE)))

  p <- show_col2(list(A = color_pro(5), B = color_iwh(5, 1)), ncol = 5)
  expect_true(inherits(p, "ggplot"))
})

# ===========================================================================
# Colour maths helpers. Source: R/Colors.R:149-188.
# These are the primitives ryb2rgb() is built from, so they also pin the
# behaviour of the exported ryb2rgb() above.
# ===========================================================================

test_that("seg() interpolates linearly between two bounds by weight", {
  # seg(range, weight): range is a 2-column matrix, weight biases towards
  # the upper bound. ratio = weight[2] / (weight[1] + weight[2]).
  expect_equal(seg(matrix(c(0, 1), nrow = 1), c(1, 1)), 0.5)
  expect_equal(seg(matrix(c(0, 3), nrow = 1), c(1, 3)), 2.25)
  expect_equal(seg(matrix(c(0, 6), nrow = 1), c(1, 1)), 3)
  # vectorised over rows
  expect_equal(seg(matrix(c(0, 6, 2, 4), nrow = 2, byrow = TRUE), c(1, 1)),
               c(3, 3))
})

test_that("seg.x() maps the RYB wheel coordinate onto an RGB triple", {
  # x runs 0..6 around the wheel; ryb_ref holds the 7 anchor colours.
  expect_equal(seg.x(0), c(1, 0.125, 0.125))
  expect_equal(seg.x(1), c(1, 0.5, 0))
  expect_equal(seg.x(6), c(1, 0.125, 0.125))     # 6 == 0, wheel closes
  # interior points interpolate and carry r/y/b names
  mid <- seg.x(2.5)
  expect_equal(unname(mid), c(0.5, 0.75, 0))
  expect_equal(names(mid), c("r", "y", "b"))
})

test_that("sca() rescales a vector onto a target range", {
  expect_equal(sca(c(1, 2, 3)), c(0, 0.5, 1))
  expect_equal(sca(c(1, 2, 3), range = c(0, 10)), c(0, 5, 10))
  # degenerate input (zero spread) collapses to the middle of the range
  expect_equal(sca(c(5, 5, 5)), c(0.5, 0.5, 0.5))
  expect_equal(sca(c(5, 5, 5), range = c(2, 4)), c(3, 3, 3))
})

test_that("ang() returns the RYB-wheel angle of a colour triple", {
  # the smallest channel picks the sector; the other two set the position
  expect_equal(ang(c(0, 0.5, 1)), 10 / 3)
  expect_equal(ang(c(1, 0, 0.5)), 16 / 3)
  expect_equal(ang(c(1, 0.5, 0)), 2 / 3)
  expect_equal(ang(c(0, 1, 1)), 3)
})

# ===========================================================================
# handle_color_scheme(): classifies a user colour scheme. Source :193-239.
# ===========================================================================

test_that("handle_color_scheme() classifies custom, brewer and single input", {
  expect_equal(handle_color_scheme(c("red", "blue")),
               list(type = "custom", colors = c("red", "blue"),
                    is_diverging = FALSE))

  div <- handle_color_scheme("RdBu")
  expect_equal(div$type, "brewer")
  expect_length(div$colors, 11)          # diverging -> full maxcolors
  expect_true(div$is_diverging)

  seq_ <- handle_color_scheme("Blues")
  expect_equal(seq_$type, "brewer")
  expect_length(seq_$colors, 7)          # sequential -> maxcolors - 2
  expect_false(seq_$is_diverging)

  single <- handle_color_scheme("not_a_palette")
  expect_equal(single$type, "single")
  expect_equal(single$colors, "not_a_palette")
})

test_that("handle_color_scheme() resolves the -rev suffix and its aliases", {
  rev_ <- handle_color_scheme("RdYlBu-rev")
  expect_equal(rev_$type, "brewer")
  expect_true(rev_$is_diverging)
  expect_equal(rev_$colors, rev(handle_color_scheme("RdYlBu")$colors))
  # "GnYlRd" is documented as an alias; it is rewritten to RdYlGn-rev
  expect_equal(handle_color_scheme("GnYlRd")$colors,
               handle_color_scheme("RdYlGn-rev")$colors)
})

test_that("handle_color_scheme() rejects qualitative palettes with a warning", {
  # KNOWN-ISSUE (L1 Colors.R): a qualitative RColorBrewer name silently falls
  # back to the default continuous scheme. A warning is emitted but the caller
  # (scale_cont_auto) discards it when the scale is built inside a plot, and the
  # user-requested palette is dropped. Source: R/Colors.R:231-234.
  expect_warning(d <- handle_color_scheme("Set1"), "Qualitative")
  expect_equal(d$type, "default")
  expect_null(d$colors)
  expect_false(d$is_diverging)
})

# ===========================================================================
# scale_cont_auto() / scale_fill_cont_auto() / scale_color_cont_auto()
# Source: R/Colors.R:241-332.
# ===========================================================================

test_that("scale_cont_auto() builds a continuous scale for each scheme shape", {
  # viridis option letter (A-H)
  expect_s3_class(scale_cont_auto("A", type = "color"), "ScaleContinuous")
  # named low/mid/high vector -> gradient2
  expect_s3_class(
    scale_cont_auto(c(low = "white", mid = "grey", high = "red"),
                    type = "color"),
    "ScaleContinuous")
  # named low/high vector -> gradient
  expect_s3_class(
    scale_cont_auto(c(low = "white", high = "red"), type = "color"),
    "ScaleContinuous")
  # a known brewer name goes through handle_color_scheme()
  expect_s3_class(scale_cont_auto("Blues", type = "color"), "ScaleContinuous")
  # NULL short-circuits before any colour work
  expect_null(scale_cont_auto(NULL))
})

test_that("the fill/color wrappers forward type and the NULL short-circuit", {
  expect_s3_class(scale_fill_cont_auto("Blues"), "ScaleContinuous")
  expect_s3_class(scale_color_cont_auto("Blues"), "ScaleContinuous")
  expect_null(scale_color_cont_auto(NULL))
  expect_null(scale_fill_cont_auto(NULL))
})

test_that("scale_cont_auto() crashes on a single colour name", {
  # KNOWN-ISSUE (L1 Colors.R, not yet in any checklist): the "single" branch
  # calls create_scale("gradient", colours) which does
  # do.call(scale_color_gradient, list(colors = ...)), but scale_*_gradient()
  # takes low/high, not colors -> "unused argument (colors = ...)".
  # Source: R/Colors.R:295-301 (create_scale) + :314 (the "single" arm).
  # Reachability: this is the documented usage shape -- a colour stored per
  # feature, e.g. save_colors(seu, list(CD3D = "red")) -- and DimPlot2() then
  # dies while plotting that gene:
  #   DimPlot2(pbmc, features = "CD3D", cols = list(CD3D = "red"))
  expect_error(scale_cont_auto("red", type = "color"), "unused argument")
  expect_error(scale_color_cont_auto("red"), "unused argument")
  expect_error(scale_fill_cont_auto("#FF0000"), "unused argument")
})

test_that("DimPlot2_SelColCont() inherits the single-colour crash", {
  # KNOWN-ISSUE (L1 Colors.R / L3 DimPlot.R): end-to-end reachability of the
  # crash above. A per-feature single colour reaches scale_color_cont_auto()
  # through DimPlot2_SelColCont(), so plotting a continuous feature with
  # cols = list(<feature> = "<colour name>") hard-errors.
  # Source: R/DimPlot.R:745-756.
  data(pbmc)
  expect_error(DimPlot2_SelColCont(pbmc, "CD3D", "red", FALSE),
               "unused argument")
  expect_error(DimPlot2(pbmc, features = "CD3D", cols = list(CD3D = "red")),
               "unused argument")
  # the same feature with the default palette still works
  expect_no_error(DimPlot2(pbmc, features = "CD3D"))
})

test_that("scale_cont_auto() guards its viridis dependency for every colour path", {
  # FIXED (L1 Colors.R, unnumbered pin): viridis is listed only under Suggests,
  # and scale_cont_auto() builds its scale_func list unconditionally at the very
  # top of the body (:252-266) -- including
  #   viridis = scale_color_viridis
  # and that symbol comes from the viridis package (ggplot2 >= 3.0 only exports
  # scale_color_viridis_c / _d). The call site used to be a bare import() whose
  # FALSE was dropped, so on a clean install the list() was evaluated with the
  # symbol unresolvable and every continuous-colour call -- and therefore
  # DimPlot2(seu, features = <gene>) -- failed with
  #   object 'scale_color_viridis' not found
  # It is now import_required("viridis"): import()'s interactive install prompt
  # is preserved, and a refused/failed/skipped install fails fast with an
  # install hint that names viridis instead of dying later at the list().
  # Source: R/Colors.R:250 (import_required).
  body_txt <- paste(deparse(body(SeuratExtend:::scale_cont_auto)),
                    collapse = "\n")
  expect_true(grepl("scale_color_viridis", body_txt, fixed = TRUE))
  expect_true(grepl('import_required("viridis")', body_txt, fixed = TRUE))
  expect_false(grepl('import("viridis")', body_txt, fixed = TRUE))

  desc <- readLines(system.file("DESCRIPTION", package = "SeuratExtend"),
                    warn = FALSE)
  imports_from <- grep("^Imports:", desc)
  suggests_at <- grep("^Suggests:", desc)
  imports <- paste(desc[imports_from:(suggests_at - 1)], collapse = " ")
  expect_false(grepl("viridis", imports, fixed = TRUE))
  expect_true(grepl("viridis", paste(desc[suggests_at:length(desc)],
                                     collapse = " "), fixed = TRUE))
})

test_that("create_scale() is the shared inner builder of scale_cont_auto()", {
  # create_scale is a local closure inside scale_cont_auto; it is the single
  # place the four scale_* constructors are dispatched from, which is why the
  # "single" arm above is broken in only one spot.
  body_txt <- paste(deparse(body(SeuratExtend:::scale_cont_auto)),
                    collapse = "\n")
  expect_true(grepl("create_scale <- function", body_txt, fixed = TRUE))
  expect_true(grepl("do.call(scale_func[[scale_type]], args)", body_txt,
                    fixed = TRUE))
})

# ===========================================================================
# scale_disc_auto() / scale_color_disc_auto() / scale_fill_disc_auto()
# Source: R/Colors.R:334-383.
# ===========================================================================

test_that("scale_disc_auto() picks manual vs brewer by scheme shape", {
  # a named brewer palette is passed straight through to scale_*_brewer
  sc_brewer <- scale_disc_auto("Set1", 3)
  expect_s3_class(sc_brewer, "ScaleDiscrete")
  expect_equal(sc_brewer$palette(3), RColorBrewer::brewer.pal(3, "Set1"))

  # "default" resolves through color_pro()
  sc_default <- suppressMessages(scale_disc_auto("default", 3))
  expect_s3_class(sc_default, "ScaleDiscrete")
  expect_equal(sc_default$palette(3),
               suppressMessages(color_pro(3, col.space = "default")))

  # an explicit colour vector becomes a manual scale
  expect_s3_class(scale_disc_auto(c("#ffffff", "#000000"), 2),
                  "ScaleDiscrete")

  # an unknown single name still produces a manual scale (one value)
  expect_s3_class(scale_disc_auto("#123456", 3), "ScaleDiscrete")

  # NULL short-circuits
  expect_null(scale_disc_auto(NULL, 3))
})

test_that("scale_disc_auto() has a documented fill/one-colour special case", {
  # n == 1 with the "light" palette and type = "fill" returns NULL on purpose
  # (a single light fill is not distinguishable). Source: R/Colors.R:358.
  expect_null(scale_disc_auto("light", 1, type = "fill"))
  # ... but the same call with 3 colours is a normal scale
  expect_s3_class(suppressMessages(scale_disc_auto("light", 3, type = "fill")),
                  "ScaleDiscrete")
})

test_that("scale_color_disc_auto() / scale_fill_disc_auto() forward type", {
  expect_s3_class(suppressMessages(scale_color_disc_auto("light", 3)),
                  "ScaleDiscrete")
  expect_s3_class(suppressMessages(scale_fill_disc_auto("default", 4)),
                  "ScaleDiscrete")
  expect_null(scale_color_disc_auto(NULL, 3))
  expect_null(scale_fill_disc_auto(NULL, 3))
})

test_that("the discrete builders also require RColorBrewer (Suggests only)", {
  # KNOWN-ISSUE (L1 Colors.R, same family as the viridis one above):
  # scale_disc_auto() opens with library(RColorBrewer) (:338) and then reads
  # rownames(brewer.pal.info) (:365); handle_color_scheme() does the same
  # (:200, :217). RColorBrewer sits in Suggests, so on a clean install the
  # library() call itself errors before any palette logic runs.
  # Verified: brewer.pal.info is NOT reachable through the package namespace
  # (lexical lookup returns FALSE) -- only the runtime library() attaches it.
  body_txt <- paste(deparse(body(SeuratExtend:::scale_disc_auto)),
                    collapse = "\n")
  expect_true(grepl("library(RColorBrewer)", body_txt, fixed = TRUE))
  expect_true(grepl("brewer.pal.info", body_txt, fixed = TRUE))

  desc <- readLines(system.file("DESCRIPTION", package = "SeuratExtend"),
                    warn = FALSE)
  imports_from <- grep("^Imports:", desc)
  suggests_at <- grep("^Suggests:", desc)
  imports <- paste(desc[imports_from:(suggests_at - 1)], collapse = " ")
  expect_false(grepl("RColorBrewer", imports, fixed = TRUE))
})
