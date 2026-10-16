# Test file for R/DotPlot.R — parameter-level regressions
#
# Deliberately kept separate from test-DotPlot.R: these three assertions guard
# the 2026-09-11 fix set and were reviewed on their own, so they can be taken
# or dropped without touching the main functional test file.
#
# Covered fixes:
#   A2 - plot.margin was declared but never passed to theme() (dead parameter:
#        passing a value changed nothing)
#   A4 - legend_position was unconditionally overridden by
#        theme(legend.position = "bottom") in the flip + annotation branch
#   A5 - element_rect(size = ) is deprecated since ggplot2 3.4.0
#
# NOTE (2026-09-16): the delivery checklist 问题清单0915 still lists all three
# of these as OPEN DotPlot.R problems ("dead parameter plot.margin",
# "legend_position silently overridden", "deprecated element_rect(size=)").
# That is STALE. Those rows were taken from the 2026-09-10 master list, and
# R/DotPlot.R was fixed on 2026-09-11 09:12. Current source:
#   DotPlot.R:392  plot.margin = plot.margin          <- consumed
#   DotPlot.R:528  if (missing(legend_position))      <- user value honoured
#   DotPlot.R:388  element_rect(fill = NA, linewidth = 0)
# Station L1-3 reached the same verdict independently: its checklist marks the
# old DP-1 / DP-2 / DP-4 / DP-6 rows as fixed. The cases below are the guards
# that keep them fixed.
#
# Data: built-in pbmc

# --- A2: plot.margin is actually consumed -----------------------------------

test_that("DotPlot2 actually consumes plot.margin (A2)", {
  data(pbmc)
  p <- DotPlot2(pbmc, features = c("CD3D", "MS4A1"),
                plot.margin = ggplot2::margin(1, 2, 3, 4))
  expect_false(is.null(p$theme$plot.margin),
               info = "plot.margin was declared but never passed to theme()")
  expect_equal(as.numeric(p$theme$plot.margin), c(1, 2, 3, 4))
})

# --- A4: legend_position is honoured with flip + annotation -----------------

test_that("DotPlot2 honours legend_position with flip + annotation (A4)", {
  data(pbmc)
  render_md5 <- function(lp) {
    args <- list(seu = pbmc, features = c("CD3D", "MS4A1"), group.by = "cluster",
                 split.by = "condition", split.by.method = "annotation", flip = TRUE)
    if (!is.null(lp)) args$legend_position <- lp
    p <- do.call(DotPlot2, args)
    f <- tempfile(fileext = ".png")
    ggplot2::ggsave(f, p, width = 7, height = 5, dpi = 72)
    as.character(tools::md5sum(f))
  }
  # before the fix all settings rendered identically (md5 equal)
  expect_false(identical(render_md5("left"), render_md5("right")))
  # and the default (no legend_position passed) must still equal the previous
  # forced "bottom" behaviour, so existing plots do not change
  expect_identical(render_md5(NULL), render_md5("bottom"))
})

# --- A5: no ggplot2 deprecation warning -------------------------------------

test_that("DotPlot2 emits no ggplot2 deprecation warning (A5)", {
  data(pbmc)
  expect_warning(DotPlot2(pbmc, features = c("CD3D", "MS4A1")), NA)
})
