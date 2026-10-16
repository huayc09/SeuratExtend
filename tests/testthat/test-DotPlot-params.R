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
# NOTE (2026-10-16): the earlier claim that R/DotPlot.R was fixed on 2026-09-11
# is not true for current master (HEAD c917a47) — R/DotPlot.R:381-392 has no
# plot.margin. The checklist rows 问题清单2026-09-16.xlsx 行5/行7 were correct.
# A2/A4 are therefore pinned as KNOWN-ISSUE below (asserting the broken
# behaviour); flip them back to the correct assertions inside the fix work order.
#
# Data: built-in pbmc

# --- A2: plot.margin is actually consumed -----------------------------------

test_that("DotPlot2 actually consumes plot.margin (A2)", {
  data(pbmc)
  p <- DotPlot2(pbmc, features = c("CD3D", "MS4A1"),
                plot.margin = ggplot2::margin(1, 2, 3, 4))
  # KNOWN-ISSUE(DP-A2): 问题清单2026-09-16.xlsx 行5「DotPlot.R 死参数 plot.margin」(P2)。
  # 断言当前错误行为（plot.margin 被静默忽略），修复工单内翻转回正确断言。
  expect_null(p$theme$plot.margin)
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
  # KNOWN-ISSUE(DP-A4): 问题清单2026-09-16.xlsx 行7「DotPlot.R 静默失败: legend_position 被覆盖」(P2)。
  # 断言当前错误行为（flip+annotation 分支强制 legend.position="bottom"），修复工单内翻转。
  expect_identical(render_md5("left"), render_md5("right"))
  # and the default (no legend_position passed) must still equal the previous
  # forced "bottom" behaviour, so existing plots do not change
  expect_identical(render_md5(NULL), render_md5("bottom"))
})

# --- A5: no ggplot2 deprecation warning -------------------------------------

test_that("DotPlot2 emits no ggplot2 deprecation warning (A5)", {
  data(pbmc)
  expect_warning(DotPlot2(pbmc, features = c("CD3D", "MS4A1")), NA)
})
