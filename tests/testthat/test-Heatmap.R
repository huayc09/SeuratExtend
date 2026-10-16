# Test file for R/Heatmap.R
# Covers: Heatmap (ggplot tile heatmap; input matrix e.g. from CalcStats)
# Regression issues:
#   #6, #15 - feature_text_subset shows the requested gene names; behaviour
#             with non-existent genes in the subset should be graceful.
# Data: built-in pbmc (500 cells x 12627 genes)

# helper: build a small zscore matrix to feed Heatmap
heatmap_input <- function() {
  data(pbmc)
  CalcStats(pbmc, features = c("CD3D", "MS4A1", "GAPDH"), method = "zscore")
}

# --- basic structure ------------------------------------------------------

test_that("Heatmap returns a ggplot with expected data", {
  score <- heatmap_input()
  p <- Heatmap(score)
  expect_s3_class(p, "ggplot")
  expect_true(length(p$layers) >= 1)
  # melt()ed data carries feature id, column variable and value
  expect_true(all(c("id", "variable", "value") %in% colnames(p$data)))
  # every cell of the input is represented
  expect_equal(nrow(p$data), nrow(score) * ncol(score))
  save_for_review(p, "Heatmap basic structure")
})

test_that("Heatmap keeps matrix row order in y axis (issue #15 pipeline)", {
  score <- heatmap_input()
  p <- Heatmap(score)
  # y axis (id) built from unique-reversed rownames of the input matrix
  expected <- rev(rownames(score))
  expect_equal(levels(p$data$id), expected)
})

# --- regression: #6 / #15 feature_text_subset -----------------------------

test_that("Heatmap feature_text_subset renders the requested genes (issue #6)", {
  score <- heatmap_input()
  p <- Heatmap(score, feature_text_subset = c("CD3D", "GAPDH"))
  expect_s3_class(p, "ggplot")
  # subsetting adds segment/text layers for the feature labels
  expect_gt(length(p$layers), 2)
  save_for_review(p, "Heatmap feature_text_subset (issue #6)")
})

test_that("Heatmap feature_text_subset crashes on a non-matching gene (issue #6)", {
  score <- heatmap_input()
  # KNOWN-ISSUE (HMP-1, P1, issue #6) -- source: R/Heatmap.R:223-229, error at :224
  #
  # This case used to be an unconditional `skip("pending fix: ...")`, which
  # meant the four assertions below it NEVER ran: the file reported a green
  # "issue #6 is pending" while actually asserting nothing at all.
  #
  # The checklist records the CURRENT behaviour as a hard error, not a warning:
  #   * :223 iterates `(nrow(text_table) - 1):1`. When the subset matches
  #     exactly one row that expression degenerates to `0:1`, so the loop body
  #     reads `text_table[0, "y1"]` -> "argument is of length zero".
  #   * a gene that matches nothing also collapses to a 1-row `text_table`
  #     (levels(ToPlot$id) %in% feature_text_subset -> character(0)), so both
  #     inputs land on the *same* error.
  # Intended fix per the checklist: `if (nrow > 1)` + `rev(seq_len(...))`.
  # Pinning the current behaviour keeps the case honest: the day the source is
  # fixed these assertions go red, which is the signal to re-file the checklist.
  expect_error(
    Heatmap(score, feature_text_subset = c("CD3D", "NOT_A_GENE")),
    "argument is of length zero"
  )
  # a subset matching EXACTLY ONE row hits the same `0:1` degeneration
  expect_error(
    Heatmap(score, feature_text_subset = "CD3D"),
    "argument is of length zero"
  )
  # ... while a subset matching >= 2 rows still renders (the working path)
  expect_s3_class(
    suppressWarnings(Heatmap(score, feature_text_subset = c("CD3D", "GAPDH"))),
    "ggplot"
  )
})

test_that("Heatmap feature_text_subset does not warn about unknown genes (issue #6)", {
  score <- heatmap_input()
  # KNOWN-ISSUE (HMP-5, P2) -- source: R/Heatmap.R:216-218
  # The checklist contrasts this with DotPlot2, which DOES warn and skip:
  #   DotPlot2(pbmc, features = c("CD3D", "NO_SUCH_GENE")) -> warning
  #   Heatmap(score, feature_text_subset = ...)            -> warning 0
  # Whatever this call does (error or render), it must not silently warn about
  # a feature name that does not exist -- so "no warning" is pinned here, and
  # the crashing path is pinned in the case above.
  expect_no_warning(try(Heatmap(score, feature_text_subset = c("CD3D", "NOT_A_GENE")),
                        silent = TRUE))
})

# --- parameter coverage (2.5) ---------------------------------------------

test_that("Heatmap color/scale parameters are accepted", {
  score <- heatmap_input()
  # color_scheme: Brewer diverging palette + custom gradient; col.min/max clip
  p <- Heatmap(score, color_scheme = "BuRd", center_color = TRUE,
               col.min = -1, col.max = 1)
  expect_s3_class(p, "ggplot")
  expect_lte(max(p$data$value), 1)
  expect_gte(min(p$data$value), -1)

  p2 <- Heatmap(score, color_scheme = c("white", "darkgreen"))
  expect_s3_class(p2, "ggplot")
  save_for_review(p, "Heatmap color scheme and col.min col.max")
})

test_that("Heatmap label/geometry parameters are accepted", {
  score <- heatmap_input()
  p <- Heatmap(
    score,
    lab_fill = "zscore",
    angle = 45, hjust = 1, vjust = 1,
    legend_position = "bottom",
    y_text_position = "right",
    border_color = "grey50",
    hide_axis_line = FALSE
  )
  expect_s3_class(p, "ggplot")
  save_for_review(p, "Heatmap labels angle border parameters")
})

test_that("Heatmap facet_col / facet_row split the heatmap", {
  score <- heatmap_input()   # 3 features x 9 clusters
  # facet_col must match the number of matrix COLUMNS (9)
  fc <- c(rep("group1", 4), rep("group2", 5))
  p <- Heatmap(score, facet_col = fc)
  expect_s3_class(p, "ggplot")
  # facet_row matches the number of ROWS (3)
  p2 <- Heatmap(score, facet_row = c("g1", "g1", "g2"))
  expect_s3_class(p2, "ggplot")
  save_for_review(p2, "Heatmap facet_row split")
})

test_that("Heatmap rejects wrong-length facet arguments", {
  score <- heatmap_input()
  expect_error(
    Heatmap(score, facet_col = c("a", "b")),
    regexp = "facet_col"
  )
  expect_error(
    Heatmap(score, facet_row = c("a", "b")),
    regexp = "facet_row"
  )
})

# --- input validation ------------------------------------------------------

test_that("Heatmap rejects a non-matrix-like input", {
  # no rownames/matrix structure -> downstream data.frame() fails loudly
  expect_error(Heatmap("not_a_matrix"))
})
