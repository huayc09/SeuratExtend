# Test file for the `external` matrix parameter (work order 2026-10-21).
#
# Contract: runs/2026-10-21-external-matrix-param/01-goal.md (v3).
#   F1  matrix rows are never trimmed - all rows participate
#   F2  features = NULL + external -> plot E only (no rownames(seu))
#   F3  features given + external -> mention order first, unmentioned E appended
#   F4  no external -> existing behaviour unchanged
#   cell-axis 4-tier alignment, rowname risk rules, sparse == dense.

ext_matrix <- function(seu, rows = c("pathA", "pathB", "pathC"), seed = 1) {
  set.seed(seed)
  matrix(
    rnorm(length(rows) * ncol(seu)),
    nrow = length(rows),
    dimnames = list(rows, colnames(seu))
  )
}

# --- 1. mixed selection values -------------------------------------------

test_that("mixed selection: matrix variable values match the external matrix per group", {
  data(pbmc)
  m <- ext_matrix(pbmc)
  cs <- CalcStats(
    pbmc,
    features = c("CD3D", "pathA"),
    external = m,
    group.by = "cluster",
    method = "mean"
  )
  expect_true(all(c("CD3D", "pathA") %in% rownames(cs)))
  for (cl in unique(as.character(pbmc$cluster))) {
    cells <- colnames(pbmc)[pbmc$cluster == cl]
    expect_equal(
      unname(unlist(cs["pathA", cl])),
      mean(m["pathA", cells]),
      tolerance = 1e-12
    )
  }
})

# --- 2. default append ----------------------------------------------------

test_that("features plus external appends every matrix row in matrix order", {
  data(pbmc)
  m <- ext_matrix(pbmc)
  cs <- CalcStats(pbmc, features = "CD3D", external = m,
                  group.by = "cluster", method = "mean")
  expect_equal(rownames(cs), c("CD3D", "pathA", "pathB", "pathC"))

  # F2: pure matrix selection never adds rownames(seu)
  cs2 <- CalcStats(pbmc, external = m, group.by = "cluster", method = "mean")
  expect_equal(rownames(cs2), rownames(m))
})

# --- 3. pure matrix equals the .default path ------------------------------

test_that("pure external selection equals CalcStats on the matrix itself", {
  data(pbmc)
  m <- ext_matrix(pbmc)
  a <- CalcStats(pbmc, external = m, group.by = "cluster", method = "mean")
  b <- CalcStats(m, f = pbmc$cluster, method = "mean")
  expect_equal(a, b)
})

# --- 4. order control -----------------------------------------------------

test_that("features controls order and unmentioned matrix rows are appended", {
  data(pbmc)
  m <- ext_matrix(pbmc)
  cs <- CalcStats(pbmc, features = c("pathC", "CD3D"), external = m,
                  group.by = "cluster", method = "mean")
  expect_equal(rownames(cs), c("pathC", "CD3D", "pathA", "pathB"))

  # duplicate features are deduplicated keeping the first occurrence
  std <- SeuratExtend:::Seu2Matr(
    pbmc,
    features = c("CD3D", "CD3D", "pathA"),
    external = m
  )
  expect_equal(colnames(std$matr), c("CD3D", "pathA", "pathB", "pathC"))
})

# --- 5. risk rules --------------------------------------------------------

test_that("duplicated external rownames are rejected with the duplicate list", {
  data(pbmc)
  m <- ext_matrix(pbmc)
  rownames(m) <- c("pathA", "pathA", "pathB")
  expect_error(SeuratExtend:::.external_validate(m, colnames(pbmc)),
               "pathA")
  expect_error(CalcStats(pbmc, external = m, group.by = "cluster"),
               "pathA")
})

test_that("cell-axis tier 1: identical column names pass silently", {
  data(pbmc)
  m <- ext_matrix(pbmc)
  expect_silent(SeuratExtend:::.external_validate(m, colnames(pbmc)))
})

test_that("cell-axis tier 2: same cell set in a different order is reordered", {
  data(pbmc)
  m <- ext_matrix(pbmc)
  mr <- m[, rev(colnames(pbmc)), drop = FALSE]
  v <- NULL
  expect_message(v <- SeuratExtend:::.external_validate(mr, colnames(pbmc)),
                 "reorder")
  expect_equal(colnames(v), colnames(pbmc))
  expect_equal(unname(v["pathA", ]), unname(m["pathA", ]))
})

test_that("cell-axis tier 3: unnamed columns of matching count align positionally with a message", {
  data(pbmc)
  m <- ext_matrix(pbmc)
  colnames(m) <- NULL
  v <- NULL
  expect_message(v <- SeuratExtend:::.external_validate(m, colnames(pbmc)),
                 "no cell names")
  expect_equal(colnames(v), colnames(pbmc))
  expect_equal(unname(v["pathA", ]), unname(m["pathA", ]))
})

test_that("cell-axis tier 4: partial name match or count mismatch is an error", {
  data(pbmc)
  m <- ext_matrix(pbmc)
  expect_error(SeuratExtend:::.external_validate(m[, -1, drop = FALSE],
                                                 colnames(pbmc)),
               "missing from 'external'")
  m2 <- m
  colnames(m2) <- NULL
  expect_error(SeuratExtend:::.external_validate(m2[, -1, drop = FALSE],
                                                 colnames(pbmc)),
               "does not match")
})

test_that("conflicting names warn and external values win", {
  data(pbmc)
  m <- ext_matrix(pbmc)
  rownames(m)[2] <- "CD3D"
  expect_warning(
    cs <- CalcStats(pbmc, features = c("CD3D", "pathA"), external = m,
                    group.by = "cluster", method = "mean"),
    "CD3D"
  )
  for (cl in unique(as.character(pbmc$cluster))) {
    cells <- colnames(pbmc)[pbmc$cluster == cl]
    expect_equal(
      unname(unlist(cs["CD3D", cl])),
      mean(m["CD3D", cells]),
      tolerance = 1e-12
    )
  }
})

test_that("validator rejects non-numeric and unnamed inputs", {
  data(pbmc)
  expect_error(SeuratExtend:::.external_validate(1:10, colnames(pbmc)),
               "matrix")
  df <- data.frame(a = 1:3, b = letters[1:3])
  expect_error(SeuratExtend:::.external_validate(df, colnames(pbmc)[1:3]),
               "b")
  m <- matrix(1:3, nrow = 3, dimnames = list(NULL, colnames(pbmc)[1]))
  expect_error(SeuratExtend:::.external_validate(m, colnames(pbmc)[1]),
               "row names")
  empty <- matrix(numeric(0), 0, 2)
  expect_error(SeuratExtend:::.external_validate(empty, colnames(pbmc)[1:2]))
})

# --- 6. cells subset alignment -------------------------------------------

test_that("cells subset by name or logical aligns external columns to the same cells", {
  data(pbmc)
  m <- ext_matrix(pbmc)
  cells <- colnames(pbmc)[seq(1, ncol(pbmc), by = 7)]
  by_name <- CalcStats(pbmc, external = m, group.by = "cluster",
                       cells = cells, method = "mean")
  by_logical <- CalcStats(pbmc, external = m, group.by = "cluster",
                          cells = colnames(pbmc) %in% cells, method = "mean")
  expect_equal(by_name, by_logical)
  for (cl in colnames(by_name)) {
    sub <- cells[as.character(pbmc$cluster[cells]) == cl]
    if (length(sub) == 0) next
    expect_equal(
      unname(unlist(by_name["pathB", cl])),
      mean(m["pathB", sub]),
      tolerance = 1e-12
    )
  }
})

# --- 7. sparse equals dense ----------------------------------------------

test_that("a dgCMatrix external gives identical results to its dense equivalent", {
  skip_if_not_installed("Matrix")
  data(pbmc)
  m <- ext_matrix(pbmc)
  sm <- Matrix::Matrix(m, sparse = TRUE)
  dense <- CalcStats(pbmc, features = c("CD3D", "pathA"), external = m,
                     group.by = "cluster", method = "mean")
  sparse <- CalcStats(pbmc, features = c("CD3D", "pathA"), external = sm,
                      group.by = "cluster", method = "mean")
  expect_equal(dense, sparse)
  expect_equal(
    SeuratExtend:::.external_validate(sm, colnames(pbmc)),
    SeuratExtend:::.external_validate(m, colnames(pbmc))
  )
})

# --- 8. plotting smoke ----------------------------------------------------

test_that("plot functions accept external and return ggplot objects", {
  data(pbmc)
  m <- ext_matrix(pbmc)

  expect_s3_class(
    VlnPlot2(pbmc, features = c("CD3D", "pathA"), external = m),
    "ggplot"
  )
  expect_s3_class(
    VlnPlot2(pbmc, external = m),
    "ggplot"
  )
  expect_s3_class(
    WaterfallPlot(pbmc, features = c("CD3D", "pathA"), external = m),
    "ggplot"
  )
  expect_s3_class(
    VolcanoPlot(pbmc, features = c("CD3D", "pathA"), external = m),
    "ggplot"
  )
  expect_s3_class(
    DimPlot2(pbmc, external = m),
    "ggplot"
  )
  expect_s3_class(
    DimPlot2(pbmc, features = c("CD3D", "pathA"), external = m),
    "ggplot"
  )
})

test_that("DimPlot2_GetData resolves external variables without the Idents fallback", {
  data(pbmc)
  m <- ext_matrix(pbmc)
  d <- DimPlot2_GetData(pbmc, external = m)
  expect_equal(colnames(d$data.var), rownames(m))
  expect_equal(unname(d$data.var[["pathA"]]), unname(m["pathA", ]))

  dm <- DimPlot2_GetData(pbmc, features = c("pathC", "CD3D"), external = m)
  expect_equal(colnames(dm$data.var), c("pathC", "CD3D", "pathA", "pathB"))
  expect_equal(unname(dm$data.var[["pathC"]]), unname(m["pathC", ]))
})

# --- 9. no external: unchanged "features required" -----------------------

test_that("features = NULL with external = NULL still errors", {
  data(pbmc)
  expect_error(VlnPlot2(pbmc), "features")
  expect_error(CalcStats(pbmc), "features")
  expect_error(WaterfallPlot(pbmc), "features")
})

# --- 10. GSEAplot untouched ----------------------------------------------

test_that("GSEAplot has no external formal (scope freeze)", {
  expect_false("external" %in% names(formals(GSEAplot)))
})

# --- repair round r1: DimPlot2 companion panels --------------------------

test_that("DimPlot2_GetData keeps group.by and cells.highlight with external", {
  data(pbmc)
  m <- ext_matrix(pbmc)

  d <- DimPlot2_GetData(pbmc, features = c("CD3D", "pathA"),
                        group.by = "cluster", external = m)
  expect_equal(colnames(d$data.var),
               c("CD3D", "pathA", "pathB", "pathC", "cluster"))
  expect_equal(as.character(d$data.var$cluster), as.character(pbmc$cluster))

  dh <- DimPlot2_GetData(pbmc, cells.highlight = colnames(pbmc)[1:10],
                         external = m)
  expect_equal(colnames(dh$data.var), c("pathA", "pathB", "pathC",
                                        "Selected_cells"))
  expect_equal(sum(dh$data.var$Selected_cells == "Selected"), 10)

  expect_s3_class(
    DimPlot2(pbmc, features = "CD3D", group.by = "cluster", external = m),
    "ggplot"
  )
})

# --- repair round r1: numeric correctness at the plot layer -------------

test_that("VlnPlot2 external values match the source matrix per cell group", {
  data(pbmc)
  m <- ext_matrix(pbmc)
  p <- VlnPlot2(pbmc, features = c("CD3D", "pathA"), external = m)
  d <- p$data
  expect_equal(sort(d$value[d$feature == "pathA"]),
               sort(unname(m["pathA", ])))
  for (grp in unique(as.character(pbmc$cluster))) {
    cells <- colnames(pbmc)[pbmc$cluster == grp]
    expect_equal(
      mean(d$value[d$feature == "pathA" & as.character(d$f) == grp]),
      mean(m["pathA", cells])
    )
  }
})

test_that("VolcanoPlot with external only uses matrix variables (F2 exception)", {
  data(pbmc)
  m <- ext_matrix(pbmc)
  p <- VolcanoPlot(pbmc, external = m)
  expect_s3_class(p, "ggplot")
  # points are exactly the matrix rows, not rownames(seu)
  expect_equal(nrow(p$data), nrow(m))
  expect_equal(rownames(p$data), rownames(m))
})

test_that("pure-matrix WaterfallPlot and VolcanoPlot smoke", {
  data(pbmc)
  m <- ext_matrix(pbmc)
  expect_s3_class(WaterfallPlot(pbmc, external = m), "ggplot")
  expect_s3_class(VolcanoPlot(pbmc, external = m), "ggplot")
})

# --- repair round r1: conflict warning in DimPlot2 -----------------------

test_that("DimPlot2 external conflict warns and keeps external values", {
  data(pbmc)
  m <- ext_matrix(pbmc)
  rownames(m)[2] <- "CD3D"
  expect_warning(
    d <- DimPlot2_GetData(pbmc, external = m),
    "CD3D"
  )
  expect_equal(unname(d$data.var[["CD3D"]]), unname(m["CD3D", ]))
})

# --- repair round r1: data.frame input and cell order --------------------

test_that("numeric data.frame external equals matrix input", {
  data(pbmc)
  m <- ext_matrix(pbmc)
  a <- CalcStats(pbmc, features = c("CD3D", "pathA"), external = m,
                 group.by = "cluster", method = "mean")
  b <- CalcStats(pbmc, features = c("CD3D", "pathA"),
                 external = as.data.frame(m),
                 group.by = "cluster", method = "mean")
  expect_equal(a, b)
})

test_that("reversed/subset cells align external values by cell name", {
  data(pbmc)
  m <- ext_matrix(pbmc)
  cells <- rev(colnames(pbmc)[seq(1, ncol(pbmc), by = 7)])
  std <- SeuratExtend:::Seu2Matr(pbmc, external = m, cells = cells)
  expect_equal(rownames(std$matr), cells)
  expect_equal(unname(std$matr[, "pathA"]), unname(m["pathA", cells]))

  d <- DimPlot2_GetData(pbmc, external = m, cells = cells)
  expect_equal(rownames(d$data.var), cells)
  expect_equal(unname(d$data.var[["pathB"]]), unname(m["pathB", cells]))
})

# --- repair round r1: input risk rules -----------------------------------

test_that("duplicated external column names and duplicated cells are rejected", {
  data(pbmc)
  m <- ext_matrix(pbmc)
  md <- m
  colnames(md)[2] <- colnames(md)[1]
  expect_error(SeuratExtend:::.external_validate(md, colnames(pbmc)),
               "duplicated column names")
  dup_cells <- rep(colnames(pbmc)[1], 2)
  expect_error(SeuratExtend:::Seu2Matr(pbmc, external = m, cells = dup_cells),
               "duplicated cell names")
})

test_that("tier-1 identical colnames emit no message through Seu2Matr", {
  data(pbmc)
  m <- ext_matrix(pbmc)
  expect_silent(SeuratExtend:::Seu2Matr(pbmc, external = m))
})
