# Test file for R/Functions.R
# Covers: feature_percent (percent-positive matrix), show_col2 (palette grid)
# Regression issues:
#   #5  - dplyr/Seurat method conflicts: grouped_df-style inputs must not
#         break the stat pipeline (original bug: VlnPlot2 + wilcox failed
#         with 'unable to find an inherited method for select/grouped_df')
#   #8  - show_col2 must be exported (was missing from NAMESPACE)
# Data: built-in pbmc

# --- regression: #8 export completeness -----------------------------------

test_that("show_col2 is exported (issue #8)", {
  expect_true("show_col2" %in% getNamespaceExports("SeuratExtend"))
})

# --- feature_percent semantics --------------------------------------------

test_that("feature_percent computes per-group positive-cell fractions", {
  data(pbmc)
  m <- feature_percent(pbmc, feature = c("CD3D", "MS4A1"),
                       group.by = "cluster")
  expect_true(is.matrix(m))
  # rows = features, cols = cluster levels
  expect_equal(rownames(m), c("CD3D", "MS4A1"))
  expect_equal(ncol(m), length(levels(pbmc$cluster)))
  # values are fractions in [0, 1]
  expect_true(all(m >= 0 & m <= 1))
  # CD3D is a T-cell marker: majority of CD4 T Naive positive, B cell low
  expect_gt(m["CD3D", "CD4 T Naive"], 0.5)
  expect_lt(m["CD3D", "B cell"], 0.2)
})

test_that("feature_percent parameter coverage (cells/above/total/ident/DefaultAssay)", {
  data(pbmc)
  # cells subset: only B cells
  b <- colnames(pbmc)[pbmc$cluster == "B cell"]
  m <- feature_percent(pbmc, feature = "MS4A1", group.by = "cluster",
                       cells = b)
  expect_true(is.matrix(m))

  # above: replaces the "expressed" cut-off (default 0). Raising it can only
  # shrink the numerator, so the fraction must be <= the looser cut-off.
  #
  # This assertion used to read `expect_true(all(m2 <= m2["CD3D", ] + 1))`,
  # which holds for every possible input (`x <= x + 1`), so it asserted
  # nothing at all. Measure it on a feature whose values actually straddle the
  # thresholds: CD3D happens to be identical at above = 0 and above = 0.5
  # (every positive value exceeds 0.5), so it cannot show the argument works.
  m1 <- feature_percent(pbmc, feature = "nCount_RNA", group.by = "cluster",
                        above = 1000)
  m2 <- feature_percent(pbmc, feature = "nCount_RNA", group.by = "cluster",
                        above = 5000)
  expect_true(all(m1["nCount_RNA", ] >= m2["nCount_RNA", ]))
  # ... and the argument is genuinely consumed (results are not identical)
  expect_false(identical(as.numeric(m1["nCount_RNA", ]),
                         as.numeric(m2["nCount_RNA", ])))

  # total = TRUE adds a "total" column
  m3 <- feature_percent(pbmc, feature = "CD3D", group.by = "cluster",
                        total = TRUE)
  expect_true("total" %in% colnames(m3))
  expect_equal(ncol(m3), length(levels(pbmc$cluster)) + 1)

  # if.expressed + min.pct -> logical matrix
  m4 <- feature_percent(pbmc, feature = "CD3D", group.by = "cluster",
                        if.expressed = TRUE, min.pct = 0.2)
  expect_true(is.logical(m4))

  # explicit DefaultAssay
  m5 <- feature_percent(pbmc, feature = "CD3D", group.by = "cluster",
                        DefaultAssay = "RNA")
  expect_true(is.matrix(m5))
})

test_that("feature_percent rejects invalid cells", {
  data(pbmc)
  expect_error(
    feature_percent(pbmc, feature = "CD3D", cells = "NOT_A_CELL"),
    regexp = "not found"
  )
})

# --- regression: #5 dplyr conflict (grouped_df) ---------------------------

test_that("VlnPlot2 stat pipeline tolerates grouped inputs (issue #5)", {
  data(pbmc)
  # Original issue #5: VlnPlot2 + stat.method failed on grouped_df via
  # dplyr::select dispatch. Regression: the call must succeed and render.
  p <- VlnPlot2(pbmc, features = c("CD3D", "MS4A1"), group.by = "cluster",
                stat.method = "wilcox.test")
  expect_s3_class(p, "ggplot")
})

# --- KNOWN-ISSUE pins (2026-09-16 recheck against 问题清单0915) ------------

test_that("feature_percent returns an all-ones row for a discrete feature", {
  data(pbmc)
  # KNOWN-ISSUE (FUN-3, P2) -- source: R/Functions.R:109 (FetchData) / :117
  # `feature` is documented as "anything that can be retrieved by FetchData",
  # which includes metadata columns. But apply() coerces the split data.frame
  # to a CHARACTER matrix, so `y > above` turns into a string comparison that
  # is always TRUE -> the "fraction" is exactly 1 in every group.
  # Zero error, zero warning: the caller gets a plausible-looking 1x9 matrix.
  m <- feature_percent(pbmc, feature = "cluster", group.by = "cluster")
  expect_true(is.matrix(m))
  expect_equal(dim(m), c(1L, length(levels(pbmc$cluster))))
  expect_true(all(m == 1))
})

test_that("feature_percent silently returns wrong values for a character 'above'", {
  data(pbmc)
  # KNOWN-ISSUE (FUN-4, P3) -- source: R/Functions.R:117 (`sum(y > above)`)
  # `above` is never type-checked, so a character value is compared
  # lexicographically and yields a WRONG but plausible number.
  # Measured on nCount_RNA (range 556..8854), first group (B cell):
  #   above = 2000   -> 0.3134
  #   above = "2000" -> 0.3731        no error, no warning
  m_num <- feature_percent(pbmc, feature = "nCount_RNA", group.by = "cluster",
                           above = 2000)
  m_chr <- feature_percent(pbmc, feature = "nCount_RNA", group.by = "cluster",
                           above = "2000")
  expect_false(identical(as.numeric(m_num), as.numeric(m_chr)))

  # min.pct has the same hole: a non-numeric value flips every cell to FALSE
  m_ok  <- feature_percent(pbmc, feature = "CD3D", group.by = "cluster",
                           if.expressed = TRUE, min.pct = 0.5)
  m_bad <- feature_percent(pbmc, feature = "CD3D", group.by = "cluster",
                           if.expressed = TRUE, min.pct = "x")
  expect_gt(sum(m_ok), 0)
  expect_equal(sum(m_bad), 0)
})

test_that("feature_percent reports an empty selection with a misleading message", {
  data(pbmc)
  # KNOWN-ISSUE (FUN-7, P3) -- source: R/Functions.R:83-84
  # `if (all(!cells %in% colnames(seu)))` fires for a ZERO-LENGTH vector
  # because all(logical(0)) is TRUE. An empty selection is therefore reported
  # as "'cells' not found in Seurat object", sending the user to look for a
  # typo in a name that was never supplied.
  expect_error(feature_percent(pbmc, feature = "CD3D", cells = character(0)),
               "not found in Seurat object")
})

test_that("show_col2 relies on symbols the package namespace cannot see", {
  # KNOWN-ISSUE (FUN-1, P1) -- source: R/Functions.R:212 (`%>%`, magrittr),
  # :213 (`list.rbind`, rlist), :220 (`NoAxes`, Seurat)
  #
  # Kept STATIC on purpose -- this case must NOT call show_col2(). In a fresh
  # session the call dies with `could not find function "%>%"`, and WHICH
  # symbol fails first depends on what the session happens to have attached,
  # so a behavioural assertion would flip green/red purely on file order.
  body_txt <- paste(deparse(body(SeuratExtend::show_col2)), collapse = "\n")
  expect_true(grepl("%>%", body_txt, fixed = TRUE))
  expect_true(grepl("list.rbind", body_txt, fixed = TRUE))

  # magrittr and rlist ARE in DESCRIPTION:Imports (so they get installed), but
  # NAMESPACE imports neither -> the namespaces are not on the search path of
  # the package namespace, which is why the call above cannot resolve them.
  ns_path <- file.path("..", "..", "NAMESPACE")
  skip_if_not(file.exists(ns_path), "NAMESPACE not reachable")
  ns <- readLines(ns_path, warn = FALSE)
  expect_false(any(grepl("magrittr", ns, fixed = TRUE)))
  expect_false(any(grepl("rlist", ns, fixed = TRUE)))

  desc_path <- file.path("..", "..", "DESCRIPTION")
  skip_if_not(file.exists(desc_path), "DESCRIPTION not reachable")
  desc_txt <- paste(readLines(desc_path, warn = FALSE), collapse = " ")
  expect_true(grepl("magrittr", desc_txt, fixed = TRUE))
  expect_true(grepl("rlist", desc_txt, fixed = TRUE))
})
