# Test file for R/WaterfallPlot.R (+ generics.R S3 dispatch)
# Covers: WaterfallPlot (Seurat + default), WaterfallPlot_v3 (Seurat alias)
# Data: built-in pbmc

# --- Seurat method ---------------------------------------------------------

test_that("WaterfallPlot returns a ggplot for marker features", {
  data(pbmc)
  # WaterfallPlot.Seurat requires 'features'; ident.1/ident.2 are optional
  p <- WaterfallPlot(pbmc, features = c("CD3D", "MS4A1", "CD79A", "LYZ"),
                     ident.1 = "B cell", ident.2 = "CD8 T cell")
  expect_true(inherits(p, "ggplot"))
  save_for_review(p, "WaterfallPlot B cell vs CD8 T cell")
})

test_that("WaterfallPlot accepts length/color parameter variants", {
  data(pbmc)
  p <- WaterfallPlot(pbmc, features = c("CD3D", "MS4A1"),
                     length = "tscore", color = "logFC",
                     flip = TRUE, style = "segment")
  expect_true(inherits(p, "ggplot"))
  save_for_review(p, "WaterfallPlot tscore segment flip")
})

test_that("WaterfallPlot's segment style uses the deprecated size aesthetic", {
  # KNOWN-ISSUE (WF-8, P3) -- R/WaterfallPlot.R:484
  #   geom_segment(aes(xend = rank, y = 0, yend = length), size = 0.8)
  # `size` on a line geom is deprecated as of ggplot2 3.4.0 (use `linewidth`),
  # and ggplot2 prints the deprecation pointing at SeuratExtend. The same row
  # is in the 问题清单0915 for this file.
  #
  # Checked STATICALLY on purpose: rlang reports this deprecation once per
  # session, so a runtime expect_warning() here is order-dependent. It passes
  # when this file runs alone and FAILS in the full suite, because
  # test-Pseudotime.R (alphabetically earlier) has already triggered the same
  # warning and it is never repeated. First attempt asserted the warning and
  # did exactly that. Same reasoning as the static pins in test-FeaturePlot.R
  # and test-Colors.R.
  src_path <- file.path("..", "..", "R", "WaterfallPlot.R")
  skip_if_not(file.exists(src_path), "source tree not reachable")
  src <- readLines(src_path, warn = FALSE)
  expect_true(any(grepl("geom_segment.*size\\s*=", src)))
})

test_that("WaterfallPlot_v3 aliases the Seurat method", {
  data(pbmc)
  p <- WaterfallPlot_v3(pbmc, features = c("CD3D", "MS4A1"))
  expect_true(inherits(p, "ggplot"))
})

# --- default (matrix) method ----------------------------------------------

test_that("WaterfallPlot.default accepts matrix + factor", {
  set.seed(42)
  matr <- matrix(abs(rnorm(240)), nrow = 4,
                 dimnames = list(c("g1", "g2", "g3", "g4"), NULL))
  f <- factor(rep(c("a", "b"), each = 30))
  p <- WaterfallPlot.default(matr, f = f)
  expect_true(inherits(p, "ggplot"))
})

# --- error handling ---------------------------------------------------------

test_that("WaterfallPlot generic requires features", {
  data(pbmc)
  # features has no default -> error when omitted
  expect_error(WaterfallPlot(pbmc), regexp = "features")
})

# --- KNOWN-ISSUE pins (2026-09-16 recheck against 问题清单0915) ------------

test_that("WaterfallPlot top.n always selects by score, order only arranges", {
  data(pbmc)
  # FIXED (WF-1) in work order runs/2026-10-22-deps-guard-silentfix:
  # R/WaterfallPlot.R used to take head(tail) of the CURRENT row order, so with
  # order = FALSE the selected bars were decided by the caller's writing order.
  # Selection is now always by value (largest top.n positives + smallest
  # top.n negatives) and `order` only controls the final arrangement.
  feat40 <- rownames(pbmc)[1:40]
  dat <- function(f, ord) WaterfallPlot(
    pbmc, features = f, ident.1 = "B cell", top.n = 3, order = ord)$data

  full <- WaterfallPlot(pbmc, features = feat40, ident.1 = "B cell")$data
  pos <- full[full$length > 0, , drop = FALSE]
  neg <- full[full$length <= 0, , drop = FALSE]
  expected <- c(rownames(pos)[order(pos$length, decreasing = TRUE)][1:3],
                rownames(neg)[order(neg$length)][1:3])

  d_nat <- dat(feat40, FALSE)
  d_rev <- dat(rev(feat40), FALSE)
  d_ord <- dat(feat40, TRUE)

  # same selected set regardless of order flag or caller row order
  expect_setequal(rownames(d_nat), expected)
  expect_setequal(rownames(d_rev), expected)
  expect_setequal(rownames(d_ord), expected)

  # order = FALSE preserves the caller's row order among the selected bars
  expect_identical(rownames(d_nat), feat40[feat40 %in% rownames(d_nat)])
  # order = TRUE is ranked by score (descending)
  expect_false(is.unsorted(rev(d_ord$length)))
})

test_that("WaterfallPlot reports a bad identity with an unrelated error", {
  data(pbmc)
  # KNOWN-ISSUE (WF-3, P2) -- :312-313 do not check the identities, so an
  # all-NA row reaches :441 max(nchar(...)) and :445 `if (NA <= 2)` blows up.
  # The message never says "that group does not exist".
  expect_error(
    WaterfallPlot(pbmc, features = c("CD3D", "MS4A1"),
                  ident.1 = "NOPE1", ident.2 = "NOPE2"),
    "missing value where TRUE/FALSE needed"
  )
})

test_that("WaterfallPlot does not validate length / color", {
  data(pbmc)
  # KNOWN-ISSUE (WF-5, P2) -- :400-401 (root cause :397-398). `style` is put
  # through match.arg while `length` and `color` go straight into the
  # filtering expression, so a typo surfaces as a maths error.
  base_args <- list(pbmc, features = c("CD3D", "MS4A1"),
                    ident.1 = "B cell", ident.2 = "CD8 T cell")
  for (arg in c("length", "color")) {
    a <- base_args
    a[[arg]] <- "bogus"
    expect_error(do.call(WaterfallPlot, a), "non-numeric argument")
  }
})

test_that("WaterfallPlot drops unknown features after Seurat's warning", {
  data(pbmc)
  # KNOWN-ISSUE (WF-6, P2) -- the lookup lives in Seu2Matr (R/Internal.R:53-57),
  # which calls Seurat::FetchData(object = seu, vars = features, ...). Seurat
  # warns and drops the name it cannot find. So the caller IS told -- but:
  #   * the wording names only the variable, never WaterfallPlot nor the
  #     `features` argument: "The following requested variables were not
  #     found: ...", and
  #   * it is Seurat's warning, so it is identical for every SeuratExtend
  #     function that routes through Seu2Matr (VlnPlot2, VolcanoPlot, CalcStats,
  #     GeneSetAnalysis, runGSVA) -- it cannot be tuned per function.
  # The 问题清单0915 row for this item words it as "静默忽略、零警告", which is
  # NOT what happens: the warning is real (measured 2026-09-16 with a
  # call-stack capture; the stack runs
  #   WaterfallPlot.Seurat -> Seu2Matr -> FetchData.Seurat -> warn()).
  # An earlier version of this case read $data only, saw nrow == 1, and
  # concluded "silently" without ever measuring the conditions.
  # The contrast with DotPlot2 still holds: DotPlot2 has its own
  # validate_features() (R/DotPlot.R:568) and words its own warning.
  warns <- character(0)
  d <- withCallingHandlers(
    WaterfallPlot(pbmc, features = c("CD3D", "NOT_A_GENE_XYZ"),
                  ident.1 = "B cell")$data,
    warning = function(w) {
      warns <<- c(warns, conditionMessage(w)); invokeRestart("muffleWarning")
    },
    message = function(m) invokeRestart("muffleMessage")
  )
  expect_equal(nrow(d), 1)                       # only CD3D survives
  expect_true(any(grepl("requested variables were not found", warns)))
  expect_true(any(grepl("NOT_A_GENE_XYZ", warns, fixed = TRUE)))
  # and the warning never says which function or argument is at fault
  expect_false(any(grepl("features", warns, fixed = TRUE)))
})

test_that("WaterfallPlot rewrites the global random stream when sampling", {
  data(pbmc)
  # KNOWN-ISSUE (WF-7, P2) -- :281 calls set.seed() when max.cells.per.ident
  # is supplied and never restores .Random.seed, so a caller's downstream
  # sampling becomes irreproducible through no action of their own.
  kmers <- c("CD3D", "MS4A1", "CD79A", "LYZ")
  set.seed(1)
  baseline <- .Random.seed
  invisible(WaterfallPlot(pbmc, features = kmers, ident.1 = "B cell",
                          ident.2 = "CD8 T cell"))
  expect_identical(.Random.seed, baseline)          # untouched by default

  invisible(WaterfallPlot(pbmc, features = kmers, ident.1 = "B cell",
                          ident.2 = "CD8 T cell", max.cells.per.ident = 20))
  expect_false(identical(.Random.seed, baseline))   # silently rewritten

  # do not leak the pollution into the rest of the suite
  assign(".Random.seed", baseline, envir = globalenv())
})
