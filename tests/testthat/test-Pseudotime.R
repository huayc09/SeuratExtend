# Test file for R/Pseudotime.R + R/RunSlingshot.R + R/Loom.R (trajectory &
# data-bridge exports)
# Covers trajectory/backend exports. Per spec 2.7, backends that are absent
# (slingshot / python) make those tests skip; R-side behaviour that works
# without the backend (Seu2Loom via hdf5r) is asserted directly.
# Also covers the trajectory helper set in R/Pseudotime.R: check_method,
# rescale, compute_gradient, trim_data_by_pseudotime, extract_pseudotime_data,
# perform_gam_analysis, build_base_plot, add_plot_layer,
# process_magic_and_set_assay and the find_matching_reduction closure inside
# DimPlot.Slingshot.

trajectory_exports <- c(
  "RunSlingshot", "DimPlot.Slingshot", "GeneTrendCurve.Slingshot",
  "GeneTrendHeatmap.Slingshot", "GeneTrendCurve.Palantir",
  "GeneTrendHeatmap.Palantir", "Palantir.RunDM", "Palantir.Pseudotime",
  "Palantir.Magic", "RunPalantirDiffusionMap", "ImportPyscenicLoom",
  "Cellrank.Compute", "Cellrank.Plot", "scVelo.Plot", "scVelo.RunBasic",
  "scVelo.SeuratToAnndata", "scVelo.SeuratToLoom", "Seu2Adata", "Seu2Loom",
  "adata.AddDR", "adata.AddMetadata", "adata.Load", "adata.LoadLoom",
  "adata.Save", "activate_python", "create_condaenv_seuratextend"
)

test_that("trajectory/python exports are all present", {
  exports <- getNamespaceExports("SeuratExtend")
  missing <- setdiff(trajectory_exports, exports)
  expect_length(missing, 0)
})

# --- Seu2Loom: pure R via hdf5r (no python needed) ------------------------

test_that("Seu2Loom writes a loom file via hdf5r", {
  data(pbmc)
  skip_if_not(requireNamespace("hdf5r", quietly = TRUE))
  loompath <- file.path(tempdir(), "test_seu.loom")
  on.exit(unlink(loompath), add = TRUE)
  # hdf5r prints progress bars to stdout; not part of the assertion
  expect_no_error(
    capture.output(
      Seu2Loom(pbmc, filename = loompath, add.normdata = FALSE,
               overwrite = TRUE),
      type = "output"
    )
  )
  expect_true(file.exists(loompath))
})

test_that("Seu2Loom refuses to overwrite without overwrite=TRUE", {
  data(pbmc)
  skip_if_not(requireNamespace("hdf5r", quietly = TRUE))
  loompath <- file.path(tempdir(), "test_overwrite.loom")
  on.exit(unlink(loompath), add = TRUE)
  invisible(capture.output(
    Seu2Loom(pbmc, filename = loompath, overwrite = TRUE),
    type = "output"))
  expect_error(Seu2Loom(pbmc, filename = loompath),
               regexp = "already exists")
})

# --- slingshot-backed (skip when package absent) ---------------------------

test_that("RunSlingshot returns a Seurat object when slingshot installed", {
  skip_if_not(requireNamespace("slingshot", quietly = TRUE))
  data(pbmc)
  # group.by has no default, so it must be supplied. reducedDim must be given
  # in the case slingshot's SCE side expects ("UMAP"): the Seurat object stores
  # the reduction as "umap", and the shipped default "PCA" uses all 50
  # dimensions, which is computationally singular on pbmc. Both observations
  # are pinned in the L3-9 audit.
  out <- RunSlingshot(pbmc, group.by = "cluster", reducedDim = "UMAP")
  expect_s4_class(out, "Seurat")
})

test_that("DimPlot.Slingshot is exported and ggplot-based when runnable", {
  data(pbmc)
  skip_if_not(requireNamespace("slingshot", quietly = TRUE))
  sling <- RunSlingshot(pbmc, group.by = "cluster", reducedDim = "UMAP")
  # The curve layer is drawn with the deprecated `size` aesthetic (KNOWN-ISSUE,
  # pinned below in "draws curves with the deprecated size aesthetic"). That is
  # a lifecycle soft warning raised once per session, so whether it fires here
  # depends on test order -- swallow it here instead of asserting on it.
  p <- suppressWarnings(DimPlot.Slingshot(sling))
  expect_true(inherits(p, "ggplot"))
})

# --- python/conda-backed (skip on headless CI) -----------------------------

test_that("python-backed Seu2Adata requires python (issue #53 family)", {
  data(pbmc)
  skip_if(py_available_for_test())
  # no python -> must fail on the R side, not hang or crash the session
  expect_error(Seu2Adata(pbmc))
})

test_that("cellrank/scvelo/palantir functions are exported", {
  exports <- getNamespaceExports("SeuratExtend")
  for (fn in c("Cellrank.Compute", "Cellrank.Plot", "scVelo.Plot",
               "scVelo.RunBasic", "scVelo.SeuratToAnndata",
               "scVelo.SeuratToLoom", "Palantir.RunDM",
               "RunPalantirDiffusionMap")) {
    expect_true(fn %in% exports, info = fn)
  }
})

# ===========================================================================
# Trajectory helpers in R/Pseudotime.R. These are the internal parts the four
# GeneTrend* functions are assembled from, so pinning them pins the plots.
# ===========================================================================

# Fixture in the shape the helpers expect:
#   df      - one row per cell: Pseudotime + one column per fate (probability)
#   df_melt - long form; perform_gam_analysis() needs Gene/cell_fate/Pseudotime/
#             Expression, trim_data_by_pseudotime() needs cell_fate/fate_prob.
pt_fixture <- function(n = 20) {
  df <- data.frame(
    Pseudotime = seq(0, 1, length.out = n),
    fate1 = seq(0.9, 0.1, length.out = n),
    fate2 = seq(0.1, 0.9, length.out = n),
    row.names = paste0("cell", seq_len(n))
  )
  df_melt <- data.frame(
    Gene = "G1",
    cell_fate = rep(c("fate1", "fate2"), each = n),
    Pseudotime = rep(df$Pseudotime, 2),
    Expression = c(seq(1, 5, length.out = n), seq(5, 1, length.out = n)),
    fate_prob = c(df$fate1, df$fate2)
  )
  list(df = df, df_melt = df_melt)
}

# --- check_method ---------------------------------------------------------

test_that("check_method() accepts gam/smooth and rejects everything else", {
  # Source: R/Pseudotime.R:415-419.
  expect_no_error(check_method("gam"))
  expect_no_error(check_method("smooth"))
  expect_error(check_method("bogus"), "Invalid method")
  expect_error(check_method(character(0)), "Invalid method")
  expect_error(check_method(NA_character_), "Invalid method")
})

test_that("check_method() validates only the first element", {
  # KNOWN-ISSUE (L3-8 PST-8): the default is written as the enumerations
  # c("smooth", "gam"), but only method[1] is inspected, so trailing junk is
  # silently accepted. Source: R/Pseudotime.R:416.
  expect_no_error(check_method(c("smooth", "bogus")))
  expect_no_error(check_method(c("gam", "bogus")))
})

# --- rescale --------------------------------------------------------------

test_that("rescale() maps a vector onto [0, 1]", {
  # Source: R/Pseudotime.R:561-563.
  # NOTE: must stay namespaced. Package source calls library(scales) at run time
  # (R/Heatmap.R:130, R/runGSVA.R:44, R/WaterfallPlot.R:430); library() attaches
  # a package at search position 2, i.e. AHEAD of package:SeuratExtend. A bare
  # rescale() therefore resolves to scales::rescale -- an S3 generic whose
  # zero_range guard returns mean(to) = 0.5 for a constant input instead of the
  # NaN this implementation returns. That is exactly how the assertion below
  # went red in a full test_dir() run while passing when run alone.
  expect_equal(SeuratExtend:::rescale(c(1, 2, 3)), c(0, 0.5, 1))
  expect_equal(SeuratExtend:::rescale(c(10, 20)), c(0, 1))
  expect_equal(SeuratExtend:::rescale(c(-1, 1)), c(0, 1))
})

test_that("rescale() has no guard for constant or NA-containing input", {
  # KNOWN-ISSUE (L3-8 PST-4): min()/max() run without na.rm and the denominator
  # can be 0. A constant row -> NaN; one NA -> NA for the whole row. Both reach
  # the heatmap as a blank line without any warning.
  expect_true(all(is.nan(SeuratExtend:::rescale(c(5, 5, 5)))))
  expect_true(all(is.nan(SeuratExtend:::rescale(7))))
  expect_true(all(is.na(SeuratExtend:::rescale(c(1, NA, 3)))))
})

# --- compute_gradient -----------------------------------------------------

test_that("compute_gradient() averages the slope around the row peak", {
  # Source: R/Pseudotime.R:565-577.
  expect_equal(compute_gradient(1:20), 1)      # monotone -> mean(diff) == 1
  expect_equal(compute_gradient(20:1), -1)     # monotone down
  # plateau then jump: peak at index 6, window = 0..10 -> a single step of 1
  # spread over 9 differences
  expect_equal(compute_gradient(c(rep(0, 5), rep(1, 5))), 1 / 9,
               tolerance = 1e-8)
})

test_that("compute_gradient() clamps its window at both ends", {
  # left = max(1, peak - n), right = min(length(row), peak + n)
  expect_equal(compute_gradient(1:3, n = 10), 1)   # window covers everything
  expect_equal(compute_gradient(1:3, n = 1), 1)
})

# --- trim_data_by_pseudotime ----------------------------------------------

test_that("trim_data_by_pseudotime() drops low-probability tail cells", {
  # Source: R/Pseudotime.R:547-559.
  fx <- pt_fixture()
  expect_equal(nrow(fx$df_melt), 40)                       # 2 fates x 20 cells
  trimmed <- trim_data_by_pseudotime(fx$df, fx$df_melt, "fate1")
  # one start cell (Pseudotime <= the 5% quantile) -> min fate1 prob is 0.9
  # -> the 19 fate1 rows below 0.9 are dropped
  expect_equal(nrow(trimmed), 21)
  # everything else is left alone
  expect_equal(sum(trimmed$cell_fate == "fate2"), 20)
})

test_that("trim_data_by_pseudotime() crashes on an NA pseudotime", {
  # KNOWN-ISSUE (L3-8 PST-2): quantile() is called without na.rm, so a single
  # NA in Pseudotime raises a generic stats error that never mentions the
  # missing pseudotime. The Slingshot path of the same package na.omit()s the
  # same data and works, so the two trajectories disagree. Source: :549.
  fx <- pt_fixture()
  fx$df$Pseudotime[3] <- NA
  expect_error(trim_data_by_pseudotime(fx$df, fx$df_melt, "fate1"),
               "missing values and NaN's not allowed")
})

test_that("trim_data_by_pseudotime() fabricates all-NA rows when a fate has NA", {
  # KNOWN-ISSUE (L3-8 PST-3): min() over the start cells runs without na.rm.
  # One NA inside the 5% start window -> min is NA -> `fate_prob < NA` is NA ->
  # `!(...)` is NA -> logical-NA indexing of a data.frame, which invents all-NA
  # rows instead of dropping them. Measured here: 40 rows in, 20 of them all-NA.
  # Source: R/Pseudotime.R:554-555.
  fx <- pt_fixture()
  fx$df$fate1[1] <- NA          # row 1 is inside the start-cell window
  out <- trim_data_by_pseudotime(fx$df, fx$df_melt, "fate1")
  expect_equal(nrow(out), 40)
  expect_equal(sum(!stats::complete.cases(out)), 20)
})

test_that("trim_data_by_pseudotime() warns when a requested fate has no column", {
  # KNOWN-ISSUE (L3-8 PST-11): GeneTrendHeatmap.Palantir passes every
  # fate_name while its df only carries the selected fate, so min(NULL) is Inf
  # (+ a warning) and that fate's trim is a silent no-op.
  fx <- pt_fixture()
  df_one <- fx$df[, c("Pseudotime", "fate1"), drop = FALSE]
  expect_warning(
    out <- trim_data_by_pseudotime(df_one, fx$df_melt, c("fate1", "fate_gone")),
    "no non-missing arguments to min")
  expect_equal(nrow(out), 21)   # the missing fate trims nothing
})

# --- extract_pseudotime_data ----------------------------------------------

test_that("extract_pseudotime_data() validates a data.frame argument", {
  # Source: R/Pseudotime.R:531-539.
  data(pbmc)
  n <- ncol(pbmc)
  df <- data.frame(Pseudotime = seq(0, 1, length.out = n),
                   Entropy = seq(1, 0, length.out = n),
                   row.names = colnames(pbmc))
  out <- extract_pseudotime_data(pbmc, df, "Palantir",
                                 required_cols = c("Pseudotime", "Entropy"))
  expect_true(is.data.frame(out))
  expect_equal(nrow(out), n)

  # a missing required column is reported by name
  expect_error(
    extract_pseudotime_data(pbmc, df, "Palantir",
                            required_cols = c("Pseudotime", "Entropy", "Nope")),
    "does not contain the required columns")

  # row names must match the cells exactly
  bad <- df
  rownames(bad) <- paste0("zzz", seq_len(n))
  expect_error(extract_pseudotime_data(pbmc, bad, "Palantir"),
               "row names")

  # anything that is neither a character nor a data frame is rejected
  expect_error(extract_pseudotime_data(pbmc, 42, "Palantir"),
               "must either be a character")
})

test_that("extract_pseudotime_data() skips both checks on the character branch", {
  # KNOWN-ISSUE (L3-8 PST-1): the two branches validate to opposite strengths.
  # The data.frame branch (:531-539) checks required_cols AND the row names,
  # the character branch (:516-530) checks neither. Downstream cbind() then
  # binds the frame by POSITION, so a frame read from seu@misc in a different
  # order lands on the wrong cells -- silently, with no warning. Contrast: the
  # same frame passed directly is rejected outright.
  data(pbmc)
  n <- ncol(pbmc)
  df <- data.frame(Pseudotime = seq(0, 1, length.out = n),
                   row.names = paste0("zzz", seq_len(n)))   # wrong rows on purpose

  # direct -> refused
  expect_error(extract_pseudotime_data(pbmc, df, "Palantir"), "row names")

  # via seu@misc -> accepted, rows still wrong
  s <- pbmc
  s@misc$Palantir <- list(pt = df)
  out <- extract_pseudotime_data(s, "pt", "Palantir",
                                 required_cols = c("Pseudotime", "Entropy"))
  expect_equal(rownames(out), rownames(df))
  expect_false(identical(rownames(out), colnames(pbmc)))

  # the Slingshot branch behaves the same way
  s2 <- pbmc
  s2@misc$slingshot <- list(umap = list(SlingPseudotime = df))
  expect_equal(rownames(extract_pseudotime_data(s2, "umap", "Slingshot")),
               rownames(df))

  # and it reports a missing Slingshot entry by name
  s3 <- pbmc
  s3@misc$slingshot <- list(umap = list(other = 1))
  expect_error(extract_pseudotime_data(s3, "umap", "Slingshot"),
               "Slingshot pseudotime output not found")
})

# --- perform_gam_analysis -------------------------------------------------

test_that("perform_gam_analysis() returns a 200-point grid per gene/fate", {
  # Source: R/Pseudotime.R:430-473.
  skip_if_not(requireNamespace("mgcv", quietly = TRUE))
  fx <- pt_fixture()
  pred <- suppressMessages(perform_gam_analysis(fx$df_melt, "G1",
                                                c("fate1", "fate2")))
  expect_equal(nrow(pred), 400)          # 2 fates x 200 grid points
  expect_equal(sort(colnames(pred)),
               sort(c("Pseudotime", "Gene", "cell_fate", "Expression",
                      "Predicted")))
  expect_equal(nlevels(pred$Gene), 1)
})

test_that("perform_gam_analysis() collapses to one unhelpful error", {
  # KNOWN-ISSUE (L3-8): every per-gene/fate failure is caught, warned about and
  # then swallowed; if none fitted, the caller only sees a stop() that blames
  # "your data" although the real causes are "fewer than 4 usable points" or
  # "all pseudotime values identical". Source: :444-447, :450-462, :466-468.
  skip_if_not(requireNamespace("mgcv", quietly = TRUE))
  fx <- pt_fixture()

  too_few <- fx$df_melt[1:3, ]
  expect_warning(
    expect_error(perform_gam_analysis(too_few, "G1", "fate1"),
                 "GAM fitting failed for all gene/fate combinations"),
    "insufficient data points")

  flat <- fx$df_melt
  flat$Pseudotime <- 0.5
  expect_warning(
    expect_error(perform_gam_analysis(flat, "G1", "fate1"),
                 "GAM fitting failed for all gene/fate combinations"),
    "fewer unique covariate combinations")
})

# --- build_base_plot / add_plot_layer -------------------------------------

test_that("build_base_plot() always facets by Gene", {
  # Source: R/Pseudotime.R:475-494.
  fx <- pt_fixture()
  with_pts <- build_base_plot(TRUE, fx$df_melt, 0.5, 0.2, 1)
  expect_s3_class(with_pts, "ggplot")
  expect_length(with_pts$layers, 1)
  expect_equal(class(with_pts$layers[[1]]$geom)[1], "GeomPoint")

  no_pts <- build_base_plot(FALSE, fx$df_melt, 0.5, 0.2, 1)
  expect_length(no_pts$layers, 0)
  expect_true(inherits(no_pts$facet, "FacetWrap"))
  expect_equal(no_pts$labels$x, "Pseudotime")
})

test_that("add_plot_layer() switches geom on method and no-ops otherwise", {
  # Source: R/Pseudotime.R:496-510.
  fx <- pt_fixture()
  p0 <- build_base_plot(FALSE, fx$df_melt, 0.5, 0.2, 1)
  pred <- data.frame(Pseudotime = fx$df$Pseudotime,
                     Predicted = fx$df$fate1,
                     cell_fate = "fate1")

  p_gam <- add_plot_layer(p0, "gam", 1.5, FALSE, pred)
  expect_equal(class(p_gam$layers[[1]]$geom)[1], "GeomLine")
  expect_equal(p_gam$layers[[1]]$aes_params$linewidth, 1.5)

  p_sm <- add_plot_layer(p0, "smooth", 1.5, TRUE, fx$df_melt)
  expect_equal(class(p_sm$layers[[1]]$geom)[1], "GeomSmooth")

  # no default branch: an unknown method returns the plot untouched. This is
  # unreachable through the public API because check_method() runs first, so it
  # is an ℹ️ exclusion rather than a defect.
  expect_length(add_plot_layer(p0, "bogus", 1, TRUE, fx$df_melt)$layers, 0)
})

# --- process_magic_and_set_assay ------------------------------------------

test_that("process_magic_and_set_assay() keeps or switches the default assay", {
  # Source: R/Pseudotime.R:421-428.
  data(pbmc)
  expect_equal(DefaultAssay(process_magic_and_set_assay(pbmc, FALSE, NULL)),
               DefaultAssay(pbmc))
  # KNOWN-ISSUE (L3-8 PST-12): magic has no type check, so a non-logical value
  # dies inside R's && with a message that says nothing about magic.
  expect_error(process_magic_and_set_assay(pbmc, NULL, NULL),
               "invalid 'x' type")
})

# --- find_matching_reduction (closure inside DimPlot.Slingshot) -----------

test_that("DimPlot.Slingshot() matches the reduction name case-insensitively", {
  # The closure find_matching_reduction() (:679-688) tries an exact match and
  # then a toupper() match. RunSlingshot stores its result under the reducedDim
  # spelling the caller passed, so "umap" vs "UMAP" is the common case.
  skip_if_not(requireNamespace("slingshot", quietly = TRUE))
  data(pbmc)
  sling <- RunSlingshot(pbmc, group.by = "cluster", reducedDim = "UMAP")
  expect_true("UMAP" %in% names(sling@misc$slingshot))

  single <- DimPlot.Slingshot(sling, reduction = "umap", combine = FALSE)
  expect_true(is.list(single))
  expect_s3_class(single[[1]], "ggplot")
  # and the auto-detection route picks the same entry without a message
  expect_s3_class(DimPlot.Slingshot(sling), "ggplot")

  # nothing matches -> the error names what is available
  expect_error(DimPlot.Slingshot(sling, reduction = "nope"),
               "Available Slingshot reductions: UMAP")
  # no trajectory at all -> a different, actionable message
  expect_error(DimPlot.Slingshot(pbmc),
               "No Slingshot results found in seu@misc\\$slingshot")

  # the closure really is defined inside the function body
  body_txt <- paste(deparse(body(SeuratExtend:::DimPlot.Slingshot)),
                    collapse = "\n")
  expect_true(grepl("find_matching_reduction <- function", body_txt,
                    fixed = TRUE))
  expect_true(grepl("toupper", body_txt, fixed = TRUE))
})

test_that("DimPlot.Slingshot() draws curves with the deprecated size aesthetic", {
  # KNOWN-ISSUE (L3-8 PST-7): the curve layer is built with
  # geom_path(size = curve.size) while the rest of the same file uses
  # linewidth (:501, :506). ggplot2 >= 3.4 deprecates it -- measured warning:
  #   Using `size` aesthetic for lines was deprecated in ggplot2 3.4.0.
  # The warning is a lifecycle soft warning (raised once per session, and not
  # catchable from inside the construction frame), so it is pinned at the
  # object level instead: the value lands in aes_params under "size".
  # Source: R/Pseudotime.R:797-804.
  skip_if_not(file.exists(file.path("..", "..", "R", "Pseudotime.R")),
              "dev-layout source inspection only (not available under R CMD check)")
  src <- readLines(file.path("..", "..", "R", "Pseudotime.R"), warn = FALSE)
  expect_true(any(grepl("size = curve.size", src, fixed = TRUE)))
  expect_true(any(grepl("linewidth = line.size", src, fixed = TRUE)))

  skip_if_not(requireNamespace("slingshot", quietly = TRUE))
  data(pbmc)
  sling <- RunSlingshot(pbmc, group.by = "cluster", reducedDim = "UMAP")
  g <- DimPlot.Slingshot(sling, combine = FALSE)[[1]]
  last <- g$layers[[length(g$layers)]]
  expect_equal(class(last$geom)[1], "GeomPath")
  expect_true("size" %in% names(last$aes_params))
  expect_false("linewidth" %in% names(last$aes_params))
})

