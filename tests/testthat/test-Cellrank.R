# Test file for R/Cellrank.R
#
# Cellrank.R has two exports - Cellrank.Compute() and Cellrank.Plot() - and
# both drive Python through reticulate, so on a machine without the
# `seuratextend` conda environment the interesting assertions are the R-side
# ones: the signatures, and the fact that both functions fail loudly instead
# of hanging or leaving the session broken.
#
# The quoting mechanism Cellrank.Plot() relies on (shQuote(type = "sh")
# embedded in a glue() string) is exercised directly in test-Python.R via
# r_vector_to_py(); see the L3-5 notes there.
#
# STYLE: characterisation tests - see test-Internal.R for the convention.

cellrank_exports <- c("Cellrank.Compute", "Cellrank.Plot")

test_that("both Cellrank entry points are exported", {
  exports <- getNamespaceExports("SeuratExtend")
  expect_length(setdiff(cellrank_exports, exports), 0)
})

test_that("Cellrank.Compute() signature is stable", {
  f <- formals(SeuratExtend::Cellrank.Compute)
  expect_identical(names(f), c("load.adata", "time_key", "conda_env"))
  expect_null(eval(f$load.adata))
  expect_identical(f$conda_env, "seuratextend")
  # time_key is required on purpose: there is no sensible default.
  # A missing formal has an empty symbol as its value, which deparses to "".
  expect_false(nzchar(deparse(f$time_key)))
})

test_that("Cellrank.Plot() signature is stable", {
  f <- formals(SeuratExtend::Cellrank.Plot)
  expect_identical(
    names(f),
    c("load.adata", "basis", "color", "groups", "palette", "alpha", "dpi",
      "legend_fontsize", "figsize", "xlim", "ylim", "save", "conda_env")
  )
  expect_null(eval(f$load.adata))
  expect_identical(eval(f$alpha), 0.15)
  expect_identical(eval(f$dpi), 300)
  expect_identical(eval(f$legend_fontsize), 9)
  expect_identical(eval(f$figsize), c(7, 5))
  expect_identical(f$conda_env, "seuratextend")
  expect_false(nzchar(deparse(f$basis)))   # basis is required
})

test_that("Cellrank.Compute() refuses to run without a time_key", {
  expect_error(suppressMessages(Cellrank.Compute()))
})

test_that("Cellrank.Plot() fails loudly when the Python stack is unusable", {
  # Environment-independent by construction: with no conda binary the failure
  # is raised by load_condaenv() ("Unable to find conda binary"); with conda
  # present but the environment/packages missing it moves on to reticulate and
  # fails there (no `adata` / no cellrank). Either way it must raise an R
  # error rather than silently returning or wedging the session.
  # Baseline recorded on: R 4.6.1 / Windows / no conda visible to reticulate.
  expect_error(suppressMessages(Cellrank.Plot(basis = "ms")))
})

test_that("Cellrank.Plot() builds its Python call from shQuote'd arguments", {
  # KNOWN-ISSUE (L3-5, L3-7 family): argument values are pasted into Python
  # source with shQuote(type = "sh"), so a value with a backslash and no
  # apostrophe reaches Python with its backslashes intact. Source:
  # R/Cellrank.R:116-134. The behaviour of shQuote itself is pinned in
  # test-Python.R; here we only pin that this function is one of its callers.
  body_txt <- deparse(body(SeuratExtend::Cellrank.Plot))
  expect_true(any(grepl("shQuote", body_txt, fixed = TRUE)))
  expect_true(any(grepl("py_run_string", body_txt, fixed = TRUE)))
})
