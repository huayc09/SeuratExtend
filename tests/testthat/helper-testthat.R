# Shared test helpers (auto-sourced by testthat before test files run)

# Python availability probe (spec 2.7): python-integration tests only run
# when a usable Python environment is present, otherwise they are skipped.
py_available_for_test <- function() {
  requireNamespace("reticulate", quietly = TRUE) &&
    reticulate::py_available(initialize = FALSE)
}

# Test fixture helpers (small, explicit, no hidden state)
#
# `SeuratObject::subset` is NOT an exported generic -- calling it raises
# "'subset' is not an exported object from 'namespace:SeuratObject'". The
# helper was never actually called by any test file, so the typo stayed
# invisible. Use the plain generic (Seurat registers the S3 method).
subset_pbmc <- function(pbmc, clusters) {
  cells <- colnames(pbmc)[pbmc$cluster %in% clusters]
  subset(pbmc, cells = cells)
}

# A few package functions (import() in particular) write an interactive
# prompt straight to stdout even when the session is non-interactive, so
# every test that reaches them would otherwise spray "Try install? y/[n]"
# through the test log. quietly() swallows that stream while leaving the
# return value and any condition (message/warning/error) untouched.
quietly <- function(expr) {
  invisible(utils::capture.output(value <- expr, type = "output"))
  value
}

# Does reticulate see a conda installation? Used to skip the conda-only
# assertions on machines without one, instead of asserting on the shape of
# whatever error reticulate happens to raise.
conda_available <- function() {
  requireNamespace("reticulate", quietly = TRUE) &&
    !is.null(tryCatch(reticulate::conda_list(), error = function(e) NULL))
}

# A `node stack overflow` escapes testthat's own condition handling: the
# handler itself needs stack space to run, so expect_error() reports it as an
# uncaught error rather than a satisfied expectation. Several functions in this
# package recurse into themselves on non-vector input (RenameGO,
# RenameReactome, the GenesymbolConversion converters), so capture the error
# with an explicit tryCatch instead. Verified with R 4.6.1.
expect_stack_overflow <- function(expr) {
  msg <- tryCatch({ force(expr); "no condition raised" },
                  error = function(e) conditionMessage(e))
  expect_match(msg, "node stack overflow")
}

# Red-line guard: tests must not write files into the package working
# directory. Some plot building (ggplotGrob / cowplot::plot_grid used by
# DimPlot2) opens the *default* device when none is active, which in
# non-interactive R writes "Rplots.pdf" into the current directory. Point the
# default device at a NULL pdf device so every implicit plot is discarded and
# leaves no file behind (verified: with this guard active, plot() creates
# nothing; it does not even create a temporary file).
#
# Do NOT gate this on capabilities("pdf"). There is no such capability:
# capabilities() reports jpeg/png/tiff/cairo/... but no "pdf", so
# capabilities("pdf") returns a zero-length logical, isTRUE() is FALSE, and the
# whole guard is silently skipped. That is exactly how Rplots.pdf kept being
# written into tests/testthat/ on every run.
options(device = function(...) grDevices::pdf(NULL))
# Belt and braces: keep one NULL pdf device open for the session.
if (length(grDevices::dev.list()) == 0) {
  grDevices::pdf(NULL)
}
