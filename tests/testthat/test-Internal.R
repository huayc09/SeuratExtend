# Test file for R/Internal.R
#
# Internal.R exports nothing: every symbol in it is a helper reached by other
# files in the package, so these tests go through `:::`.
#
# SCOPE: Seu2Matr() input validation, CheckIdent() contrast resolution,
# `%||%` / is_empty() semantics, check_spe() / check_pkg_version() guards,
# message_only_once() once-per-session behaviour, the SeuratObject v4/v5
# GetAssayData shim, and MinMax_internal() clipping.
#
# STYLE: these are *characterisation* tests - they pin the behaviour the
# shipped source actually has, so that a refactor which changes it goes red.
# Where the pinned behaviour is a known defect it is tagged
# `# KNOWN-ISSUE <id>` together with the source line, and the id refers to the
# audit checklist for this package.

internal_helpers <- c(
  "Seu2Matr", "CheckIdent", "%||%", "is_empty", "check_spe",
  "check_pkg_version", "message_only_once", ".get_assay_data_compat",
  "MinMax_internal"
)

test_that("Internal.R helpers are internal, not exported", {
  exports <- getNamespaceExports("SeuratExtend")
  leaked <- intersect(internal_helpers, exports)
  expect_length(leaked, 0)
  # ...and they really do exist in the namespace
  ns <- asNamespace("SeuratExtend")
  expect_true(all(internal_helpers %in% ls(ns, all.names = TRUE)))
})

# --- `%||%` ---------------------------------------------------------------

test_that("`%||%` replaces only NULL, not other falsy values", {
  null_coalesce <- SeuratExtend:::`%||%`
  expect_identical(null_coalesce(NULL, "fallback"), "fallback")
  expect_identical(null_coalesce("value", "fallback"), "value")
  # NULL is the *only* trigger - 0 / FALSE / "" must pass through unchanged
  expect_identical(null_coalesce(0, "fallback"), 0)
  expect_identical(null_coalesce(FALSE, "fallback"), FALSE)
  expect_identical(null_coalesce("", "fallback"), "")
})

test_that("is_empty() is length-based and NA-agnostic", {
  expect_true(SeuratExtend:::is_empty(character(0)))
  expect_true(SeuratExtend:::is_empty(NULL))
  expect_false(SeuratExtend:::is_empty(""))
  expect_false(SeuratExtend:::is_empty(NA))
  expect_false(SeuratExtend:::is_empty(list(1)))
})

# --- guards ---------------------------------------------------------------

test_that("check_spe() rejects NULL and passes anything else", {
  expect_error(SeuratExtend:::check_spe(NULL), "Species not defined")
  expect_true(SeuratExtend:::check_spe("human"))
  expect_true(SeuratExtend:::check_spe("mouse"))
})

test_that("check_pkg_version() compares against the installed version", {
  expect_null(SeuratExtend:::check_pkg_version("SeuratObject", "4.0.0"))
  expect_error(SeuratExtend:::check_pkg_version("SeuratObject", "99.0.0"),
               "is required")
})

test_that("message_only_once() fires once per session and then goes quiet", {
  opt <- "seuratextend_warning_internal_test_key"
  on.exit(options(structure(list(NULL), names = opt)), add = TRUE)
  options(structure(list(NULL), names = opt))  # clear any earlier state

  expect_message(
    SeuratExtend:::message_only_once("internal_test_key", "bench message"),
    "bench message"
  )
  expect_false(getOption(opt))
  # second call must be silent
  expect_silent(SeuratExtend:::message_only_once("internal_test_key", "bench message"))
})

# --- SeuratObject v4/v5 shim ---------------------------------------------

test_that(".get_assay_data_compat() returns the same slab as GetAssayData", {
  data(pbmc)
  counts <- SeuratExtend:::.get_assay_data_compat(pbmc, slot = "counts")
  expect_identical(dim(counts), dim(pbmc))
  ref <- if (utils::packageVersion("SeuratObject") >= "5.0.0") {
    SeuratObject::GetAssayData(pbmc, layer = "counts")
  } else {
    SeuratObject::GetAssayData(pbmc, slot = "counts")
  }
  expect_identical(counts, ref)
})

# --- MinMax_internal ------------------------------------------------------

test_that("MinMax_internal() clips only the bounds it is given", {
  mm <- SeuratExtend:::MinMax_internal
  x <- c(1, 5, 10)
  expect_identical(mm(x, min = 3), c(3, 5, 10))
  expect_identical(mm(x, max = 7), c(1, 5, 7))
  expect_identical(mm(x, min = 3, max = 7), c(3, 5, 7))
  expect_identical(mm(x), x)
  expect_identical(mm(x, min = NULL, max = NULL), x)
})

test_that("MinMax_internal() works on a matrix in place", {
  m <- matrix(c(1, 5, 10, 20), nrow = 2)
  expect_identical(SeuratExtend:::MinMax_internal(m, min = 3, max = 15),
                   matrix(c(3, 5, 10, 15), nrow = 2))
})

# --- Seu2Matr: input validation ------------------------------------------

test_that("Seu2Matr() returns the matrix plus the ident / split factors", {
  data(pbmc)
  out <- SeuratExtend:::Seu2Matr(pbmc, features = rownames(pbmc)[1:3])
  expect_identical(names(out), c("matr", "f", "f2"))
  expect_identical(dim(out$matr), c(ncol(pbmc), 3L))
  expect_s3_class(out$f, "factor")
  expect_length(out$f, ncol(pbmc))
  expect_identical(levels(out$f), levels(Seurat::Idents(pbmc)))
  expect_null(out$f2)          # no split.by -> NULL, not a factor
})

test_that("Seu2Matr() rejects a 'cells' value the object does not have", {
  data(pbmc)
  f <- function(...) SeuratExtend:::Seu2Matr(pbmc, features = rownames(pbmc)[1:3], ...)
  expect_error(f(cells = "NOT_A_CELL"), "'cells' not found in Seurat object")
  expect_error(f(cells = c(colnames(pbmc)[1], "NOT_A_CELL")), "not found in Seurat object")
})

test_that("Seu2Matr() rejects a logical 'cells' of the wrong length", {
  data(pbmc)
  expect_error(
    SeuratExtend:::Seu2Matr(pbmc, features = rownames(pbmc)[1:3],
                             cells = c(TRUE, FALSE)),
    "should be the same length as cells in Seurat object"
  )
})

test_that("Seu2Matr() rejects an unknown group.by / split.by", {
  data(pbmc)
  f <- function(...) SeuratExtend:::Seu2Matr(pbmc, features = rownames(pbmc)[1:3], ...)
  expect_error(f(group.by = "NOT_A_COLUMN"), "Cannot find 'NOT_A_COLUMN' in meta.data")
  expect_error(f(split.by = "NOT_A_COLUMN"), "Cannot find 'NOT_A_COLUMN' in meta.data")
})

test_that("Seu2Matr() builds the split factor from meta.data", {
  data(pbmc)
  out <- SeuratExtend:::Seu2Matr(pbmc, features = rownames(pbmc)[1:3],
                                 split.by = "condition")
  expect_s3_class(out$f2, "factor")
  expect_identical(names(out$f2), colnames(pbmc))
  expect_identical(levels(out$f2), levels(factor(pbmc$condition)))
})

test_that("Seu2Matr() warns when a feature name is also a meta.data column", {
  data(pbmc)
  # priority = 'expr' only warns when a requested feature is *both* a row of
  # the expression matrix and a column of meta.data, which never happens with
  # the shipped pbmc object - so plant the collision deliberately.
  collided <- pbmc
  collided$MALAT1 <- 1
  expect_warning(
    try(SeuratExtend:::Seu2Matr(collided, features = "MALAT1", priority = "expr"),
        silent = TRUE),
    "meta.data"
  )
})

# --- CheckIdent: contrast resolution -------------------------------------

test_that("CheckIdent() resolves ident.1 by name or by position", {
  f <- factor(c("a", "a", "b", "c"))
  ci <- SeuratExtend:::CheckIdent
  by_name <- ci(f, ident.1 = "b")
  expect_identical(by_name$name.1, "b")
  expect_identical(by_name$cell.1, c(FALSE, FALSE, TRUE, FALSE))
  by_pos <- ci(f, ident.1 = 2)
  expect_identical(by_pos$name.1, "b")
  expect_identical(by_pos$cell.1, by_name$cell.1)
})

test_that("CheckIdent() defaults ident.2 to everything else", {
  f <- factor(c("a", "a", "b", "c"))
  out <- SeuratExtend:::CheckIdent(f, ident.1 = "a")
  expect_identical(out$name.2, "non-a")
  expect_identical(out$cell.2, c(FALSE, FALSE, TRUE, TRUE))
  expect_false(any(out$cell.1 & out$cell.2))
})

test_that("CheckIdent() warns when both sides resolve to the same group", {
  f <- factor(c("a", "a", "b"))
  expect_warning(SeuratExtend:::CheckIdent(f, ident.1 = "a", ident.2 = "a"),
                 "are the same")
})

test_that("CheckIdent() rejects an unknown ident and lists the alternatives", {
  f <- factor(c("a", "a", "b"))
  expect_error(SeuratExtend:::CheckIdent(f, ident.1 = "zzz"), "not found")
  expect_error(SeuratExtend:::CheckIdent(f, ident.1 = "a", ident.2 = "zzz"),
               "not found")
})
