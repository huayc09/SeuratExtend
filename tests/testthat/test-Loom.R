# Test file for R/Loom.R
# Covers: Seu2Loom (Seurat object -> loom / HDF5 exporter; the only top-level
#         object in this file)
# Regression issues: L3-3 checklist (9 items: 3x P2 / 5x P3 / 1 pending)
# Data: built-in pbmc (500 cells)
#
# NOTE every test writes its loom file into tempdir() and unlinks it on exit --
# nothing is written into the package working directory. hdf5r is an optional
# runtime dependency of Seu2Loom (attached by a bare library() inside the
# function), so the whole file skips when it is absent.

loom_tmp <- function(name) file.path(tempdir(), name)

h5names <- function(file, group = NULL) {
  h5 <- hdf5r::H5File$new(file, mode = "r")
  on.exit(h5$close_all(), add = TRUE)
  if (is.null(group)) names(h5) else names(h5[[group]])
}

# --- happy path -----------------------------------------------------------

test_that("Seu2Loom() writes the documented loom layout", {
  # Source: R/Loom.R:49-108.
  skip_if_not_installed("hdf5r")
  data(pbmc)
  f <- loom_tmp("test-seu2loom-basic.loom")
  unlink(f)
  on.exit(unlink(f), add = TRUE)

  expect_no_error(Seu2Loom(pbmc, f))
  expect_true(file.exists(f))
  expect_gt(file.size(f), 0)

  expect_setequal(h5names(f),
                  c("attrs", "col_attrs", "col_graphs", "layers",
                    "matrix", "row_attrs", "row_graphs"))
  # one dataset per meta.data column plus CellID
  expect_true("CellID" %in% h5names(f, "col_attrs"))
  expect_true(all(colnames(pbmc@meta.data) %in% h5names(f, "col_attrs")))
  expect_true("Gene" %in% h5names(f, "row_attrs"))
})

test_that("Seu2Loom() add.normdata = TRUE moves counts into /layers", {
  # Source: R/Loom.R:66-73.
  skip_if_not_installed("hdf5r")
  data(pbmc)
  f <- loom_tmp("test-seu2loom-norm.loom")
  unlink(f)
  on.exit(unlink(f), add = TRUE)

  expect_no_error(Seu2Loom(pbmc, f, add.normdata = TRUE))
  expect_true("counts" %in% h5names(f, "layers"))
})

test_that("Seu2Loom() add.metadata = FALSE writes only CellID", {
  # Source: R/Loom.R:84-98 (the whole block is gated on add.metadata).
  skip_if_not_installed("hdf5r")
  data(pbmc)
  f <- loom_tmp("test-seu2loom-nometa.loom")
  unlink(f)
  on.exit(unlink(f), add = TRUE)

  expect_no_error(Seu2Loom(pbmc, f, add.metadata = FALSE))
  expect_identical(h5names(f, "col_attrs"), "CellID")
})

test_that("Seu2Loom() writes LOOM_SPEC_VERSION as a scalar file attribute", {
  # CellWeave A1 / audit D1-D2: loom 3.0 puts the version at the FILE level,
  # and loompy 3.0.6 get_loom_spec_version reads /attrs/LOOM_SPEC_VERSION
  # FIRST. A plain h5attr()-style write lands as a 1x1x1 SIMPLE dataspace,
  # which hdf5r folds back into a length-1 character but loompy sees as an
  # ndarray, so compare_loom_spec_version()/ds.ca writes still raise
  # AttributeError. Assert the real dataspace: 0 == H5S_SCALAR.
  # Source: R/Loom.R:52-55 (create_group + create_attr(space = H5S$new("scalar"))).
  skip_if_not_installed("hdf5r")
  data(pbmc)
  f <- loom_tmp("test-seu2loom-specversion.loom")
  unlink(f)

  expect_no_error(Seu2Loom(pbmc, f))
  expect_true("attrs" %in% h5names(f))
  expect_false("LOOM_SPEC_VERSION" %in% h5names(f, "attrs"))

  h5 <- hdf5r::H5File$new(f, mode = "r")
  spec <- h5$attr_open_by_name("LOOM_SPEC_VERSION", ".")
  on.exit(spec$close(), add = TRUE)
  # Single handler, close before unlink (Windows-friendly; audit D5).
  on.exit({ h5$close_all(); unlink(f) }, add = TRUE)

  # hdf5r returns a factor_ext here, so compare by its HDF5 name.
  expect_identical(as.character(spec$get_space()$get_simple_extent_type()),
                   "H5S_SCALAR")

  attrs <- hdf5r::h5attributes(h5)
  expect_true("LOOM_SPEC_VERSION" %in% names(attrs))
  expect_type(attrs$LOOM_SPEC_VERSION, "character")
  expect_length(attrs$LOOM_SPEC_VERSION, 1)
  expect_identical(unname(attrs$LOOM_SPEC_VERSION), "3.0.0")
})

test_that("Seu2Loom() output satisfies loompy's version/attribute contract", {
  # Consumer-side regression for CellWeave A1 / audit D3: connect(validate =
  # TRUE) alone passes even on the non-scalar defect (LoomValidator only
  # compares), so assert the semantics pyscenic actually exercises:
  # get_loom_spec_version() returns a str, compare_loom_spec_version() == 0,
  # and a mode = "r+" session can append a col_attr (the AUCell path).
  # No conda env / no loompy -> SKIP; those skips go to the settlement
  # uncovered list.
  skip_if_not_installed("hdf5r")
  skip_if_not_installed("reticulate")
  if (!reticulate::py_available(initialize = FALSE)) {
    reticulate::use_condaenv("seuratextend", required = FALSE)
  }
  if (!reticulate::py_available(initialize = TRUE)) skip("no usable Python")
  skip_if_not(reticulate::py_module_available("loompy"))
  data(pbmc)
  f <- loom_tmp("test-seu2loom-loompy.loom")
  unlink(f)
  expect_no_error(Seu2Loom(pbmc, f))

  h5py <- reticulate::import("h5py")
  loompy <- reticulate::import("loompy")
  np <- reticulate::import("numpy")

  fh <- h5py$File(f, "r")
  on.exit(fh$close(), add = TRUE)
  ver <- loompy$utils$get_loom_spec_version(fh)
  expect_type(ver, "character")
  expect_identical(ver, "3.0.0")
  expect_equal(loompy$utils$compare_loom_spec_version(fh, "3.0.0"), 0L)
  # h5py holds a read lock: release it before opening the same file r+.
  fh$close()

  ds <- loompy$connect(f, mode = "r+")
  on.exit(ds$close(), add = TRUE)
  # Close the Python connections before unlinking (audit D5).
  on.exit(unlink(f), add = TRUE)
  # loom stores genes x cells; reticulate hands the tuple back as a list.
  expect_equal(as.integer(unlist(ds$shape)), c(nrow(pbmc), ncol(pbmc)))
  expect_no_error({
    ds$ca[["RegulonsAUC"]] <- np$zeros(c(ncol(pbmc), 3L))
  })
})

# --- KNOWN-ISSUE: the layer guard -----------------------------------------

test_that("Seu2Loom() guards a NAMED layer whose dims are wrong", {
  # Contrast case for the two issues below: for a named matrix the guard does
  # fire. Source: R/Loom.R:35-43.
  skip_if_not_installed("hdf5r")
  data(pbmc)
  f <- loom_tmp("test-seu2loom-badlayer.loom")
  unlink(f)
  on.exit(unlink(f), add = TRUE)

  bad <- matrix(1:6, nrow = 2)
  expect_error(Seu2Loom(pbmc, f, layers = list(bad = bad)),
               "do not match the transposed dimensions")
})

test_that("Seu2Loom() lets a NON-MATRIX layer through the dimension guard", {
  # KNOWN-ISSUE (LOM-1, P2): the guard is `all(seu_dim == rev(dim(layer)))`.
  # For anything without a dim() -- a plain vector, a list, a sparse Matrix --
  # dim() is NULL, the comparison is logical(0), and all() of an empty vector
  # is TRUE, so the malformed layer is written as-is with no error.
  # Source: R/Loom.R:39 (:35-43 is the whole guard); write happens at :103.
  skip_if_not_installed("hdf5r")
  data(pbmc)
  f <- loom_tmp("test-seu2loom-veclayer.loom")
  unlink(f)
  on.exit(unlink(f), add = TRUE)

  expect_no_error(Seu2Loom(pbmc, f, layers = list(weird = 1:10)))
  expect_true("weird" %in% h5names(f, "layers"))
})

test_that("Seu2Loom() silently drops unnamed and empty layer lists", {
  # KNOWN-ISSUE (LOM-2, P2): the write loop iterates over names(layers), so a
  # positional (unnamed) list is dropped without a word -- and so is list().
  # There is no warning even though the caller asked for extra layers.
  # Source: R/Loom.R:35 (check) / :101-102 (write loop).
  skip_if_not_installed("hdf5r")
  data(pbmc)
  f <- loom_tmp("test-seu2loom-unnamed.loom")
  unlink(f)
  on.exit(unlink(f), add = TRUE)

  expect_silent(Seu2Loom(pbmc, f, layers = list(matrix(1:6, nrow = 2))))
  expect_length(h5names(f, "layers"), 0)

  f2 <- loom_tmp("test-seu2loom-emptylayers.loom")
  unlink(f2)
  on.exit(unlink(f2), add = TRUE)
  expect_silent(Seu2Loom(pbmc, f2, layers = list()))
  expect_length(h5names(f2, "layers"), 0)
})

# --- KNOWN-ISSUE: booleans, extension, metadata typing --------------------

test_that("Seu2Loom() overwrite only accepts the literal TRUE", {
  # KNOWN-ISSUE (LOM-4, P3): all three flags go through isTRUE(), so a truthy
  # non-TRUE value (1, "yes", NA) behaves exactly like FALSE.
  # Source: R/Loom.R:23 (overwrite), :66 (add.normdata), :84 (add.metadata).
  skip_if_not_installed("hdf5r")
  data(pbmc)
  f <- loom_tmp("test-seu2loom-ow.loom")
  unlink(f)
  on.exit(unlink(f), add = TRUE)

  expect_no_error(Seu2Loom(pbmc, f))
  expect_error(Seu2Loom(pbmc, f, overwrite = FALSE), "already exists")
  expect_error(Seu2Loom(pbmc, f, overwrite = 1), "already exists")
  expect_error(Seu2Loom(pbmc, f, overwrite = "yes"), "already exists")
  expect_warning(Seu2Loom(pbmc, f, overwrite = TRUE),
                 "Overwriting previous file")
})

test_that("Seu2Loom() appends .loom case-sensitively", {
  # KNOWN-ISSUE (LOM-8, P3): the suffix test is grepl("^loom$", file_ext(x)),
  # which is case sensitive and never replaces an existing extension -- so
  # "x.LOOM" becomes "x.LOOM.loom" and "x.txt" becomes "x.txt.loom".
  # Source: R/Loom.R:17-19.
  skip_if_not_installed("hdf5r")
  data(pbmc)

  upper <- loom_tmp("test-seu2loom-upper.LOOM")
  made_upper <- paste0(upper, ".loom")
  unlink(c(upper, made_upper))
  on.exit(unlink(c(upper, made_upper)), add = TRUE)
  expect_no_error(Seu2Loom(pbmc, upper))
  expect_true(file.exists(made_upper))        # appended, original name unused
  expect_false(file.exists(upper))

  lower <- loom_tmp("test-seu2loom-lower.loom")
  unlink(c(lower, paste0(lower, ".loom")))
  on.exit(unlink(c(lower, paste0(lower, ".loom"))), add = TRUE)
  expect_no_error(Seu2Loom(pbmc, lower))
  expect_true(file.exists(lower))
})

test_that("Seu2Loom() degrades Date meta.data to a bare numeric column", {
  # KNOWN-ISSUE (LOM-5, P3): only factor columns are converted (to character);
  # everything else is handed straight to create_dataset(), and hdf5r stores a
  # Date as its underlying day number -- the type is lost with no warning.
  # Source: R/Loom.R:88-89 (the only conversion) / :92-96 (everything else).
  skip_if_not_installed("hdf5r")
  data(pbmc)
  p2 <- pbmc
  p2$datestamp <- as.Date("2020-01-01") + seq_len(ncol(p2))
  f <- loom_tmp("test-seu2loom-date.loom")
  unlink(f)
  on.exit(unlink(f), add = TRUE)

  expect_no_error(Seu2Loom(p2, f))
  h5 <- hdf5r::H5File$new(f, mode = "r")
  on.exit(h5$close_all(), add = TRUE)
  stored <- h5[["col_attrs/datestamp"]]$read()
  expect_true(is.numeric(stored))
  expect_false(inherits(stored, "Date"))
})

test_that("Seu2Loom() returns the already-closed H5File handle", {
  # KNOWN-ISSUE (LOM-6, P3): the last statement is file$close_all(), so the
  # value that escapes is an R6 handle whose file has been closed -- not NULL
  # and not the "loom file" the adata.Load docs promise.
  # Source: R/Loom.R:108 (last line of the function body).
  skip_if_not_installed("hdf5r")
  data(pbmc)
  f <- loom_tmp("test-seu2loom-return.loom")
  unlink(f)
  on.exit(unlink(f), add = TRUE)

  v <- Seu2Loom(pbmc, f)
  expect_false(is.null(v))
  expect_true(inherits(v, "H5File"))
})

# --- KNOWN-ISSUE: structural facts pinned at the source level -------------

test_that("Seu2Loom() has no on.exit/tryCatch around the H5 handle", {
  # KNOWN-ISSUE (LOM-3, P2): the file is opened at :49 and closed only at the
  # very last line (:108). A failure anywhere in between leaves a half-written
  # .loom on disk and, on Windows, a locked handle.
  # Source: R/Loom.R:49 (open) / :108 (the single close).
  body_txt <- paste(deparse(body(SeuratExtend::Seu2Loom)), collapse = "\n")
  expect_false(grepl("on.exit", body_txt, fixed = TRUE))
  expect_false(grepl("tryCatch", body_txt, fixed = TRUE))
  expect_true(grepl("close_all", body_txt, fixed = TRUE))
})

test_that("Seu2Loom() attaches three packages with a bare library()", {
  # KNOWN-ISSUE (LOM-7, P3): the whole package's bare-library family (X5);
  # these three attach onto the caller's search path at position 2.
  # Source: R/Loom.R:12 (hdf5r) / :13 (Seurat) / :14 (tools).
  body_txt <- paste(deparse(body(SeuratExtend::Seu2Loom)), collapse = "\n")
  for (pkg in c("hdf5r", "Seurat", "tools")) {
    expect_true(grepl(sprintf("library(%s)", pkg), body_txt, fixed = TRUE))
  }
})
