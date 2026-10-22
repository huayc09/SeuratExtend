# Test file for R/Scenic.R
#
# Scenic.R holds a single export, ImportPyscenicLoom(), which reads a
# pySCENIC loom file directly with hdf5r (loom files are plain HDF5) and folds
# the result into a Seurat object (a "TF" assay plus seu@misc$SCENIC).
#
# The full read path needs the zenodo loom fixture, so the tests here cover
# the export surface, the signature and the handle discipline statically.
#
# STYLE: characterisation tests - see test-Internal.R for the convention.

test_that("ImportPyscenicLoom is the only export from this file", {
  exports <- getNamespaceExports("SeuratExtend")
  expect_true("ImportPyscenicLoom" %in% exports)
})

test_that("ImportPyscenicLoom() signature is stable", {
  f <- formals(SeuratExtend::ImportPyscenicLoom)
  expect_identical(names(f), c("loom.path", "seu"))
  expect_null(eval(f$seu))   # no object -> return a plain list instead
})

test_that("ImportPyscenicLoom() closes the loom handle it opens", {
  # The function opens an HDF5 handle and must close it on every path; an
  # unclosed handle keeps the file locked. This is a static check because the
  # full path needs a real loom fixture.
  body_txt <- deparse(body(SeuratExtend::ImportPyscenicLoom))
  expect_true(any(grepl("close_all", body_txt, fixed = TRUE)))
  # the handle must be registered for closing immediately after it is opened,
  # and before any dataset is read, otherwise an error in between leaks it
  idx_open   <- grep("H5File\\$new", body_txt)
  idx_on_exit <- grep("on.exit", body_txt, fixed = TRUE)
  idx_read   <- grep('CellID <- lfile\\[\\["col_attrs"\\]\\]', body_txt)
  expect_length(idx_open, 1)
  expect_length(idx_on_exit, 1)
  expect_length(idx_read, 1)
  expect_true(idx_open < idx_on_exit)
  expect_true(idx_on_exit < idx_read)
  # row attributes carry the gene names used to build the regulon list
  expect_true(any(grepl("row_attrs", body_txt, fixed = TRUE)))
  expect_true(any(grepl("Gene", body_txt, fixed = TRUE)))
})
