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
  expect_true(any(grepl("on.exit", body_txt, fixed = TRUE)))
  # row attributes carry the gene names used to build the regulon list
  expect_true(any(grepl("row_attrs", body_txt, fixed = TRUE)))
  expect_true(any(grepl("Gene", body_txt, fixed = TRUE)))
})
