# Test file for R/Scenic.R
#
# Scenic.R holds a single export, ImportPyscenicLoom(), which reads a
# pySCENIC loom file through loomR and folds the result into a Seurat object
# (a "TF" assay plus seu@misc$SCENIC).
#
# loomR is not a declared dependency of this package - it is only reachable
# through the internal import() helper - so on a stock install the function
# dies part-way through. That path is what the behavioural test below pins.
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
  # full path needs loomR plus a real loom fixture.
  body_txt <- deparse(body(SeuratExtend::ImportPyscenicLoom))
  expect_true(any(grepl("close_all", body_txt, fixed = TRUE)))
  # row attributes carry the gene names used to build the regulon list
  expect_true(any(grepl("row.attrs$Gene", body_txt, fixed = TRUE)))
})

test_that("ImportPyscenicLoom() reports the missing backend before it fails", {
  skip_if(requireNamespace("loomR", quietly = TRUE),
          "loomR is installed, so the missing-dependency path is not reachable")
  # KNOWN-ISSUE (Scenic / L3-6 family): import("loomR") only *messages* and
  # returns FALSE, and the return value is discarded, so execution continues
  # into connect() and dies with a bare "could not find function". The user
  # has to reason backwards from an unrelated-looking error to the real cause.
  # Source: R/Scenic.R:27-30, R/import.R:87-134.
  expect_message(
    expect_error(
      quietly(ImportPyscenicLoom("definitely-not-a-loom-file.loom")),
      "could not find function"
    ),
    "Required package"
  )
})
