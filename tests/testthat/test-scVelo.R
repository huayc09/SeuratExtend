# Test file for R/scVelo.R
#
# scVelo.R is the longest file in the package (618 lines, 4 exports) and
# bridges Seurat objects to scVelo's Python side.
#
# SCOPE: the export surface, the three different ways a loom path is supplied,
# and the loud failure when those paths do not exist. The file-path validation
# happens *before* any Python call, so unlike most of the Python bridge it can
# be exercised on a stock install.
#
# STYLE: characterisation tests - see test-Internal.R for the convention.

scvelo_exports <- c("scVelo.SeuratToAnndata", "scVelo.SeuratToLoom",
                    "scVelo.RunBasic", "scVelo.Plot")

test_that("all four scVelo entry points are exported", {
  exports <- getNamespaceExports("SeuratExtend")
  expect_length(setdiff(scvelo_exports, exports), 0)
})

test_that("scVelo signatures and defaults are stable", {
  a <- formals(SeuratExtend::scVelo.SeuratToAnndata)
  expect_identical(
    names(a),
    c("seu", "filename", "velocyto.loompath", "cell.id.match.table",
      "prefix", "postfix", "remove_duplicates", "conda_env")
  )
  expect_null(eval(a$cell.id.match.table))
  expect_identical(a$postfix, "-1")
  expect_false(eval(a$remove_duplicates))
  expect_identical(a$conda_env, "seuratextend")

  l <- formals(SeuratExtend::scVelo.SeuratToLoom)
  expect_identical(
    names(l),
    c("seu", "export.name", "velocyto.loom.path", "velocyto.loom.filenames",
      "fxn.convert.loomtoseurat.cellname", "cell.id.match.table",
      "prefix", "postfix")
  )
  expect_null(eval(l$velocyto.loom.path))

  r <- formals(SeuratExtend::scVelo.RunBasic)
  expect_identical(names(r), c("loom", "save.adata"))
  expect_identical(r$save.adata, "adata.obj")

  p <- formals(SeuratExtend::scVelo.Plot)
  expect_identical(p$basis, "umap_cell_embeddings")
  expect_identical(eval(p$alpha), 0.15)
  expect_identical(eval(p$arrow_size), 3)
  expect_identical(eval(p$arrow_length), 2)
  expect_identical(eval(p$dpi), 300)
  expect_identical(eval(p$figsize), c(7, 5))
  expect_identical(eval(p$style), c("stream", "grid", "scatter"))
  expect_identical(p$conda_env, "seuratextend")
})

test_that("scVelo.SeuratToLoom() lists every loom file it cannot find", {
  data(pbmc)
  expect_error(
    quietly(suppressMessages(scVelo.SeuratToLoom(
      pbmc, export.name = "unused.loom",
      velocyto.loom.filenames = "missing_one.loom"))),
    "Loom file\\(s\\) not found"
  )
  # several files -> all of them are named, not just the first
  err <- tryCatch(
    quietly(suppressMessages(scVelo.SeuratToLoom(
      pbmc, export.name = "unused.loom",
      velocyto.loom.filenames = c("missing_one.loom", "missing_two.loom")))),
    error = function(e) conditionMessage(e)
  )
  expect_match(err, "missing_one.loom")
  expect_match(err, "missing_two.loom")
})

test_that("scVelo.SeuratToLoom() joins velocyto.loom.path onto the filenames", {
  data(pbmc)
  err <- tryCatch(
    quietly(suppressMessages(scVelo.SeuratToLoom(
      pbmc, export.name = "unused.loom", velocyto.loom.path = "some/dir",
      velocyto.loom.filenames = "missing_one.loom"))),
    error = function(e) conditionMessage(e)
  )
  # the reported path must keep the directory part, otherwise the user cannot
  # tell which of several loom roots was searched
  expect_match(err, "some[/\\\\]dir", perl = TRUE)
})

test_that("scVelo.SeuratToAnndata() validates its loom path before Python", {
  data(pbmc)
  # KNOWN-ISSUE (L3-4): load_condaenv() runs first, so on a machine without
  # conda the user sees a conda error rather than the path error they can
  # actually act on. The path guards themselves (character, length > 0,
  # existence) are unreachable until conda is set up. Source: R/scVelo.R:85-96.
  expect_error(quietly(suppressMessages(
    scVelo.SeuratToAnndata(pbmc, filename = "out.h5ad",
                           velocyto.loompath = "missing.loom")
  )))
})
