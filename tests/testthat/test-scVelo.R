# Test file for R/scVelo.R
#
# scVelo.R is the longest file in the package (511 lines, 4 exports) and
# bridges Seurat objects to scVelo's Python side.
#
# SCOPE: the defunct contract of scVelo.SeuratToLoom (it must fail loudly and
# name its replacement), the export surface, the signatures of all four entry
# points, and scVelo.SeuratToAnndata's pre-Python failure behavior (conda-env
# check precedes the loom-path guards — see the L3-4 pin below).
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

test_that("scVelo.SeuratToLoom() is defunct and names its replacement", {
  # The function body is now .Defunct(): the loom backend dependency used to
  # build an integrated loom was removed, so every call must fail loudly and
  # point at the supported alternative instead of reaching for Python.
  data(pbmc)
  expect_error(
    quietly(scVelo.SeuratToLoom(
      pbmc, export.name = "unused.loom",
      velocyto.loom.filenames = "missing_one.loom")),
    "defunct"
  )
  err <- tryCatch(
    quietly(scVelo.SeuratToLoom(
      pbmc, export.name = "unused.loom",
      velocyto.loom.filenames = "missing_one.loom")),
    error = function(e) conditionMessage(e)
  )
  expect_match(err, "scVelo\\.SeuratToAnndata")
  expect_match(err, "loom backend")
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
