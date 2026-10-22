# Test file for R/RunBasicSeurat.R
# Covers: RunBasicSeurat (standard QC -> clustering -> UMAP pipeline)
# Data: built-in pbmc
#
# RunBasicSeurat re-runs the full pipeline on a (possibly raw) Seurat object.
# On an already-processed object it must still complete and return a Seurat.

test_that("RunBasicSeurat returns a Seurat object with clusters", {
  data(pbmc)
  out <- RunBasicSeurat(pbmc)
  expect_s4_class(out, "Seurat")
  # pipeline produced identities and reductions
  expect_gt(nlevels(Idents(out)), 0)
  expect_true("umap" %in% Reductions(out))
})

test_that("RunBasicSeurat accepts cluster resolution parameter", {
  data(pbmc)
  out <- RunBasicSeurat(pbmc, resolution = 0.8)
  expect_s4_class(out, "Seurat")
})

# --- KNOWN-ISSUE pin (2026-09-16 recheck against 问题清单0915) ------------

test_that("RunBasicSeurat guards its harmony dependency explicitly", {
  # FIXED (RBS-1, P1): R/RunBasicSeurat.R no longer drops import()'s return
  # value. The call site is now import_required("harmony"), which keeps
  # import()'s interactive install prompt and then fails fast with an install
  # hint when the package is absent (script / knitr / CI sessions), instead of
  # dying later at RunHarmony() with `could not find function "RunHarmony"`.
  # harmony is declared in DESCRIPTION Suggests by this work order.
  root <- file.path("..", "..")
  skip_if_not(dir.exists(file.path(root, "R")), "source tree not reachable")
  desc <- paste(readLines(system.file("DESCRIPTION", package = "SeuratExtend"),
                          warn = FALSE), collapse = " ")
  expect_true(grepl("harmony", desc, fixed = TRUE))
  src <- paste(readLines(file.path(root, "R", "RunBasicSeurat.R"),
                         warn = FALSE), collapse = "\n")
  expect_true(grepl('import_required("harmony")', src, fixed = TRUE))
  expect_false(grepl('import("harmony")', src, fixed = TRUE))
  expect_true(grepl("RunHarmony(", src, fixed = TRUE))
})
