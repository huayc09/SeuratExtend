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

test_that("RunBasicSeurat reaches harmony through an undeclared import()", {
  # KNOWN-ISSUE (RBS-1, P1) -- R/RunBasicSeurat.R:98 calls import("harmony")
  # and discards the return value, then :111 calls RunHarmony() directly.
  # harmony is in NO field of DESCRIPTION, and import() needs an interactive
  # answer when the package is missing -- so in a script / knitr / CI session
  # the "auto install" branch cannot run and the call dies with
  # `could not find function "RunHarmony"`, naming nothing the user can fix.
  # (This is the same family as GeneSetAnalysis's DelayedArray hole.)
  root <- file.path("..", "..")
  skip_if_not(dir.exists(file.path(root, "R")), "source tree not reachable")
  desc <- paste(readLines(file.path(root, "DESCRIPTION"), warn = FALSE),
                collapse = " ")
  expect_false(grepl("harmony", desc, fixed = TRUE))
  src <- paste(readLines(file.path(root, "R", "RunBasicSeurat.R"),
                         warn = FALSE), collapse = "\n")
  expect_true(grepl('import("harmony")', src, fixed = TRUE))
  expect_true(grepl("RunHarmony(", src, fixed = TRUE))
})
