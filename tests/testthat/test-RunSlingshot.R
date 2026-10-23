# Test file for R/RunSlingshot.R
# Covers: RunSlingshot (Seurat <-> slingshot trajectory wrapper; the only
#         top-level object in this file)
# Regression issues: L3-9 checklist (14 items: 3x P2 / 11x P3)
# Data: built-in pbmc (500 cells; reductions pca + umap, cluster has 9 levels)
#
# NOTE the reduction NAME that works is whatever as.SingleCellExperiment()
# stored, i.e. upper case ("UMAP"), NOT Seurat's lower case -- see RSS-14.
# All calls are wrapped so the slingshot attach banner and the function's own
# "Time usage of slingshot" message do not leak into the test log.

skip_if_no_slingshot <- function() {
  skip_if_not_installed("slingshot")
  skip_if_not_installed("SingleCellExperiment")
}

rs_quiet <- function(expr) {
  suppressMessages(suppressPackageStartupMessages(expr))
}

# --- main path ------------------------------------------------------------

test_that("RunSlingshot() runs on a reduction named the way SCE stores it", {
  # Source: R/RunSlingshot.R:43 (as.SingleCellExperiment) / :45 (slingshot()).
  skip_if_no_slingshot()
  data(pbmc)
  out <- rs_quiet(RunSlingshot(pbmc, group.by = "cluster", reducedDim = "UMAP"))

  expect_true(inherits(out, "Seurat"))
  expect_true("slingshot" %in% names(out@misc))
  expect_true("UMAP" %in% names(out@misc$slingshot))
  expect_setequal(names(out@misc$slingshot$UMAP),
                  c("SlingshotDataSet", "f", "SlingPseudotime"))
  expect_true(is.factor(out@misc$slingshot$UMAP$f))
  expect_length(out@misc$slingshot$UMAP$f, ncol(pbmc))
  expect_true(all(grepl("^slingPseudotime_\\d+$",
                        colnames(out@misc$slingshot$UMAP$SlingPseudotime))))
})

test_that("RunSlingshot() writes a 'f' key that nothing in the package reads", {
  # KNOWN-ISSUE (RSS-6, P3): :51 stores the group labels under "f"; grep over
  # R/ finds no consumer of that key.
  # Source: R/RunSlingshot.R:51.
  skip_if_no_slingshot()
  data(pbmc)
  out <- rs_quiet(RunSlingshot(pbmc, group.by = "cluster", reducedDim = "UMAP"))
  expect_true("f" %in% names(out@misc$slingshot$UMAP))
})

# --- KNOWN-ISSUE: the reduction name / default ----------------------------

test_that("RunSlingshot() rejects Seurat's lower-case reduction names", {
  # KNOWN-ISSUE (RSS-14, P3): as.SingleCellExperiment upper-cases the entries
  # of reducedDims, so "pca"/"umap" (the names Seurat itself uses) are not
  # found. The error comes from slingshot and does not explain the casing.
  # Source: R/RunSlingshot.R:32 (default "PCA") / :45 (pass-through).
  skip_if_no_slingshot()
  data(pbmc)
  expect_error(rs_quiet(RunSlingshot(pbmc, group.by = "cluster", reducedDim = "umap")),
               "not found in reducedDims")
  expect_error(rs_quiet(RunSlingshot(pbmc, group.by = "cluster", reducedDim = "pca")),
               "not found in reducedDims")
})

test_that("RunSlingshot() default reducedDim = 'PCA' dies on pbmc", {
  # KNOWN-ISSUE (RSS-3, P2): with the documented default the call reaches
  # slingshot and fails numerically -- pbmc's 50-PC matrix is singular -- and
  # the only knob the user has is to pass a different whole reduction.
  # Source: R/RunSlingshot.R:32 / :45.
  skip_if_no_slingshot()
  data(pbmc)
  expect_error(rs_quiet(RunSlingshot(pbmc, group.by = "cluster")),
               "computationally singular")
  expect_error(rs_quiet(RunSlingshot(pbmc, group.by = "cluster", reducedDim = "PCA")),
               "computationally singular")
})

test_that("RunSlingshot() cannot accept a matrix reducedDim", {
  # KNOWN-ISSUE (RSS-2, P2): slingshot's own API accepts a matrix here, but
  # :49 uses reducedDim as a list key, so a matrix dies with an indexing error
  # that says nothing about reducedDim.
  # Source: R/RunSlingshot.R:45 (pass-through) / :49 (used as a key).
  skip_if_no_slingshot()
  data(pbmc)
  expect_error(rs_quiet(RunSlingshot(pbmc, group.by = "cluster",
                                     reducedDim = Embeddings(pbmc, "umap"))),
               "no such index")
})

# --- KNOWN-ISSUE: meta.data net loss --------------------------------------

test_that("RunSlingshot() deletes slingPseudotime_* columns and writes none back", {
  # KNOWN-ISSUE (RSS-1, P2): :40-41 strips every slingPseudotime_<n> column
  # from meta.data up front, and this function only ever returns the SCE
  # result inside @misc -- it never repopulates those columns. A user who had
  # them (e.g. from an earlier run) silently loses the data.
  # Source: R/RunSlingshot.R:40-41 (delete) / :54 (return, no write-back).
  skip_if_no_slingshot()
  data(pbmc)
  p2 <- pbmc
  p2$slingPseudotime_1 <- seq_len(ncol(p2))
  p2$slingPseudotime_2 <- seq_len(ncol(p2))
  p2$unrelated <- "keep"
  expect_true(all(c("slingPseudotime_1", "slingPseudotime_2") %in%
                    colnames(p2@meta.data)))

  out <- rs_quiet(RunSlingshot(p2, group.by = "cluster", reducedDim = "UMAP"))

  expect_false("slingPseudotime_1" %in% colnames(out@meta.data))
  expect_false("slingPseudotime_2" %in% colnames(out@meta.data))
  expect_true("unrelated" %in% colnames(out@meta.data))
  # nothing of the sort is written back, even though this run produced 3
  expect_length(grep("slingPseudotime", colnames(out@meta.data), value = TRUE), 0)
  expect_length(colnames(out@misc$slingshot$UMAP$SlingPseudotime), 3)
})

test_that("RunSlingshot() leaves the caller's object untouched", {
  # The function is documented as "Modifies the Seurat object ...", but R copy
  # semantics mean the input keeps its columns and the caller must use the
  # return value. Hence the wording drift pinned as RSS-11.
  # Source: R/RunSlingshot.R:8 (doc) / :54 (return(Seu)).
  skip_if_no_slingshot()
  data(pbmc)
  p2 <- pbmc
  p2$slingPseudotime_1 <- seq_len(ncol(p2))
  invisible(rs_quiet(RunSlingshot(p2, group.by = "cluster", reducedDim = "UMAP")))
  expect_true("slingPseudotime_1" %in% colnames(p2@meta.data))
})

# --- KNOWN-ISSUE: parameter validation ------------------------------------

test_that("RunSlingshot() silently ignores an unknown start.clus", {
  # KNOWN-ISSUE (RSS-4, P3): start.clus is forwarded with no check against the
  # cluster levels, and slingshot does not complain either -- the run simply
  # proceeds as if no start had been given.
  # Source: R/RunSlingshot.R:34 (formal, unchecked) / :45 (pass-through).
  skip_if_no_slingshot()
  data(pbmc)
  out <- rs_quiet(RunSlingshot(pbmc, group.by = "cluster", reducedDim = "UMAP",
                               start.clus = "NoSuchCluster"))
  expect_true(inherits(out, "Seurat"))
})

test_that("RunSlingshot() multi-column group.by error does not mention 'one column'", {
  # KNOWN-ISSUE (RSS-12, P3): group.by is documented as "a metadata column"
  # but never checked; two names reach slingshot and come back as a length
  # mismatch that never mentions group.by.
  # Source: R/RunSlingshot.R:33 (formal, unchecked) / :45 (pass-through).
  skip_if_no_slingshot()
  data(pbmc)
  expect_error(rs_quiet(RunSlingshot(pbmc, group.by = c("cluster", "condition"),
                                     reducedDim = "UMAP")),
               "must equal length")
})

# --- source-level pins (RSS-7 fixed by the deps-guard work order, RSS-8 open)

test_that("RunSlingshot() guards its slingshot dependency explicitly", {
  # FIXED (RSS-7, P3): import()'s FALSE used to be discarded here, so a script
  # / knitr / CI session continued past a refused or skipped install and died
  # later inside slingshot()/SlingshotDataSet() with `could not find function`.
  # The call site is now import_required("slingshot"), which keeps import()'s
  # interactive install prompt and then fails fast with an install hint when
  # the package is absent.
  # Source: R/RunSlingshot.R:37 (import_required).
  body_txt <- paste(deparse(body(SeuratExtend::RunSlingshot)), collapse = "\n")
  expect_true(grepl('import_required("slingshot")', body_txt, fixed = TRUE))
  expect_false(grepl('import("slingshot")', body_txt, fixed = TRUE))
})

test_that("RunSlingshot() uses bare library() to attach Seurat and dplyr", {
  # KNOWN-ISSUE (RSS-8, P3): the two bare library() calls attach Seurat and
  # dplyr onto the caller's search path (Seurat is a real functional dependency
  # here, dplyr provides %>%). Deliberately unchanged by the deps-guard work
  # order, so this part of the old combined pin still asserts current behaviour.
  # Source: R/RunSlingshot.R:38-39 (library).
  body_txt <- paste(deparse(body(SeuratExtend::RunSlingshot)), collapse = "\n")
  expect_true(grepl("library(Seurat)", body_txt, fixed = TRUE))
  expect_true(grepl("library(dplyr)", body_txt, fixed = TRUE))
})

test_that("RunSlingshot() shadows base::message locally", {
  # Noted as harmless in the checklist (exclusion e1): :47 assigns a character
  # vector to a local called `message` and :48 then calls it. R resolves the
  # call position to base::message, so the message still prints -- bad smell
  # only, and it is the same family as the return/quit shadowing elsewhere.
  # Source: R/RunSlingshot.R:47-48.
  body_txt <- paste(deparse(body(SeuratExtend::RunSlingshot)), collapse = "\n")
  expect_true(grepl("message <- paste0", body_txt, fixed = TRUE))
  expect_true(grepl("message(message)", body_txt, fixed = TRUE))
})

test_that("RunSlingshot() has neither verbose nor an on.exit / tryCatch", {
  # KNOWN-ISSUE (RSS-10, P3): the timing message prints unconditionally --
  # there is no verbose switch -- and nothing guards the intermediate state.
  # Source: R/RunSlingshot.R:47-48 (message) / the function body throughout.
  f <- formals(SeuratExtend::RunSlingshot)
  expect_false("verbose" %in% names(f))
  body_txt <- paste(deparse(body(SeuratExtend::RunSlingshot)), collapse = "\n")
  expect_false(grepl("on.exit", body_txt, fixed = TRUE))
})
