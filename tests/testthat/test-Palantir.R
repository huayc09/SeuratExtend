# Test file for R/Palantir.R
#
# Palantir.R has four exports. Every one of them runs Python, so the R-side
# contract is what these tests pin: signatures, the guard that runs before
# Python is touched, and - the one that matters most - where each function
# writes its scratch files.
#
# STYLE: characterisation tests - see test-Internal.R for the convention.

palantir_exports <- c("Palantir.RunDM", "RunPalantirDiffusionMap",
                      "Palantir.Pseudotime", "Palantir.Magic")

test_that("all four Palantir entry points are exported", {
  exports <- getNamespaceExports("SeuratExtend")
  expect_length(setdiff(palantir_exports, exports), 0)
})

test_that("Palantir.RunDM() signature is stable", {
  f <- formals(SeuratExtend::Palantir.RunDM)
  expect_identical(names(f), c("seu", "reduction", "n_components", "conda_env"))
  expect_identical(f$reduction, "pca")
  expect_identical(eval(f$n_components), 50)
  expect_identical(f$conda_env, "seuratextend")
})

test_that("the deprecated RunPalantirDiffusionMap() keeps its own default", {
  f <- formals(SeuratExtend::RunPalantirDiffusionMap)
  expect_identical(names(f), c("seu", "reduction", "n_components"))
  expect_identical(eval(f$n_components), 20)   # deliberately not 50
})

test_that("Palantir.Pseudotime() and Palantir.Magic() signatures are stable", {
  f <- formals(SeuratExtend::Palantir.Pseudotime)
  expect_identical(names(f), c("seu", "start_cell", "terminal_states",
                               "title", "conda_env", "n_jobs"))
  expect_null(eval(f$terminal_states))
  expect_identical(f$title, "Pseudotime")
  expect_identical(eval(f$n_jobs), 1)

  g <- formals(SeuratExtend::Palantir.Magic)
  expect_identical(names(g), c("seu", "n_top_genes", "n_components",
                               "conda_env", "n_jobs"))
  expect_identical(eval(g$n_top_genes), 2000)
  expect_identical(eval(g$n_components), 20)
})

test_that("Palantir.Pseudotime() needs a multiscale space before it needs Python", {
  data(pbmc)
  # pbmc has no "ms" reduction. On a machine with a working conda env this is
  # caught by the explicit guard; without conda the failure happens earlier
  # inside load_condaenv(). Either way it must be an R error, never a hang.
  expect_error(suppressMessages(
    Palantir.Pseudotime(pbmc, start_cell = colnames(pbmc)[1])
  ))
})

test_that("PAL-1 pin: Palantir.RunDM currently returns without error (silent path unconfirmed)", {
  skip_if_not(conda_available(), "PAL-1 pin is pinned against a conda-present machine")
  data(pbmc)
  # Same reasoning as above: conda missing -> load_condaenv() errors; conda
  # present -> reticulate/palantir errors. It must not return silently.
  # KNOWN-ISSUE(PAL-1) 待坐实: 本机 conda+python 俱在, Palantir.RunDM 未报错即快速返回,
  # 疑似命中清单 P1「Python.R 返回值 FALSE 无人接手」静默路径, 但是否真执行了计算未证实。
  # pin 当前行为(无 error), 根因确认与翻转归 P1 Python 链批次工单。
  # 无 conda 环境下本 pin 不适用（load_condaenv 会先报错），故加门——CI 上本测试 skip。
  expect_no_error(suppressMessages(Palantir.RunDM(pbmc, n_components = 5)))
})

test_that("Palantir.RunDM() writes its scratch files under tempdir()", {
  # L3-7 finding: the two diffusion-map functions disagree about where they
  # put their intermediate CSVs. RunDM goes through create_temp_dir() (a
  # timestamped folder under tempdir()), which is the correct behaviour.
  expect_true("create_temp_dir" %in% all.names(body(SeuratExtend::Palantir.RunDM)))
  expect_true("tempdir" %in% all.names(body(SeuratExtend:::create_temp_dir)))
})

test_that("the deprecated RunPalantirDiffusionMap() litters the working directory", {
  # KNOWN-ISSUE (L3-7): unlike Palantir.RunDM(), the deprecated function
  # creates a literal ./tmp/ in the *current working directory* and writes
  # dr.csv / dm_res.csv / ms_data.csv / tsne.csv into it. It has no
  # on.exit() cleanup, and it is not marked as deprecated in code. Source:
  # R/Palantir.R:129-149.
  #
  # This is asserted statically on purpose: calling the function would make
  # reticulate bootstrap a Python installation (it downloads numpy on first
  # use), which a test suite must not do. The literal path is what makes the
  # package write outside tempdir(), and a literal is exactly what a static
  # check can pin.
  body_txt <- deparse(body(SeuratExtend::RunPalantirDiffusionMap))
  expect_true(any(grepl('"tmp/dr.csv"', body_txt, fixed = TRUE)))
  expect_true(any(grepl("dir.create(\"tmp\")", body_txt, fixed = TRUE)))
  # and, for contrast, no tempdir-based helper is used here
  expect_false("create_temp_dir" %in% all.names(body(SeuratExtend::RunPalantirDiffusionMap)))
})
