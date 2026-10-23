# Test file for R/import.R
#
# import.R holds the internal `import()` helper plus its `package_data`
# catalogue of optional Python/omics backends. Neither is exported.
#
# SCOPE: the catalogue's shape, the three paths of import() (already
# installed / not installed and refused / detach requested), and the exact
# return values other files in the package rely on - several call sites
# discard the return value, which is the subject of a known issue.
#
# STYLE: characterisation tests - see test-Internal.R for the convention.

test_that("package_data describes every optional backend", {
  pd <- SeuratExtend:::package_data
  expect_type(pd, "list")
  expect_true(length(pd) >= 8)
  expect_true(all(c("SeuratDisk", "slingshot",
                    "DelayedArray", "DelayedMatrixStats") %in% names(pd)))
  # every entry must carry at least a website and an installer
  expect_true(all(vapply(pd, function(x) "website" %in% names(x), TRUE)))
  expect_true(all(vapply(pd, function(x) "install" %in% names(x), TRUE)))
  # FIXED (D4, audit 03b-audit.md) in fix round of work order
  # runs/2026-10-22-deps-guard-silentfix: `slingshot` used to be the one entry
  # without `install_info`, so import_required() would fall back to
  # install.packages("slingshot") for a Bioconductor-only package. It now
  # carries the same BiocManager hint as the DelayedArray entry.
  expect_true("install_info" %in% names(pd[["slingshot"]]))
  expect_identical(pd[["slingshot"]]$install_info,
                   'BiocManager::install("slingshot")')
  # registry completeness: every Bioconductor entry must carry an
  # install_info, otherwise the helper's CRAN fallback would be wrong.
  bioc <- vapply(pd, function(x) grepl("bioconductor.org", x$website,
                                       fixed = TRUE), logical(1))
  expect_true(all(vapply(pd[bioc], function(x) "install_info" %in% names(x),
                         logical(1))))
})

test_that("import() returns TRUE for a package that is already installed", {
  # `parallel` is an R-recommended package that is not attached by default,
  # which is what makes it usable as a detach fixture below.
  was_attached <- "package:parallel" %in% search()
  expect_true(quietly(SeuratExtend:::import("stats")))
  expect_true(quietly(SeuratExtend:::import("parallel")))
  # leave the search path the way we found it so the detach test can run
  if (!was_attached && "package:parallel" %in% search()) {
    detach("package:parallel", unload = FALSE)
  }
})

test_that("import() detaches again when detach = TRUE", {
  skip_if("package:parallel" %in% search(), "parallel is already attached")
  # KNOWN-ISSUE (L2-7 AUC10 family): the detach branch passes unload = TRUE,
  # which warns whenever another loaded package imports the namespace
  # ("namespace 'parallel' is imported by 'future', ..."). The detach itself
  # still happens, so the warning is silenced here - and pinned separately
  # below - rather than being allowed to fail an unrelated assertion.
  res <- suppressWarnings(quietly(SeuratExtend:::import("parallel", detach = TRUE)))
  expect_true(res)
  expect_false("package:parallel" %in% search())
})

test_that("import(detach = TRUE) really does ask for a namespace unload", {
  body_txt <- deparse(body(SeuratExtend:::import))
  expect_true(any(grepl("unload = TRUE", body_txt, fixed = TRUE)))
})

test_that("import() of a missing package reports and returns FALSE when refused", {
  # Non-interactive sessions answer readline() with "" -> the "try install?"
  # prompt is declined. Never run this interactively: the prompt would block.
  skip_if(interactive(), "would block on readline()")
  pkg <- "seuratextend_nonexistent_probe_pkg"
  expect_message(res <- quietly(SeuratExtend:::import(pkg)), "not installed")
  expect_false(res)
})

test_that("import() of a missing package announces the skip", {
  skip_if(interactive(), "would block on readline()")
  expect_message(
    quietly(SeuratExtend:::import("seuratextend_nonexistent_probe_pkg")),
    "installation skipped"
  )
})

# KNOWN-ISSUE (L2-7 AUC10, P1; same root cause as L3-1 PY-1 / L3-4):
# `import()` signals failure only through messages and a return value that
# several call sites drop on the floor. The assertions above pin that
# behaviour so that a fix which starts throwing (or starts returning TRUE)
# shows up as a test change rather than a silent difference.
test_that("import() failure is signalled by a return value, not by an error", {
  skip_if(interactive(), "would block on readline()")
  expect_no_error(
    quietly(suppressMessages(
      SeuratExtend:::import("seuratextend_nonexistent_probe_pkg")))
  )
})

# --- import_required(): the call-site guard ---------------------------------
#
# Batch-1 of runs/2026-10-22-deps-guard-silentfix replaced the call sites that
# dropped import()'s return value with import_required(), which keeps the
# interactive install prompt and then fails fast with an install hint.

test_that("import_required() returns TRUE for an installed package", {
  res <- quietly(SeuratExtend:::import_required("ggpubr"))
  expect_identical(res, TRUE)
})

test_that("import_required() fails fast with an install hint when missing", {
  # import() answers readline() with "" in a non-interactive session, so this
  # exercises the exact path a script / knitr / CI run takes.
  skip_if(interactive(), "would block on readline()")
  pkg <- "seuratextend_nonexistent_probe_pkg"

  expect_error(
    quietly(SeuratExtend:::import_required(pkg)),
    "required but not installed"
  )

  msg <- tryCatch(
    { quietly(SeuratExtend:::import_required(pkg)); NULL },
    error = function(e) conditionMessage(e)
  )
  # non-registry packages get the truthful CRAN route, naming the package
  expect_match(msg, sprintf('install.packages\\("%s"\\)', pkg))

  # contrast: the bare helper still returns FALSE rather than throwing
  expect_false(suppressMessages(quietly(SeuratExtend:::import(pkg))))
})
