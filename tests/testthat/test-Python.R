# Test file for R/Python.R
#
# Python.R is the conda/reticulate bridge: load_condaenv(),
# create_condaenv_seuratextend(), activate_python() (all exported) plus the
# string/path helpers is_windows(), normalizePathWindows(), create_temp_dir()
# and r_vector_to_py() that the Python-facing wrappers build their calls with.
#
# SCOPE: the pure-R half - which is where the interesting failure modes live,
# because every Python-facing wrapper pastes R strings into Python source.
# The conda-touching half is skipped unless a conda installation is found.
#
# KNOWN-ISSUE references: the two-mode quoting rule below was established for
# the L3-5 audit (Cellrank) and applies to every consumer of
# r_vector_to_py() / shQuote(type = "sh") in this package.

py_helpers_internal <- c("load_condaenv", "is_windows", "normalizePathWindows",
                         "create_temp_dir", "r_vector_to_py")

test_that("only the user-facing conda helpers are exported", {
  exports <- getNamespaceExports("SeuratExtend")
  expect_true("create_condaenv_seuratextend" %in% exports)
  expect_true("activate_python" %in% exports)
  expect_length(intersect(py_helpers_internal, exports), 0)
})

test_that("conda helper signatures and defaults are stable", {
  expect_identical(names(formals(SeuratExtend::create_condaenv_seuratextend)), "force")
  expect_false(eval(formals(SeuratExtend::create_condaenv_seuratextend)$force))
  expect_identical(
    names(formals(SeuratExtend::activate_python)),
    c("conda_env", "verbose", "packages")
  )
  expect_identical(formals(SeuratExtend::activate_python)$conda_env, "seuratextend")
  expect_identical(eval(formals(SeuratExtend::activate_python)$packages), "os")
})

# --- r_vector_to_py: the R-to-Python source builder ----------------------

test_that("r_vector_to_py() renders numeric and character vectors", {
  rv <- SeuratExtend:::r_vector_to_py
  expect_identical(rv(c(1, 2, 3)), "[1, 2, 3]")
  expect_identical(rv(c("a", "b")), "['a', 'b']")
  expect_identical(rv(c(1, 2), type = "tuple"), "(1, 2)")
  expect_identical(rv(c("a", "b"), type = "tuple"), "('a', 'b')")
  expect_null(rv(NULL))
  expect_null(rv(NULL, type = "tuple"))
})

test_that("a value with neither quotes nor backslashes gets plain single quotes", {
  rv <- SeuratExtend:::r_vector_to_py
  out <- rv("umap")
  expect_identical(out, "['umap']")
  expect_true(startsWith(out, "['"))
})

test_that("a Windows-style path keeps its backslashes unescaped (single-quote mode)", {
  # KNOWN-ISSUE (L3-5, P2 family): shQuote(type = "sh") switches to
  # single-quote mode whenever the value contains no single quote, and in that
  # mode it copies backslashes through verbatim. Python then reads "\\d" etc.
  # as escape sequences. The trigger is exactly "backslash and no apostrophe".
  rv <- SeuratExtend:::r_vector_to_py
  out <- rv("C:\\data\\a")
  expect_true(startsWith(out, "['"))              # single-quote mode
  expect_false(grepl("\\\\", out, fixed = TRUE))  # backslashes NOT doubled
  expect_true(grepl("\\data", out, fixed = TRUE))
})

test_that("a value containing an apostrophe switches to double-quote mode", {
  rv <- SeuratExtend:::r_vector_to_py
  out <- rv("a'b")
  expect_true(startsWith(out, '["'))   # opened with a double quote...
  expect_true(endsWith(out, '"]'))     # ...and closed with one
})

test_that("apostrophe + backslash is the safe combination", {
  # The sufficient-and-necessary condition for the L3-5 defect is
  # "contains a backslash AND contains no single quote". Adding an apostrophe
  # flips shQuote into double-quote mode, which doubles the backslashes.
  rv <- SeuratExtend:::r_vector_to_py
  risky <- rv("C:\\data\\a")     # no apostrophe -> backslashes verbatim
  safe  <- rv("C:\\data\\a'b")   # apostrophe    -> backslashes doubled
  expect_false(grepl("\\\\", risky, fixed = TRUE))
  expect_true(grepl("\\\\", safe, fixed = TRUE))
  expect_true(startsWith(safe, '["'))
})

# --- path / temp-dir helpers ---------------------------------------------

test_that("is_windows() agrees with the running platform", {
  expect_type(SeuratExtend:::is_windows(), "logical")
  expect_length(SeuratExtend:::is_windows(), 1)
  # the source returns the named Sys.info() comparison unchanged, so compare
  # like for like rather than assuming a bare logical
  expect_identical(unname(SeuratExtend:::is_windows()),
                   unname(Sys.info()["sysname"]) == "Windows")
})

test_that("normalizePathWindows() hands back a forward-slash path", {
  p <- SeuratExtend:::normalizePathWindows("C:\\a\\b")
  expect_type(p, "character")
  expect_false(grepl("\\", p, fixed = TRUE))
})

test_that("create_temp_dir() makes a timestamped Palantir scratch directory", {
  d <- SeuratExtend:::create_temp_dir()
  on.exit(unlink(d, recursive = TRUE), add = TRUE)
  expect_true(dir.exists(d))
  expect_true(startsWith(basename(d), "PalantirPs_"))
  # it must live under tempdir(), not in the working directory
  expect_true(startsWith(normalizePath(dirname(d), winslash = "/"),
                         normalizePath(tempdir(), winslash = "/")))
})

test_that("the bundled conda environment recipes are present", {
  ex <- system.file("extdata", package = "SeuratExtend")
  for (yml in c("environment-windows.yml", "environment-linux.yml",
                "environment-mac.yml", "environment-mac-silicon.yml")) {
    expect_true(file.exists(file.path(ex, yml)), info = yml)
  }
})

# --- conda-touching paths (skipped without a conda installation) ----------

test_that("load_condaenv() refuses a conda_env that does not exist", {
  skip_if_not(conda_available(), "no conda installation visible to reticulate")
  expect_error(
    SeuratExtend:::load_condaenv("seuratextend_probe_env_that_does_not_exist"),
    "does not exist"
  )
})

test_that("every Python-facing wrapper defaults to the same conda_env name", {
  wrapped <- c("Palantir.RunDM", "Palantir.Pseudotime", "Palantir.Magic",
               "Cellrank.Compute", "Cellrank.Plot", "scVelo.SeuratToAnndata")
  for (fn in wrapped) {
    f <- get(fn, envir = asNamespace("SeuratExtend"))
    expect_identical(formals(f)$conda_env, "seuratextend", info = fn)
  }
})
