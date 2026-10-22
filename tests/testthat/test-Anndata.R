# Test file for R/Anndata.R — Python-integration entry points
# Covers the R-side validation of the Seurat<->AnnData bridge.
#
# Regression issue #53: Seu2Adata()/adata.AddDR() must give a *clear R error*
# when the requested reduction is missing — not a cryptic Python NameError
# ("name 'r_umap' is not defined") raised later inside the Python runtime.
#
# NOTE on the gate: both cases below used to skip whenever a bare Python
# interpreter existed, which made them skip on exactly the machines where the
# bug they guard is observable -- the gate and the defect lived in the same
# place. They now run everywhere and branch on the environment instead
# (checklist AND-2). The environment branch below probes whether the
# 'seuratextend' conda env exists and scanpy imports inside it
# (python_module_available_for_test()), not whether this session happens to
# have an initialized interpreter.
#
# The heavy python round trips (Seu2Adata full run) are still NOT executed
# here; only the R-side guard is exercised. (conda_available() comes from
# helper-testthat.R)

test_that("adata.AddDR validates the Python environment before it validates 'dr'", {
  data(pbmc)
  # KNOWN-ISSUE (AND-2, P1) -- the issue #53 guard sits at R/Anndata.R:134,
  # but load_condaenv() runs earlier, at :125. When the Python environment is
  # unusable the call therefore dies with an environment error raised by
  # load_condaenv() before the #53 guard is reached -- an error that never
  # mentions the reduction the caller asked for, so the friendly
  # "Dimension reduction(s) ... not found" message is unreachable exactly
  # where it would be most useful.
  err <- tryCatch({
    adata.AddDR(pbmc, dr = "r_fake")
    "no condition raised"
  }, error = function(e) conditionMessage(e))

  if (conda_available() && python_module_available_for_test("scanpy")) {
    # environment is ready -> the issue #53 guard is reached and speaks up
    expect_match(err, "not found")
  } else {
    # the guard is shadowed by the environment check
    expect_match(err, "conda")
  }
})

test_that("adata.AddDR never leaks a raw Python traceback for a bad reduction", {
  data(pbmc)
  # Whichever branch the environment takes (see above), the failure must be an
  # R condition that names the problem -- not a Python NameError about r_umap.
  err <- tryCatch({
    adata.AddDR(pbmc, dr = "X_not_a_reduction")
    "no condition raised"
  }, error = function(e) conditionMessage(e))
  expect_false(grepl("r_umap", err, fixed = TRUE))
  expect_match(err, "conda|not found")
})
