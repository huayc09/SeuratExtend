# Test file for R/SearchDatabase.R
# Covers: SearchDatabase, SearchPathways (offline GO/Reactome search)
#
# The search runs against the *local* SeuratExtendData database (no network
# access in these tests: only installed tables are queried).
#
# NOTE on the gates: this file used to be guarded by
#   search_db_available() -> requireNamespace("SeuratExtendData")
# which is a DEAD GATE -- SeuratExtendData is in DESCRIPTION:Depends, so the
# condition is always TRUE and the `skip_if_not(...)` calls never fired.
# Removed (checklist L2-6).
#
# NOTE on coverage: the original file had ZERO Reactome assertions and never
# called SearchPathways (it only checked that the name is exported), even
# though all three share the same grepl()-based code path. Both are covered at
# the bottom.

test_that("SearchDatabase finds pathways by gene (local table)", {
  options(spe = "human")
  res <- SearchDatabase("CD3D", type = "gene", database = "GO",
                        spe = "human", return = "ID")
  # CD3D is a T-cell marker present in immune GO terms
  expect_true(is.character(res))
  expect_true(length(res) > 0)
  expect_true(all(grepl("^GO:", res)))
})

test_that("SearchDatabase export.to.data.frame returns a data.frame", {
  options(spe = "human")
  df <- SearchDatabase("CD3D", database = "GO", spe = "human",
                       export.to.data.frame = TRUE)
  expect_s3_class(df, "data.frame")
  expect_gt(nrow(df), 0)
})

test_that("SearchPathways is exported", {
  expect_true("SearchPathways" %in% getNamespaceExports("SeuratExtend"))
})

# --- KNOWN-ISSUE pins (2026-09-16 recheck against 问题清单0915) ------------

test_that("SearchDatabase treats the search term as a regular expression", {
  options(spe = "human")
  # KNOWN-ISSUE (SD-1, P1) -- SearchDatabase.R:87, :91 (GO), :103, :107
  # (Reactome), :178 (SearchPathways) all call grepl(item, ...) WITHOUT
  # fixed = TRUE, so the user's search term is compiled as a regex. Measured
  # 2026-09-16; three separate consequences, all of them silent.

  # (a) a lone metacharacter matches EVERYTHING -- the whole database comes
  #     back with no indication that nothing was actually searched for.
  #     22144 is every human GO id; revisit if SeuratExtendData changes.
  hits <- suppressMessages(quietly(
    SearchDatabase("+", database = "GO", spe = "human", return = "ID")))
  expect_length(hits, 22144)

  # (b) a real pathway name copied verbatim out of the database returns 0 hits,
  #     because "(" opens a regex group. Its own prefix does match, which is
  #     what makes the failure indistinguishable from "not in the database".
  by_name <- function(x) suppressMessages(quietly(
    SearchDatabase(x, type = "SetName", database = "GO", spe = "human",
                   return = "ID")))
  exact <- "acylglycerone-phosphate reductase (NADP+) activity"
  expect_length(by_name(exact), 0)
  expect_length(by_name("acylglycerone-phosphate reductase"), 1)
  # Scale of the trap: of the 22144 human GO names, 483 contain "(", 54 contain
  # "[" and 126 contain "+". (NB count over the human GO ids -- GO_ontology$name
  # alone holds 48035 rows spanning species and would give 1881/178/331.)

  # (c) an unbalanced bracket aborts with a raw regex-engine message that names
  #     neither the argument nor the function -- "报错指向不明". Note it DOES
  #     echo the pattern, so the term is visible; what is missing is which
  #     argument to blame.
  tre <- character(0)
  bad <- tryCatch({
    withCallingHandlers(
      suppressMessages(quietly(
        SearchDatabase("metabolic(", type = "SetName", database = "GO",
                       spe = "human", return = "ID"))),
      warning = function(w) {
        tre <<- c(tre, conditionMessage(w)); invokeRestart("muffleWarning")
      },
      message = function(m) invokeRestart("muffleMessage"))
    "no condition raised"
  }, error = function(e) conditionMessage(e))
  expect_match(bad, "invalid regular expression")
  expect_false(grepl("\\bitem\\b", bad))
  expect_false(grepl("SearchDatabase", bad, fixed = TRUE))
  # the regex engine's own complaint also arrives as a warning
  expect_true(any(grepl("TRE pattern", tre)))
})

test_that("an unrecognised database is silently ignored", {
  options(spe = "human")
  # KNOWN-ISSUE (SD-4, P2) -- :79 `if("GO" %in% database)` and :95
  # `if("Reactome" %in% database)`: when neither matches there is no error, no
  # warning and no message, the call just falls through to the return. A user
  # who typos the database name (or writes "go") gets an empty result that is
  # indistinguishable from a genuine no-hit search.
  expect_null(SearchDatabase("CD3D", database = "KEGG", spe = "human",
                             return = "ID"))
  expect_equal(length(SearchDatabase("CD3D", database = "go", spe = "human",
                                     return = "all")), 0)
})

test_that("SearchDatabase also works against Reactome", {
  # The original file had zero Reactome coverage even though :103/:107 is the
  # same code path as :87/:91 (checklist L2-6).
  options(spe = "human")
  res <- suppressMessages(quietly(
    SearchDatabase("CD3D", type = "gene", database = "Reactome",
                   spe = "human", return = "ID")))
  expect_true(is.character(res))
  expect_true(length(res) > 0)
  expect_true(all(grepl("^R-HSA-", res)))
})

test_that("SearchPathways actually filters the supplied gene sets", {
  options(spe = "human")
  data(hall50)
  res <- suppressMessages(quietly(
    SearchPathways(genesets = hall50$human, item = "CD3D", type = "gene")))
  # a named list of the gene sets that contain CD3D
  expect_true(length(res) >= 1)
  expect_true(all(grepl("HALLMARK", names(res))))
})

test_that("SearchPathways crashes when only one gene set is supplied", {
  data(hall50)
  # KNOWN-ISSUE (SD-3, P2) -- SearchDatabase.R:178. `sapply(item, ...)` returns
  # a plain vector when the target has length 1, and apply() rejects atomic
  # vectors of any length (`dim(X) must have a positive length`). With the
  # default type (which includes "SetName") a single-gene-set list therefore
  # aborts; measured: 6/6 combinations fail, while length 2/3/50 all work.
  one <- hall50$human[1]
  expect_error(
    suppressMessages(quietly(SearchPathways(genesets = one, item = "CD3D"))),
    "dim\\(X\\) must have a positive length"
  )
})
