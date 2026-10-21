# Test file for R/generics.R
# Covers: the four S3 generics defined in this file and the roxygen blocks that
#         build their shared help pages -- VlnPlot2 (:63) / CalcStats (:108) /
#         WaterfallPlot (:184) / VolcanoPlot (:286). The implementations live in
#         R/Vlnplots.R, R/CalcStats.R, R/WaterfallPlot.R, R/VolcanoPlot.R; this
#         file only owns the generic layer plus the documentation.
# Regression issues: combined-audit checklist (9 items: 5x P2 / 4x P3)
# Data: built-in pbmc
#
# Unless noted otherwise, the documentation assertions read the INSTALLED Rd
# database (tools::Rd_db), so they verify what users actually see, not the
# roxygen source text. GEN-5 is the exception: it asserts source-tree doc sync
# by reading ../../man/*.Rd directly, because the installed package lags the
# working tree under devtools::load_all()/test().

rd_tag_of <- function(x) {
  t <- attr(x, "Rd_tag")
  if (is.null(t)) "" else t
}

rd_sections <- function(rd_name) {
  db <- tryCatch(tools::Rd_db("SeuratExtend"), error = function(e) NULL)
  if (is.null(db) || !rd_name %in% names(db)) {
    skip(paste("installed Rd database has no", rd_name))
  }
  tags <- vapply(db[[rd_name]], rd_tag_of, "")
  list(rd = db[[rd_name]], tags = tags)
}

rd_flat <- function(rd_name, tag) {
  s <- rd_sections(rd_name)
  seg <- s$rd[s$tags == tag]
  if (length(seg) == 0) return(NA_character_)
  out <- character(0)
  walk <- function(z) {
    if (is.character(z)) out <<- c(out, z)
    else if (is.list(z)) for (e in z) walk(e)
    invisible(NULL)
  }
  walk(seg)
  paste(out, collapse = "")
}

rd_arg_items <- function(rd_name) {
  # named list: argument -> its description text
  s <- rd_sections(rd_name)
  seg <- s$rd[s$tags == "\\arguments"]
  if (length(seg) == 0) return(list())
  out <- list()
  for (item in seg[[1]]) {
    if (identical(rd_tag_of(item), "\\item") && length(item) >= 2) {
      nm <- paste0(unlist(item[[1]]), collapse = "")
      out[[nm]] <- paste0(unlist(item[[2]]), collapse = "")
    }
  }
  out
}

rd_arg_names <- function(rd_name) {
  names(rd_arg_items(rd_name))
}

# --- the generic / method argument mismatch -------------------------------

test_that("the four generics take `object` while their methods take `seu`", {
  # KNOWN-ISSUE (GEN-1, P2): every generic is function(object, ...) but every
  # S3 method starts with seu (or matr for VlnPlot2.default). UseMethod passes
  # arguments BY NAME, so the two names never meet.
  # Source: R/generics.R:63/:108/:184/:286 (signatures), :64/:109/:185/:287
  # (UseMethod); methods at R/Vlnplots.R:20, R/CalcStats.R:15,
  # R/WaterfallPlot.R:18, R/VolcanoPlot.R:17.
  ns <- asNamespace("SeuratExtend")
  for (f in c("VlnPlot2", "CalcStats", "WaterfallPlot", "VolcanoPlot")) {
    expect_identical(names(formals(get(f, envir = ns)))[1], "object", info = f)
  }
  for (m in c("VlnPlot2.Seurat", "CalcStats.Seurat",
              "WaterfallPlot.Seurat", "VolcanoPlot.Seurat")) {
    expect_true("seu" %in% names(formals(get(m, envir = ns))), info = m)
    expect_false("object" %in% names(formals(get(m, envir = ns))), info = m)
  }
  expect_identical(names(formals(get("VlnPlot2.default", envir = ns)))[1], "matr")
})

test_that("naming the first argument after the docs or the method both fail", {
  # KNOWN-ISSUE (GEN-1, P2) -- behavioural half. Positional calls work; both
  # documented spellings do not. The docs themselves disagree: the usage line
  # says VlnPlot2(object, ...) while the method entry says (seu, ...).
  # Source: R/generics.R:63/:108.
  data(pbmc)
  expect_no_error(VlnPlot2(pbmc, features = "CD3D"))
  expect_no_error(CalcStats(pbmc, features = "CD3D"))

  expect_error(VlnPlot2(object = pbmc, features = "CD3D"),
               'argument "seu" is missing')
  expect_error(VlnPlot2(seu = pbmc, features = "CD3D"),
               'argument "object" is missing')
  expect_error(CalcStats(object = pbmc, features = "CD3D"),
               "unused argument")
})

# --- the VolcanoPlot documentation block ----------------------------------

test_that("VolcanoPlot's help page documents neither `object` nor `...`", {
  # KNOWN-ISSUE (GEN-2, P2): the roxygen block at R/generics.R:189-285 carries
  # no @param at all. The Rd does have an \arguments section (the methods in
  # R/VolcanoPlot.R contribute theirs) but it starts at `seu` -- the generic's
  # own two formals are undocumented, which is what checkDocFiles reports as
  # "missing: object, ...".
  args <- rd_arg_names("VolcanoPlot.Rd")
  expect_gt(length(args), 0)                 # the section exists
  expect_false("object" %in% args)
  expect_false("..." %in% args)
  expect_identical(args[1], "seu")
})

test_that("VolcanoPlot's help page has no \\value section at all", {
  # KNOWN-ISSUE (GEN-3, P2): neither R/generics.R:189-285 nor R/VolcanoPlot.R
  # carries an @return, so the merged Rd never gets a \value section.
  # Source: R/generics.R:189-285.
  s <- rd_sections("VolcanoPlot.Rd")
  expect_false("\\value" %in% s$tags)
  # contrast: the other three all have one
  for (nm in c("VlnPlot2.Rd", "CalcStats.Rd", "WaterfallPlot.Rd")) {
    expect_true("\\value" %in% rd_sections(nm)$tags, info = nm)
  }
})

test_that("CalcStats documents `object` but not `...`", {
  # KNOWN-ISSUE (GEN-4, P2): VlnPlot2 / WaterfallPlot also omit @param ... in
  # the generic block, but their default methods have `...` and their method
  # @param ... is merged into the same Rd, which fills the gap. CalcStats has
  # no `...` in either method, so nothing does -- and it shows.
  # Source: R/generics.R:67-107 (only @param object and @return) / :108.
  args <- rd_arg_names("CalcStats.Rd")
  expect_true("object" %in% args)
  expect_false("..." %in% args)
  # contrast: the two that do get filled in
  expect_true("..." %in% rd_arg_names("VlnPlot2.Rd"))
  expect_true("..." %in% rd_arg_names("WaterfallPlot.Rd"))
})

# --- the VlnPlot2 help page ------------------------------------------------

# Source-tree Rd reader for GEN-5: dev-layout gated, environment-independent.
rd_source_usage <- function(rd_name) {
  path <- file.path("..", "..", "man", rd_name)
  skip_if_not(file.exists(path),
              paste("source Rd not available:", path))
  lines <- readLines(path, warn = FALSE)
  start <- grep("^\\\\usage\\{", lines)
  if (length(start) == 0) return(NA_character_)
  start <- start[1]
  end <- which(seq_along(lines) > start & lines == "}")[1]
  paste(lines[start:end], collapse = "\n")
}

test_that("VlnPlot2's usage section now includes the `sort` argument", {
  # RESOLVED (GEN-5, P2): the Rd drift is gone -- man/VlnPlot2.Rd and
  # man/StackedViolin.Rd were regenerated in the external-matrix work order and
  # now list `sort`, matching the method signature.
  # Source: man/VlnPlot2.Rd + man/StackedViolin.Rd usage blocks.
  usage <- rd_source_usage("VlnPlot2.Rd")
  expect_true(grepl("\\bsort\\b", usage))
  # the neighbours it is sandwiched between are still present
  expect_true(grepl("\\bviolin\\b", usage))
  expect_true(grepl("\\bbox\\b", usage))

  # the StackedViolin page embeds the VlnPlot2 method usage too, and is synced
  stacked_usage <- rd_source_usage("StackedViolin.Rd")
  expect_true(grepl("\\bsort\\b", stacked_usage))

  # contrast: the code really has it
  expect_true("sort" %in% names(formals(
    get("VlnPlot2.Seurat", envir = asNamespace("SeuratExtend")))))
})

test_that("VlnPlot2's value section is written twice", {
  # KNOWN-ISSUE (GEN-6, P3): two @return tags (:4 and :5) collapse into one
  # \value section holding both sentences.
  # Source: R/generics.R:4-5.
  v <- rd_flat("VlnPlot2.Rd", "\\value")
  expect_true(grepl("ggplot object", v))
  expect_true(grepl("A ggplot object.", v, fixed = TRUE))
  expect_true(grepl("ggplot object", v))
  # contrast: the others say it once
  expect_identical(gsub("\\s+", " ", trimws(rd_flat("CalcStats.Rd", "\\value"))),
                   "A matrix.")
})

test_that("VlnPlot2's description and its object argument disagree about input", {
  # KNOWN-ISSUE (GEN-7, P3): the description accepts "a Seurat object or a
  # data frame" while @param object says "a Seurat object or matrix", and the
  # method-side @param matr says "matrix or data frame" -- three wordings.
  # Source: R/generics.R:2 vs :3.
  d <- rd_flat("VlnPlot2.Rd", "\\description")
  items <- rd_arg_items("VlnPlot2.Rd")
  expect_true(grepl("data frame", d))
  expect_true(grepl("matrix", items[["object"]]))
  # the generic's own @param object never mentions data frames, while the
  # method-side @param matr does -- hence the three-way wording drift
  expect_false(grepl("data frame", items[["object"]]))
  expect_true(grepl("data frame", items[["matr"]]))
})

test_that("the four generics ship examples with no \\donttest", {
  # KNOWN-ISSUE (GEN-8, P3): R CMD check runs these by default. One of them
  # (VolcanoPlot(pbmc) at R/generics.R:231) ranks all 12627 genes by default,
  # and GeneSetAnalysis(pbmc, hall50$human) measures ~6.5 s cold.
  # Source: R/generics.R:55/:98/:139/:231/:270.
  for (nm in c("VlnPlot2.Rd", "CalcStats.Rd", "WaterfallPlot.Rd", "VolcanoPlot.Rd")) {
    ex <- rd_flat(nm, "\\examples")
    expect_false(is.na(ex), info = nm)
    expect_false(grepl("donttest", ex), info = nm)
    expect_false(grepl("dontrun", ex), info = nm)
  }
})

# --- the `cell` / `cells` inconsistency ------------------------------------

test_that("cell= reaches the Seurat methods by partial matching only", {
  # KNOWN-ISSUE (GEN-9, P3): R/generics.R uses `cells = cells` at :43 and
  # `cell = cells` at :46/:50 in the same example block. It works today
  # because VlnPlot2.Seurat has `cells` (so `cell=` partially matches) and
  # WaterfallPlot.Seurat/VolcanoPlot.Seurat have a real `cell` formal -- a new
  # method argument starting with "cell" would turn the first case ambiguous.
  # Source: R/generics.R:43/:46/:50.
  data(pbmc)
  cells <- colnames(pbmc)[1:50]
  expect_no_error(VlnPlot2(pbmc, features = "CD3D", cells = cells))
  expect_no_error(VlnPlot2(pbmc, features = "CD3D", cell = cells))
  expect_true("cells" %in% names(formals(
    get("VlnPlot2.Seurat", envir = asNamespace("SeuratExtend")))))
  expect_false("cell" %in% names(formals(
    get("VlnPlot2.Seurat", envir = asNamespace("SeuratExtend")))))
})
