# Test file for R/GenesymbolConversion.R
# Covers: HumanToMouseGenesymbol, MouseToHumanGenesymbol, EnsemblToGenesymbol,
#         GenesymbolToEnsembl, UniprotToGenesymbol
# Data: built-in pbmc gene space; conversion tables ship with the package
#       (data/mouse_human_genesymbols.rda etc.). local.mode=TRUE uses these
#       offline tables (no BioMart / no network) by default.

# --- human <-> mouse orthologues ------------------------------------------

test_that("HumanToMouseGenesymbol converts known orthologues (local mode)", {
  out <- HumanToMouseGenesymbol(c("CD14", "CD3D", "CD79A"))
  # default match=TRUE -> data.frame of the conversion table rows
  expect_s3_class(out, "data.frame")
  expect_equal(sort(out$HGNC.symbol), sort(c("CD14", "CD3D", "CD79A")))
  expect_equal(sort(out$MGI.symbol), sort(c("Cd14", "Cd3d", "Cd79a")))
})

test_that("MouseToHumanGenesymbol converts known orthologues (local mode)", {
  out <- MouseToHumanGenesymbol(c("Cd14", "Cd3d", "Cd79a"))
  expect_s3_class(out, "data.frame")
  expect_equal(sort(out$MGI.symbol), sort(c("Cd14", "Cd3d", "Cd79a")))
  expect_equal(sort(out$HGNC.symbol), sort(c("CD14", "CD3D", "CD79A")))
})

test_that("MouseToHumanGenesymbol match=FALSE returns bare vector", {
  out <- MouseToHumanGenesymbol(c("Cd14", "Cd3d", "Cd79a"), match = FALSE)
  expect_true(is.character(out))
  expect_setequal(out, c("CD14", "CD3D", "CD79A"))
})

# --- Ensembl conversions (local tables) -----------------------------------

test_that("GenesymbolToEnsembl returns one-to-one IDs with keep.seq", {
  out <- GenesymbolToEnsembl(c("CD14", "CD3D"), spe = "human",
                             keep.seq = TRUE)
  expect_type(out, "character")
  expect_named(out, c("CD14", "CD3D"))
  expect_true(all(grepl("^ENSG", out)))
})

test_that("EnsemblToGenesymbol converts back to symbols", {
  ids <- c("ENSG00000170458", "ENSG00000167286")  # CD14, CD3D
  # default: data.frame with ensembl id + symbol columns
  out <- EnsemblToGenesymbol(ids, spe = "human")
  expect_s3_class(out, "data.frame")
  expect_true(all(c("ensembl_gene_id", "hgnc_symbol") %in% colnames(out)))
  expect_true(all(c("CD14", "CD3D") %in% out$hgnc_symbol))
  # keep.seq = TRUE -> named character, one-to-one
  out2 <- EnsemblToGenesymbol(ids[1], spe = "human", keep.seq = TRUE)
  expect_type(out2, "character")
  expect_equal(unname(out2), "CD14")
})

# --- Uniprot conversion ----------------------------------------------------

test_that("UniprotToGenesymbol maps known human entries", {
  out <- UniprotToGenesymbol(c("P08571", "P04234"), spe = "human")
  expect_true(is.character(out))
  # P08571 = CD14, P04234 = CD3D
  expect_setequal(unname(out), c("CD14", "CD3D"))
})

test_that("UniprotToGenesymbol reports unknown ids", {
  expect_message(
    UniprotToGenesymbol("P00000_UNKNOWN", spe = "human"),
    regexp = "not found"
  )
})

# --- KNOWN-ISSUE pins (2026-09-16 recheck against 问题清单0915) ------------

test_that("the converters overflow the stack on a factor", {
  # KNOWN-ISSUE (GSC-1, P1) -- R/GenesymbolConversion.R:240 dispatches on
  # is.vector(item); is.vector(factor(...)) is FALSE, so a factor takes the
  # else branch at :286-288 and recurses on rownames(item), which is NULL.
  # A factor is exactly what you get from a data.frame column, and
  # as.matrix(<vector>) hits the same path. (expect_stack_overflow lives in
  # helper-testthat.R.)
  expect_stack_overflow(HumanToMouseGenesymbol(factor("CD14")))
})

test_that("the converters overflow the stack on an unnamed matrix", {
  expect_stack_overflow(HumanToMouseGenesymbol(as.matrix(1:4)))
})

test_that("the converters work on a named matrix (control for the two above)", {
  m <- matrix(1:4, 2, dimnames = list(c("CD14", "CD3D"), NULL))
  expect_no_error(MouseToHumanGenesymbol(m, match = FALSE))
})

test_that("keep.seq = TRUE leaves NAs behind when the input has duplicates", {
  # KNOWN-ISSUE (GSC-3, P2) -- R/GenesymbolConversion.R:266-283 (especially
  # :278). The assignment `genelist2[glf[[col.source]]] <- glf[[col.target]]`
  # indexes by NAME, and `x["name"] <- v` writes only the FIRST match, so every
  # duplicate after the first is left as NA. The documentation at :10 promises
  # a one-to-one mapping "when precise correspondence is required".
  out <- HumanToMouseGenesymbol(c("CD14", "CD14", "CD3D"), keep.seq = TRUE)
  expect_length(out, 3)
  expect_equal(as.character(out[1]), "Cd14")
  expect_true(is.na(out[2]))               # the duplicate silently lost
  expect_equal(as.character(out[3]), "Cd3d")
})

test_that("an illegal 'spe' is not validated", {
  # KNOWN-ISSUE (GSC-2, P2) -- R/GenesymbolConversion.R:124, :139, :168, :183.
  # check_spe() only tests is.null(spe), so a mistyped species name reaches
  # `mappings[[spe]]` -> NULL -> match.table[logical(0), ] -> nrow(NULL) ->
  # the opaque "argument is of length zero".
  expect_error(EnsemblToGenesymbol("ENSG00000170458", spe = "bogus"),
               "argument is of length zero")
})
