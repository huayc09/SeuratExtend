#' @include GeneSetAnalysis.R
#'
NULL

#' @param parent ID or name of parent (top-level) geneset in GO/Reactome database,
#' so that only the subset of genesets will be calculated. Default: NULL
#' @param dataset (GO) Alias of 'parent', Default: NULL
#' @param root (GO) Which root to use (BP - Biological process, MF - Molecular Function,
#' CC - Cellular Component), Default: 'BP'
#' @param spe Species (human or mouse), Default: getOption("spe")
#' @param only.end.terms If only the end-level pathways/genesets are used for calculation, Default: F
#' @rdname GeneSetAnalysis
#' @export

GeneSetAnalysisGO <- function(
    seu = NULL,
    dataset = NULL, parent = NULL, root = "BP", spe = getOption("spe"),
    ratio = 0.4, n.min = 1, n.max = Inf, only.end.terms = F,
    slot = "counts", assay = "RNA", nCores = 1,
    aucMaxRank = NULL,
    title = NULL,
    export_to_matrix = F, verbose = TRUE,
    n.items.part = NULL){

  check_spe(spe)
  library(SeuratObject)
  library(SeuratExtendData)

  if(is.null(seu)){
    message("Commonly used datasets:\n  ",
            paste(DatabaseList, names(DatabaseList), sep = ": " ,collapse = "\n  "))
    return(DatabaseList)
  }

  parent <- union(dataset, parent)
  if(verbose) message(Sys.time(), " Retrieve GO database")
  GenesetNames <- getGOdatabase(
    parent = parent,
    root = root,
    spe = spe)
  if(is.null(title)) {
    if(is.null(parent)) title <- root[1] else title <- parent[1]
  }

  GenesetList <- SeuratExtendData::GO_Data[[spe]]$GO2Gene[GenesetNames]
  DefaultAssay(seu) <- assay
  GenesetList <- FilterGenesets(
    genes_in_data = rownames(seu),
    genesets = GenesetList,
    ratio = ratio,
    n.min = n.min,
    n.max = n.max,
    verbose = verbose)

  if(only.end.terms) {
    GenesetNames <- GetEndTermsGO(names(GenesetList), spe)
    GenesetList <- GenesetList[GenesetNames]
  }

  seu <- BuildAUCRank(
    seu,
    slot = slot,
    assay = assay,
    verbose = verbose)

  if(is.null(aucMaxRank)) aucMaxRank <- ceiling(0.05 * nrow(seu))
  AUC_matrix <- calcAUC_matrix(
    GenesetList = GenesetList,
    rankings = seu@misc$AUCell$cells_rankings,
    nCores = nCores,
    aucMaxRank = aucMaxRank,
    verbose = verbose,
    n.items.part = n.items.part)

  if(export_to_matrix) return(AUC_matrix)
  seu@misc[["AUCell"]][["GO"]][[title]] <- AUC_matrix
  return(seu)
}

#' @title Filter GO/Reactome genesets
#' @description Filter GO/Reactome genesets IDs
#' @param term GO/Reactome IDs, can be character or matrix (ID in rownames)
#' @param spe Species (human or mouse), Default: getOption("spe")
#' @param n.min Min number of genes in the geneset, Default: 1
#' @param n.max Max number of genes in the geneset, Default: Inf
#' @param only.end.terms If only the end-level pathways/genesets are used for calculation, Default: F
#' @param change.name If to change names, Default: F
#' @param parent ID or name of parent (top-level) geneset in GO/Reactome database, Default: NULL
#' @return character or matrix, depends on the input
#' @examples
#' options(max.print = 100, spe = "human")
#'
#' # Filter GO terms
#' terms <- names(SeuratExtendData::GO_Data$human$GO2Gene)
#' terms
#' terms <- FilterGOTerms(terms, n.max = 1000, n.min = 10, only.end.terms = T, parent = "GO:0002376")
#' RenameGO(terms)
#'
#' # Filter Reactome terms
#' terms <- names(SeuratExtendData::Reactome_Data$human$Path2Gene)
#' terms
#' terms <- FilterReactomeTerms(terms, n.max = 1000, n.min = 10, only.end.terms = T, parent = "R-HSA-168256")
#' RenameReactome(terms)
#' @rdname FilterGOTerms
#' @export

FilterGOTerms <- function(
    term,
    spe = getOption("spe"),
    n.min = 1,
    n.max = Inf,
    only.end.terms = T,
    change.name = F,
    parent = NULL){
  check_spe(spe)
  library(SeuratExtendData)
  if(is.vector(term)){
    filter <-
      SeuratExtendData::GO_Data[[spe]]$GO2Gene[term] %>%
      sapply(function(x) length(x) >= n.min & length(x) <= n.max)
    term <- term[filter]
    if(!is.null(parent)){
      if(all(parent %in% SeuratExtendData::GO_Data[[spe]]$GO_ontology$id)) {
        term <- intersect(term, GetAllChilrenGO(parent, spe))
      }else{
        message <- paste0(
          "Parent item(s) not found: ",
          paste0(setdiff(parent, SeuratExtendData::GO_Data[[spe]]$GO_ontology$id),
                 collapse = ", "))
        warning(message)
      }
    }
    if(only.end.terms) term <- GetEndTermsGO(term, spe)
    if(length(term) == 0) stop("No item(s) passed the filter")
  }else{
    term <- term[FilterGOTerms(rownames(term), spe = spe, n.min = n.min, n.max = n.max,
                               only.end.terms = only.end.terms, change.name = F, parent = parent),]
  }
  if(change.name) term <- RenameGO(term)
  return(term)
}

#' @title Add full name to GO/Reactome pathway ID
#' @description Rename GO/Reactome ID to full geneset/pathway name
#' @param term GO/Reactome IDs, can be character or matrix (ID in rownames)
#' @param add_id If show IDs, Default: T
#' @param add_n_gene If show how many genes in each geneset, Default: T
#' @param spe Species (human or mouse), Default: getOption("spe")
#' @return character or matrix, depends on the input
#' @details Also see examples in \code{\link[SeuratExtend:GeneSetAnalysis]{GeneSetAnalysis()}}
#' @examples
#' RenameGO(c("GO:0002376","GO:0050896"), spe = "human")
#'
#' RenameReactome(c("R-HSA-109582","R-HSA-112316"), spe = "human")
#' @rdname RenameGO
#' @export

RenameGO <- function(
    term,
    add_id = T,
    add_n_gene = T,
    spe = getOption("spe")){
  if(is.vector(term)){
    check_spe(spe)
    library(SeuratExtendData)
    renamed <- SeuratExtendData::GO_Data[[spe]]$GO_ontology$name[term]
    if(add_id) renamed <- paste(term, renamed)
    if(add_n_gene) renamed <- paste0(
      renamed, " (",
      sapply(SeuratExtendData::GO_Data[[spe]]$GO2Gene[term], length), "g)")
    return(renamed)
  }else{
    library(magrittr)
    library(dplyr)
    term <- term %>% set_rownames(RenameGO(rownames(.), add_id = add_id, add_n_gene = add_n_gene))
    return(term)
  }
}

# Internal ----------------------------------------------------------------

GetChilrenGO <- function(term, spe = getOption("spe")){
  check_spe(spe)
  Child <- union(term, unlist(SeuratExtendData::GO_Data[[spe]]$GO_ontology$children[term]))
  return(Child)
}

GetAllChilrenGO <- function(term, spe = getOption("spe")){
  Child <- GetChilrenGO(term, spe)
  if(length(Child) > length(term)) return(GetAllChilrenGO(Child, spe)) else return(term)
}

GetEndTermsGO <- function(term, spe = getOption("spe")){
  library(dplyr)
  check_spe(spe)
  anc <- SeuratExtendData::GO_Data[[spe]]$GO_ontology$ancestors[term] %>%
    lapply(function(x) x[-length(x)]) %>%
    unlist() %>% unique()
  return(setdiff(term, anc))
}

DatabaseList <- c(
  "immune_system_process" = "GO:0002376",
  "response_to_stimulus" = "GO:0050896",
  "signaling" = "GO:0023052",
  "metabolic_process" = "GO:0008152",
  "regulation_of_vasculature_development" = "GO:1901342",
  "signal_transduction" = "GO:0007165")

getGOdatabase <- function(
    parent = NULL,
    root = "BP",
    spe = getOption("spe")) {

  roots_go <- c("BP" = "GO:0008150","MF" = "GO:0003674","CC" = "GO:0005575")
  if(!any(root %in% names(roots_go))) {
    stop("'root' database should be in 'BP', 'MF' and 'CC'")
  }
  GenesetNames <- GetAllChilrenGO(roots_go[root], spe = spe)
  if(is.null(parent)) return(GenesetNames) else {
    if(any(!parent %in% c(GenesetNames, names(DatabaseList)))) {
      ol <- setdiff(parent, c(GenesetNames, names(DatabaseList)))
      stop(length(ol),
           " ‘parent' GO term(s) not found in the current root database (",
           paste0(root, collapse = ", "),
           "): ", ol[1], "...")
    }
  }
  if(any(parent %in% names(DatabaseList)) & !"BP" %in% root) {
    stop("Term(s) '", paste0(intersect(parent, names(DatabaseList)), collapse = "', '"),
         "' are under 'BP' root database. Please add 'BP' to 'root' parameter")
  }
  if(any(parent %in% names(DatabaseList)) ) {
    parent <- unique(plyr::revalue(parent,DatabaseList,warn_missing = F))
  }
  GenesetNames <- GetAllChilrenGO(parent, spe = spe)
  return(GenesetNames)
}
