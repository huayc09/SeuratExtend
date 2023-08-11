#' @include GeneSetAnalysis.R
#'
NULL

#' @rdname GeneSetAnalysis
#' @export

GeneSetAnalysisReactome <- function(
    seu = NULL,
    parent = "All", spe = getOption("spe"),
    ratio = 0.4, n.min = 1, n.max = Inf, only.end.terms = F,
    slot = "counts", assay = "RNA", nCores = 1,
    aucMaxRank = NULL,
    title = NULL,
    export_to_matrix = F, verbose = TRUE,
    n.items.part = NULL) {

  check_spe(spe)
  library(SeuratObject)
  library(SeuratExtendData)

  if(is.null(seu)){
    rd <- SeuratExtendData::Reactome_Data[[spe]]
    DatabaseList <- rd$Roots
    message("Commonly used datasets:\n  ",
            paste(format(names(DatabaseList)), DatabaseList, sep = ": " ,collapse = "\n  "))
    return(DatabaseList)
  }

  if(verbose) message(Sys.time(), " Retrieve Reactome database")
  GenesetNames <- getReactomedatabase(parent = parent, spe = spe)
  if(is.null(title)) {
    if(is.null(parent)) title <- "All" else title <- parent[1]
  }

  GenesetList <- SeuratExtendData::Reactome_Data[[spe]]$Path2Gene[GenesetNames]
  DefaultAssay(seu) <- assay
  GenesetList <- FilterGenesets(
    genes_in_data = rownames(seu),
    genesets = GenesetList,
    ratio = ratio,
    n.min = n.min,
    n.max = n.max,
    verbose = verbose)

  if(only.end.terms) {
    GenesetNames <- GetEndTermsReactome(names(GenesetList), spe)
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
  seu@misc[["AUCell"]][["Reactome"]][[title]] <- AUC_matrix
  return(seu)
}

#' @rdname FilterGOTerms
#' @export

FilterReactomeTerms <-
  function(term, spe = getOption("spe"), n.min = 1, n.max = Inf,
           only.end.terms = F, change.name = F, parent = NULL){
    check_spe(spe)
    library(SeuratExtendData)
    if(is.vector(term)){
      filter <-
        SeuratExtendData::Reactome_Data[[spe]]$Path2Gene[term] %>%
        sapply(function(x) length(x) >= n.min & length(x) <= n.max)
      term <- term[filter]
      if(!is.null(parent)){
        if(all(parent %in% SeuratExtendData::Reactome_Data[[spe]]$Ontology$id)) {
          term <- intersect(term, GetAllChilrenReactome(parent, spe))
        }else{
          message <- paste0(
            "Parent item(s) not found: ",
            paste0(setdiff(parent, SeuratExtendData::Reactome_Data[[spe]]$Ontology$id), collapse = ", "))
          warning(message)
        }
      }
      if(only.end.terms) term <- GetEndTermsReactome(term, spe)
      if(length(term) == 0) stop("No item(s) passed the filter")
    }else{
      term <- term[FilterReactomeTerms(rownames(term), spe = spe, n.min = n.min, n.max = n.max,
                                       only.end.terms = only.end.terms, change.name = F, parent = parent),]
    }
    if(change.name) term <- RenameReactome(term)
    return(term)
  }

#' @rdname RenameGO
#' @export

RenameReactome <- function(term, add_id = T, add_n_gene = T, spe = getOption("spe")){
  check_spe(spe)
  library(SeuratExtendData)
  if(is.vector(term)){
    check_spe(spe)
    renamed <- SeuratExtendData::Reactome_Data[[spe]]$Ontology$name[term]
    if(add_id) renamed <- paste(term, renamed)
    if(add_n_gene) renamed <- paste0(
      renamed, " (",
      sapply(SeuratExtendData::Reactome_Data[[spe]]$Path2Gene[term], length), "g)")
    return(renamed)
  }else{
    library(magrittr)
    library(dplyr)
    term <- term %>% set_rownames(RenameReactome(rownames(.), add_id = add_id, add_n_gene = add_n_gene))
    return(term)
  }
}

# Internal ----------------------------------------------------------------

GetChilrenReactome <- function(term, spe = getOption("spe")){
  library(dplyr)
  check_spe(spe)
  Child <- union(term, SeuratExtendData::Reactome_Data[[spe]]$Ontology$children[term] %>% unlist)
  return(Child)
}

GetAllChilrenReactome <- function(term, spe = getOption("spe")){
  check_spe(spe)
  Child <- GetChilrenReactome(term, spe)
  if(length(Child) > length(term)) return(GetAllChilrenReactome(Child, spe)) else return(term)
}

GetEndTermsReactome <- function(term, spe = getOption("spe")){
  library(dplyr)
  check_spe(spe)
  anc <- SeuratExtendData::Reactome_Data[[spe]]$Ontology$ancestors[term] %>%
    lapply(function(x) x[-length(x)]) %>%
    unlist() %>% unique()
  return(setdiff(term, anc))
}

getReactomedatabase <- function(
    parent = NULL,
    spe = getOption("spe")) {

  check_spe(spe)
  rd <- SeuratExtendData::Reactome_Data[[spe]]
  DatabaseList <- rd$Roots
  all_term_names <- na.omit(rd$Ontology$name)
  parent <- parent %||% "All"

  if("All" %in% parent){
    GenesetNames <- names(rd$Path2Gene)
    return(GenesetNames)
  }else if(any(!parent %in% c(names(rd$Path2Gene), all_term_names))) {
    message("Commonly used datasets:\n  ",
            paste(format(names(DatabaseList)), DatabaseList, sep = ": " ,collapse = "\n  "))
    stop("Term(s) '", paste0(head(setdiff(parent, c(names(rd$Path2Gene), all_term_names))), collapse = "', '"),
         "'... not found in Reactome database. Please enter Reactome pathway ID or pathway name. ")
  }
  if(any(parent %in% all_term_names)) {
    parent <- unique(plyr::revalue(
      parent,
      setNames(names(all_term_names), all_term_names),
      warn_missing = F))
  }
  GenesetNames <- GetAllChilrenReactome(parent, spe = spe)
  return(GenesetNames)
}
