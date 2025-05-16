#' @include GeneSetAnalysis.R
#'
NULL

#' @param parent ID or name of the parent (top-level) gene set in the GO/Reactome database.
#' This restricts the analysis to a subset of gene sets. Default: NULL.
#' @param dataset (GO) Alias for 'parent'. Default: NULL.
#' @param root (GO) Specifies which root category to use:
#' 'BP' (Biological Process), 'MF' (Molecular Function), or 'CC' (Cellular Component). Default: 'BP'.
#' @param spe Species, either 'human' or 'mouse'. Default: getOption("spe").
#' @param only.end.terms Boolean indicating whether only the end-level pathways or gene sets
#' should be used for calculation. Default: FALSE.
#' @rdname GeneSetAnalysis
#' @export

GeneSetAnalysisGO <- function(
    seu = NULL,
    dataset = NULL,
    parent = NULL,
    root = "BP",
    spe = getOption("spe"),
    ratio = 0.4,
    n.min = 1,
    n.max = Inf,
    only.end.terms = F,
    slot = "counts",
    assay = "RNA",
    nCores = 1,
    aucMaxRank = NULL,
    title = NULL,
    export_to_matrix = F,
    verbose = TRUE,
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

  GenesetList <- GO_Data[[spe]]$GO2Gene[GenesetNames]
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

#' @title Filter GO/Reactome Gene Sets
#' @description Subset gene sets based on GO/Reactome IDs.
#' @param term GO/Reactome IDs, either as a character vector or matrix (with IDs in rownames). If NULL, the function will filter from the entire database.
#' @param spe Species, either 'human' or 'mouse'. Default: getOption("spe").
#' @param n.min Minimum number of genes in the gene set. Default: 1.
#' @param n.max Maximum number of genes in the gene set. Default: Inf.
#' @param only.end.terms Boolean indicating whether only return the end-level pathways or gene sets. Default: FALSE.
#' @param change.name Boolean indicating whether to change the gene set names. Default: FALSE.
#' @param parent ID or name of the parent (top-level) gene set in the GO/Reactome database. Default: NULL.
#' @return Returns a character vector or matrix, depending on the input.
#' @examples
#' library(SeuratExtend)
#' options(max.print = 10, spe = "human")
#'
#' # Filter GO terms. For instance, select pathways under the category GO:0002376 (immune system process).
#' terms <- FilterGOTerms(parent = "GO:0002376")
#' RenameGO(terms)
#'
#' # If you want to further restrict the number of genes in a pathway to be between 10 and 1000,
#' # you can adjust the n.min and n.max parameters.
#' terms2 <- FilterGOTerms(term = terms, n.min = 10, n.max = 1000)
#' RenameGO(terms2)
#'
#' # To only return the end-level pathways, set the only.end.terms parameter to TRUE.
#' terms3 <- FilterGOTerms(term = terms, only.end.terms = TRUE)
#' RenameGO(terms3)
#'
#' # Filter Reactome terms. For example, select pathways under the category R-HSA-168256 (Immune System).
#' terms <- FilterReactomeTerms(parent = "R-HSA-168256")
#' RenameReactome(terms)
#' @rdname FilterGOTerms
#' @export

FilterGOTerms <- function(
    term = NULL,
    spe = getOption("spe"),
    n.min = 1,
    n.max = Inf,
    only.end.terms = F,
    change.name = F,
    parent = NULL
){
  check_spe(spe)
  library(SeuratExtendData)
  library(dplyr)
  if(is.null(term)) term <- names(GO_Data[[spe]]$GO2Gene)
  if(is.vector(term)){
    filter <-
      GO_Data[[spe]]$GO2Gene[term] %>%
      sapply(function(x) length(x) >= n.min & length(x) <= n.max)
    term <- term[filter]
    if(!is.null(parent)){
      if(all(parent %in% GO_Data[[spe]]$GO_ontology$id)) {
        term <- intersect(term, GetAllChilrenGO(parent, spe))
      }else{
        message <- paste0(
          "Parent item(s) not found: ",
          paste0(setdiff(parent, GO_Data[[spe]]$GO_ontology$id),
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

#' @title Convert GO/Reactome Pathway IDs to Full Names
#' @description Converts GO or Reactome IDs to their respective full geneset or pathway names.
#' @param term GO/Reactome IDs, which can be either a character vector or a matrix (with IDs in rownames).
#' @param add_id Logical: Should the IDs be displayed? Default: TRUE.
#' @param add_n_gene Logical: Should the number of genes in each geneset be displayed? Default: TRUE.
#' @param spe Species, either 'human' or 'mouse'. Default: getOption("spe").
#' @return Returns either a character vector or a matrix, depending on the input.
#' @details For more examples, see \code{\link[SeuratExtend:GeneSetAnalysis]{GeneSetAnalysis()}}.
#' @examples
#' library(SeuratExtend)
#'
#' RenameGO(c("GO:0002376","GO:0050896"), spe = "human")
#'
#' RenameReactome(c("R-HSA-109582","R-HSA-112316"), spe = "human")
#' @rdname RenameGO
#' @export

RenameGO <- function(
    term,
    add_id = T,
    add_n_gene = T,
    spe = getOption("spe")
) {
  if(is.vector(term)){
    check_spe(spe)
    library(SeuratExtendData)
    filter <- (term %in% GO_Data[[spe]]$GO_ontology$id)
    term2 <- term[filter]
    renamed <- GO_Data[[spe]]$GO_ontology$name[term2]
    if(add_id) renamed <- paste(term2, renamed)
    if(add_n_gene) renamed <- paste0(
      renamed, " (",
      sapply(GO_Data[[spe]]$GO2Gene[term2], length), "g)")
    term[filter] <- renamed
    return(term)
  }else{
    rownames(term) <- RenameGO(rownames(term), add_id = add_id, add_n_gene = add_n_gene, spe = spe)
    return(term)
  }
}

# Internal ----------------------------------------------------------------

GetChilrenGO <- function(term, spe = getOption("spe")){
  check_spe(spe)
  Child <- union(term, unlist(GO_Data[[spe]]$GO_ontology$children[term]))
  return(Child)
}

GetAllChilrenGO <- function(term, spe = getOption("spe")){
  check_spe(spe)
  Child <- GetChilrenGO(term, spe)
  if(length(Child) > length(term)) return(GetAllChilrenGO(Child, spe)) else return(term)
}

GetEndTermsGO <- function(term, spe = getOption("spe")){
  library(dplyr)
  check_spe(spe)
  anc <- GO_Data[[spe]]$GO_ontology$ancestors[term] %>%
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

  check_spe(spe)
  roots_go <- c("BP" = "GO:0008150","MF" = "GO:0003674","CC" = "GO:0005575")
  if(!any(root %in% names(roots_go))) {
    stop("'root' database should be in 'BP', 'MF' and 'CC'")
  }
  GenesetNames <- GetAllChilrenGO(roots_go[root], spe = spe)
  if(is.null(parent)) return(GenesetNames) else {
    if(any(!parent %in% c(GenesetNames, names(DatabaseList)))) {
      ol <- setdiff(parent, c(GenesetNames, names(DatabaseList)))
      stop(length(ol),
           " 'parent' GO term(s) not found in the current root database (",
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
