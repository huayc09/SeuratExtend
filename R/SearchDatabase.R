#' @title Search Pathways in GO/Reactome Database
#' @description Search for pathways in the GO or Reactome database using a gene name, pathway ID (SetID), or pathway name (SetName).
#' @param item A gene name, pathway ID, or pathway name.
#' @param type Types of search criteria. Can be one or a combination of "gene", "SetID", "SetName". Default: c("gene", "SetID", "SetName").
#' @param database The database to search in: "GO", "Reactome", or both. Default: c("GO", "Reactome").
#' @param spe Species, either 'human' or 'mouse'. Default: getOption("spe").
#' @param export.to.data.frame If TRUE, the results will be returned as a data frame, which can be easily exported to csv or xlsx formats. Default: FALSE.
#' @param n.min Minimum number of genes required in the gene set. Default: 1.
#' @param n.max Maximum number of genes allowed in the gene set. Default: Inf.
#' @param only.end.terms If TRUE, only the end-level pathways or gene sets will be returned. Default: FALSE.
#' @param return Specifies the format of the result. If set to "all" (default), the result will include all information (pathway ID, pathway name, genes). If set to "ID", it will only return a vector of pathway IDs. If set to "genelist", it will return a list where the names are pathway IDs and the contents are the corresponding genes.
#' @return Depending on the specified parameters, the function returns:
#' \itemize{
#'  \item A full list (with pathway ID, pathway name, and genes) if \code{return = "all"}.
#'  \item A vector of pathway IDs if \code{return = "ID"}.
#'  \item A named list of genes with pathway IDs as the list names if \code{return = "genelist"}.
#'  \item A data frame if \code{export.to.data.frame = TRUE}.
#' }
#' @details Further details and usage can be seen in the provided examples.
#' @examples
#' library(SeuratExtend)
#' library(dplyr)
#' options(max.print = 10)
#'
#' # Set the default database search to 'human'
#' options(spe = "human")
#'
#' # The 'item' parameter can be anything like a gene, pathway ID, or pathway name.
#' # For example, search in the GO/Reactome database for pathways that contain the "CD3D" gene
#' # or pathway names that include "metabolic".
#' result <- SearchDatabase(c("CD3D","metabolic"))
#' names(result)
#' glimpse(head(result, 5))
#'
#' # If you only want to search based on gene names, you can specify this in the 'type' parameter.
#' result <- SearchDatabase("CD3D", type = "gene")
#' names(result)
#'
#' # If you want to search in a specific database, you can specify it using the 'database' parameter.
#' result <- SearchDatabase("CD3D", database = "Reactome")
#' names(result)
#'
#' # You can specify the database to be either 'human' or 'mouse' using the 'spe' parameter.
#' result <- SearchDatabase("Cd3d", spe = "mouse")
#' glimpse(head(result, 5))
#'
#' # If you only want the result to be a list of pathway IDs, which could be useful for subsequent analysis,
#' # you can set the 'return' parameter accordingly.
#' result <- SearchDatabase("CD3D", return = "ID")
#' result
#'
#' # If you want the output as a gene list, with the IDs as names, suitable for inputs like GeneSetAnalysis,
#' # you can adjust the 'return' parameter.
#' result <- SearchDatabase("CD3D", return = "genelist")
#' glimpse(head(result, 5))
#'
#' # If you prefer the output to be a data frame, which is suitable for exporting to Excel or CSV formats,
#' # you can set 'export.to.data.frame' to TRUE.
#' result <- SearchDatabase("CD3D", export.to.data.frame = TRUE)
#' glimpse(result)
#' @rdname SearchDatabase
#' @export

SearchDatabase <- function(
    item,
    type = c("gene", "SetID", "SetName"),
    database = c("GO", "Reactome"),
    spe = getOption("spe"),
    export.to.data.frame = F,
    n.min = 1,
    n.max = Inf,
    only.end.terms = F,
    return = c("all","ID","genelist")
) {
  check_spe(spe)
  ID_GO <- c()
  ID_Reactome <- c()

  if("GO" %in% database) {
    message("Search in GO database")
    ID_GO <- FilterGOTerms(spe = spe, n.min = n.min, n.max = n.max, only.end.terms = only.end.terms)
    filter <- rep(FALSE, length(ID_GO))
    if("gene" %in% type) {
      filter <- filter | sapply(GO_Data[[spe]]$GO2Gene[ID_GO], function(x) any(x %in% item))
    }
    if("SetID" %in% type) {
      filter <- filter | apply(sapply(item, function(x) grepl(x, ID_GO)), 1, any)
    }
    if("SetName" %in% type) {
      name_GO <- RenameGO(term = ID_GO, add_id = F, add_n_gene = F, spe = spe)
      filter <- filter | apply(sapply(item, function(x) grepl(tolower(x), tolower(name_GO))), 1, any)
    }
    ID_GO <- ID_GO[filter]
  }
  if("Reactome" %in% database) {
    message("Search in Reactome database")
    ID_Reactome <- FilterReactomeTerms(spe = spe, n.min = n.min, n.max = n.max, only.end.terms = only.end.terms)
    filter <- rep(FALSE, length(ID_Reactome))
    if("gene" %in% type) {
      filter <- filter | sapply(Reactome_Data[[spe]]$Path2Gene[ID_Reactome], function(x) any(x %in% item))
    }
    if("SetID" %in% type) {
      filter <- filter | apply(sapply(item, function(x) grepl(x, ID_Reactome)), 1, any)
    }
    if("SetName" %in% type) {
      name_Reactome <- RenameReactome(term = ID_Reactome, add_id = F, add_n_gene = F, spe = spe)
      filter <- filter | apply(sapply(item, function(x) grepl(tolower(x), tolower(name_Reactome))), 1, any)
    }
    ID_Reactome <- ID_Reactome[filter]
  }
  return <- return[1]
  if(return == "all" | export.to.data.frame) {
    SearchResult <- data.frame()
    if(!is_empty(ID_GO)){
      tmp <- data.frame(
        SetID = ID_GO,
        SetName = GO_Data[[spe]]$GO_ontology$name[ID_GO]
      )
      tmp$Genes <- GO_Data[[spe]]$GO2Gene[ID_GO]
      rownames(tmp) <- tmp$SetID
      tmp <- RenameGO(tmp, spe = spe)
      SearchResult <- rbind(SearchResult, tmp)
    }
    if(!is_empty(ID_Reactome)){
      tmp <- data.frame(
        SetID = ID_Reactome,
        SetName = Reactome_Data[[spe]]$Ontology$name[ID_Reactome]
      )
      tmp$Genes <- Reactome_Data[[spe]]$Path2Gene[ID_Reactome]
      rownames(tmp) <- tmp$SetID
      tmp <- RenameReactome(tmp, spe = spe)
      SearchResult <- rbind(SearchResult, tmp)
    }
    if(export.to.data.frame) {
      SearchResult$Genes <- sapply(SearchResult$Genes, function(x) paste0(x, collapse = ","))
      return(SearchResult)
    } else {
      SearchResult <- apply(SearchResult, 1, function(x) x)
      return(SearchResult)
    }
  } else if(return == "ID") {
    return(c(ID_GO, ID_Reactome))
  } else if(return == "genelist") {
    return(c(
      GO_Data[[spe]]$GO2Gene[ID_GO],
      Reactome_Data[[spe]]$Path2Gene[ID_Reactome]
    ))
  } else {
    stop("'return' should be either 'all', 'ID', or 'genelist'")
  }
}

#' @title Search Pathways in a Given Geneset List
#' @description This function filters a provided geneset list based on a specified gene name, or pathway name.
#' @param genesets A list of genesets where the names are pathway names and the contents are genes.
#' @param item A vector containing gene names or (parts of) pathway names to search for.
#' @param type Specifies the types of search criteria. This can be "gene", "SetName", or both. Default is c("gene", "SetName").
#' @return A subset of the provided geneset list containing the specified items.
#' @details For more details and use-cases, see the examples section.
#' @examples
#' library(SeuratExtend)
#' library(dplyr)
#' options(max.print = 10)
#'
#' # Check within the hallmark 50 database for pathways that contain the "CD3D" gene
#' # or have pathway names that include "interferon".
#' SearchPathways(genesets = hall50$human, item = c("CD3D", "interferon"))
#' @rdname SearchPathways
#' @export
SearchPathways <- function(
    genesets,
    item,
    type = c("gene", "SetName")
){
  library(dplyr)
  filter <- rep(FALSE, length(genesets))
  if("SetName" %in% type) {
    filter <- filter | (sapply(item, function(x) grepl(tolower(x), tolower(names(genesets)))) %>% apply(1, any))
  }
  if("gene" %in% type) {
    filter <- filter | sapply(genesets, function(x) any(x %in% item))
  }
  return(genesets[filter])
}


