#' @title FUNCTION_TITLE
#' @description FUNCTION_DESCRIPTION
#' @param item PARAM_DESCRIPTION
#' @param type PARAM_DESCRIPTION, Default: c("gene", "SetID", "SetName")
#' @param database PARAM_DESCRIPTION, Default: c("GO", "Reactome")
#' @param spe PARAM_DESCRIPTION, Default: getOption("spe")
#' @param export.to.data.frame PARAM_DESCRIPTION, Default: F
#' @param n.min PARAM_DESCRIPTION, Default: 1
#' @param n.max PARAM_DESCRIPTION, Default: Inf
#' @param only.end.terms PARAM_DESCRIPTION, Default: F
#' @return OUTPUT_DESCRIPTION
#' @details DETAILS
#' @examples
#' \dontrun{
#' if(interactive()){
#'  #EXAMPLE1
#'  }
#' }
#' @rdname SearchDatabase
#' @export

SearchDatabase <-
  function(item, type = c("gene", "SetID", "SetName"), database = c("GO", "Reactome"),
           spe = getOption("spe"), export.to.data.frame = F,
           n.min = 1, n.max = Inf, only.end.terms = F){
    check_spe(spe)
    library(dplyr)
    library(rlist)
    ID_GO <- c()
    ID_Rct <- c()
    if("gene" %in% type){
      message("Check if item(s) are gene symbol")
      if("GO" %in% database) {
        message("Check in GO database")
        Allgenes <- GO_Data[[spe]]$GO2Gene %>% unlist() %>% unique()
        if(any(item %in% Allgenes)){
          message <- paste0(intersect(item, Allgenes), collapse = ", ")
          message <- paste0("Gene symbol(s) found in GO database: " ,message)
          message(message)
          for (i in intersect(item, Allgenes)) {
            filter <- GO_Data[[spe]]$GO2Gene %>% sapply(function(x) i %in% x)
            message <- paste(i, "is found in", sum(filter), "GO sets")
            message(message)
            ID_GO <- union(ID_GO, names(GO_Data[[spe]]$GO2Gene)[filter])
          }
        }else{
          message("Gene symbol(s) not found in GO database")
        }
      }
      if("Reactome" %in% database) {
        message("Check in Reactome database")
        Allgenes <- Reactome_Data[[spe]]$Path2Gene %>% unlist() %>% unique()
        if(any(item %in% Allgenes)){
          message <- paste0(intersect(item, Allgenes), collapse = ", ")
          message <- paste0("Gene symbol(s) found in Reactome database: " ,message)
          message(message)
          for (i in intersect(item, Allgenes)) {
            filter <- Reactome_Data[[spe]]$Path2Gene %>% sapply(function(x) i %in% x)
            message <- paste(i, "is found in", sum(filter), "Reactome sets")
            message(message)
            ID_Rct <- union(ID_Rct, names(Reactome_Data[[spe]]$Path2Gene)[filter])
          }
        }else{
          message("Gene symbol(s) not found in Reactome database")
        }
      }
    }
    if("SetID" %in% type){
      message("\nCheck if item(s) match Set ID")
      if("GO" %in% database) {
        message("Check in GO database")
        AllsetID <- names(GO_Data[[spe]]$GO2Gene)
        for (i in item) {
          filter <- grepl(i, AllsetID)
          if(any(filter)){
            Set_i <- AllsetID[filter]
            message <- paste0("Set ID(s) matched to ",i, ": ", paste0(Set_i, collapse = ", "))
            message(message)
            ID_GO <- union(ID_GO, Set_i)
          }else{
            message <- paste0(i, ": no set IDs matched in GO database")
            message(message)
          }
        }
      }
      if("Reactome" %in% database) {
        message("Check in Reactome database")
        AllsetID <- names(Reactome_Data[[spe]]$Path2Gene)
        for (i in item) {
          filter <- grepl(i, AllsetID)
          if(any(filter)){
            Set_i <- AllsetID[filter]
            message <- paste0("Set ID(s) matched to ",i, ": ", paste0(Set_i, collapse = ", "))
            message(message)
            ID_Rct <- union(ID_Rct, Set_i)
          }else{
            message <- paste0(i, ": no set IDs matched in Reactome database")
            message(message)
          }
        }
      }
    }
    if("SetName" %in% type){
      message("\nCheck if item(s) match Set Names")
      if("GO" %in% database) {
        message("Check in GO database")
        AllsetName <- GO_Data[[spe]]$GO_ontology$name
        for (i in item) {
          filter <- grepl(tolower(i), tolower(AllsetName))
          if(any(filter)){
            name_i <- AllsetName[filter]
            message <- paste0(sum(filter), " set name(s) matched to ",i, ":\n",
                              paste0(name_i, collapse = "\n"))
            message(message)
            ID_GO <- union(ID_GO, names(name_i))
          }else{
            message <- paste0(i, ": no set name(s) matched in GO database")
            message(message)
          }
        }
      }
      if("Reactome" %in% database) {
        message("Check in Reactome database")
        AllsetName <- Reactome_Data[[spe]]$Ontology$name
        for (i in item) {
          filter <- grepl(tolower(i), tolower(AllsetName))
          if(any(filter)){
            name_i <- AllsetName[filter]
            message <- paste0(sum(filter), " set name(s) matched to ",i, ":\n",
                              paste0(name_i, collapse = "\n"))
            message(message)
            ID_Rct <- union(ID_Rct, names(name_i))
          }else{
            message <- paste0(i, ": no set name(s) matched in Reactome database")
            message(message)
          }
        }
      }
    }
    message("\nFiltering outputs")
    SearchResult <- list()
    if("GO" %in% database & !is.null(ID_GO)){
      ID_GO <- FilterGOTerms(ID_GO, spe = spe, n.min = n.min, n.max = n.max,
                             only.end.terms = only.end.terms)
      for (i in ID_GO) {
        SetName <- GO_Data[[spe]]$GO_ontology$name[[i]]
        Genes <- GO_Data[[spe]]$GO2Gene[[i]]
        ListName <- paste0(i, ": ", SetName, " (", length(Genes),"g)")
        SearchResult[[ListName]]$SetID <- i
        SearchResult[[ListName]]$SetName <- SetName
        SearchResult[[ListName]]$Genes <- Genes
      }
    }
    if("Reactome" %in% database & !is.null(ID_Rct)){
      ID_Rct <- FilterReactomeTerms(ID_Rct, spe = spe, n.min = n.min, n.max = n.max,
                                    only.end.terms = only.end.terms)
      for (i in ID_Rct) {
        SetName <- Reactome_Data[[spe]]$Ontology$name[[i]]
        Genes <- Reactome_Data[[spe]]$Path2Gene[[i]]
        ListName <- paste0(i, ": ", SetName, " (", length(Genes),"g)")
        SearchResult[[ListName]]$SetID <- i
        SearchResult[[ListName]]$SetName <- SetName
        SearchResult[[ListName]]$Genes <- Genes
      }
    }
    if(export.to.data.frame){
      SearchResult <- SearchResult %>%
        lapply(function(x){
          data.frame(SetID = x$SetID, SetName = x$SetName,
                     Genes = paste0(x$Genes, collapse = ","))
          }) %>%
        list.rbind()
    }
    return(SearchResult)
  }

# item = c("1234", "vegf", "Selp", "Sele")
# item = "vegf"
# type = c("gene", "SetID", "SetName")
# database = c("GO", "Reactome")
# spe = getOption("spe")
# export.to.data.frame = F
# SearchDatabase(item) %>% View

