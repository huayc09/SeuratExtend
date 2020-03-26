#' @title FUNCTION_TITLE
#' @description FUNCTION_DESCRIPTION
#' @param term PARAM_DESCRIPTION
#' @param add_id PARAM_DESCRIPTION, Default: T
#' @param add_n_gene PARAM_DESCRIPTION, Default: T
#' @param spe PARAM_DESCRIPTION, Default: getOption("spe")
#' @return OUTPUT_DESCRIPTION
#' @details DETAILS
#' @examples
#' \dontrun{
#' if(interactive()){
#'  #EXAMPLE1
#'  }
#' }
#' @rdname RenameReactome
#' @export

RenameReactome <- function(term, add_id = T, add_n_gene = T, spe = getOption("spe")){
  check_spe(spe)
  if(is.vector(term)){
    check_spe(spe)
    renamed <- Reactome_Data[[spe]]$Ontology$name[term]
    if(add_id) renamed <- paste(term, renamed)
    if(add_n_gene) renamed <- paste0(renamed, " (", sapply(Reactome_Data[[spe]]$Path2Gene[term], length), "g)")
    return(renamed)
  }else{
    library(magrittr)
    library(dplyr)
    term <- term %>% set_rownames(RenameReactome(rownames(.)))
    return(term)
  }
}

# term = names(Reactome_Data[[spe]]$Path2Gene)[1:10]
# RenameReactome(term)
# getOption("spe")
# options(spe = "mouse")
