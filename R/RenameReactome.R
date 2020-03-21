#' @title FUNCTION_TITLE
#' @description FUNCTION_DESCRIPTION
#' @param term PARAM_DESCRIPTION
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

RenameReactome <- function(term, spe = getOption("spe")){
  check_spe(spe)
  if(is.vector(term)){
    return(Reactome_Data[[spe]]$Ontology$name[term])
  }else{
    library(magrittr)
    library(dplyr)
    term <- term %>% set_rownames(RenameReactome(rownames(.)))
    return(term)
  }
}

# term = names(Reactome_Data[[spe]]$Path2Gene)
# RenameReactome(term)
# getOption("spe")
# options(spe = "mouse")
