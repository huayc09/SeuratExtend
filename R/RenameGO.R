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
#' @rdname RenameGO
#' @export

RenameGO <- function(term, add_id = T, add_n_gene = T, spe = getOption("spe")){
  if(is.vector(term)){
    check_spe(spe)
    renamed <- GO_Data[[spe]]$GO_ontology$name[term]
    if(add_id) renamed <- paste(term, renamed)
    if(add_n_gene) renamed <- paste0(renamed, " (", sapply(GO_Data[[spe]]$GO2Gene[term], length), "g)")
    return(renamed)
  }else{
    library(magrittr)
    library(dplyr)
    term <- term %>% set_rownames(RenameGO(rownames(.)))
    return(term)
  }
}

# term = GO_Data[[spe]]$GO_ontology$name[1:10] %>% names
# RenameGO(term)
