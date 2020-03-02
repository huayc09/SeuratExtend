#' @title FUNCTION_TITLE
#' @description FUNCTION_DESCRIPTION
#' @param term PARAM_DESCRIPTION
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
RenameGO <- function(term){
  if(is.vector(term)){
    return(GO_ontology$name[term])
  }else{
    library(magrittr)
    library(dplyr)
    term <- term %>% set_rownames(RenameGO(rownames(.)))
    return(term)
  }
}
