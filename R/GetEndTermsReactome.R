#' @title FUNCTION_TITLE
#' @description FUNCTION_DESCRIPTION
#' @param term PARAM_DESCRIPTION
#' @param spe PARAM_DESCRIPTION
#' @return OUTPUT_DESCRIPTION
#' @details DETAILS
#' @examples
#' \dontrun{
#' if(interactive()){
#'  #EXAMPLE1
#'  }
#' }
#' @rdname GetEndTermsReactome
#' @export

GetEndTermsReactome <- function(term, spe = getOption("spe")){
  check_spe(spe)
  anc <- Reactome_Data[[spe]]$Ontology$ancestors[term] %>%
    lapply(function(x) x[-length(x)]) %>%
    unlist() %>% unique()
  return(setdiff(term, anc))
}

# term <- names(Reactome_Data[[spe]]$Path2Gene)
# term
# GetEndTermsReactome(term)
