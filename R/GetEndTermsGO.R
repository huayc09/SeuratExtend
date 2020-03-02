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
#' @rdname GetEndTermsGO
#' @export
GetEndTermsGO <- function(term, spe){
  anc <- GO_Data[[spe]]$GO_ontology$ancestors[term] %>%
    lapply(function(x) x[-length(x)]) %>%
    unlist() %>% unique()
  return(setdiff(term, anc))
}

# term = GO_Data[[spe]]$GO_ontology$ancestors[c("GO:0000083", "GO:0000098")] %>% unlist()
# spe = "mouse"
# GetEndTermsGO(term, spe)
