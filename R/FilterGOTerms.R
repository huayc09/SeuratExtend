#' @title FUNCTION_TITLE
#' @description FUNCTION_DESCRIPTION
#' @param term PARAM_DESCRIPTION
#' @param spe PARAM_DESCRIPTION, Default: 'mouse'
#' @param n.min PARAM_DESCRIPTION, Default: 1
#' @param n.max PARAM_DESCRIPTION, Default: Inf
#' @param only.end.terms PARAM_DESCRIPTION, Default: T
#' @param change.name PARAM_DESCRIPTION, Default: F
#' @param parent PARAM_DESCRIPTION, Default: NULL
#' @return OUTPUT_DESCRIPTION
#' @details DETAILS
#' @examples
#' \dontrun{
#' if(interactive()){
#'  #EXAMPLE1
#'  }
#' }
#' @rdname FilterGOTerms
#' @export

FilterGOTerms <- function(term, spe = getOption("spe"), n.min = 1, n.max = Inf,
                          only.end.terms = T, change.name = F, parent = NULL){
  check_spe(spe)
  if(is.vector(term)){
    filter <-
      GO_Data[[spe]]$GO2Gene[term] %>%
      sapply(function(x) length(x) >= n.min & length(x) <= n.max)
    term <- term[filter]
    if(!is.null(parent)){
      if(all(parent %in% GO_Data[[spe]]$GO_ontology$id)) {
        term <- intersect(term, GetAllChilrenGO(parent, spe))
      }else{
        message <- paste0("Parent item(s) not found: ", paste0(setdiff(parent, GO_Data[[spe]]$GO_ontology$id), collapse = ", "))
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

# term = GO_Data[[spe]]$GO_ontology$ancestors[c("GO:0000083", "GO:0000098")] %>% unlist()
# term = seu@misc$AUCell$GO$All
# spe = "mouse"
# n.min = 10
# n.max = 500
# change.name = T
# only.end.terms = T
# FilterGOTerms(term, change.name = F)[1:5,1:5]
# parent = "GO:0008150"
#
# FilterGOTerms(term, parent = "GO:0000083", change.name = T)[,1:5]

