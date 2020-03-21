#' @title FUNCTION_TITLE
#' @description FUNCTION_DESCRIPTION
#' @param term PARAM_DESCRIPTION
#' @param spe PARAM_DESCRIPTION, Default: getOption("spe")
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
#' @rdname FilterReactomeTerms
#' @export

FilterReactomeTerms <-
  function(term, spe = getOption("spe"), n.min = 1, n.max = Inf,
           only.end.terms = T, change.name = F, parent = NULL){
    check_spe(spe)
    if(is.vector(term)){
      filter <-
        Reactome_Data[[spe]]$Path2Gene[term] %>%
        sapply(function(x) length(x) >= n.min & length(x) <= n.max)
      term <- term[filter]
      if(!is.null(parent)){
        if(all(parent %in% Reactome_Data[[spe]]$Ontology$id)) {
          term <- intersect(term, GetAllChilrenReactome(parent, spe))
        }else{
          message <- paste0("Parent item(s) not found: ", paste0(setdiff(parent, GO_Data[[spe]]$GO_ontology$id), collapse = ", "))
          warning(message)
        }
      }
      if(only.end.terms) term <- GetEndTermsReactome(term, spe)
      if(length(term) == 0) stop("No item(s) passed the filter")
    }else{
      term <- term[FilterReactomeTerms(rownames(term), spe = spe, n.min = n.min, n.max = n.max,
                                       only.end.terms = only.end.terms, change.name = F, parent = parent),]
    }
    if(change.name) term <- RenameReactome(term)
    return(term)
  }

# term = names(Reactome_Data[[spe]]$Path2Gene)
# term = seu@misc$AUCell
# spe = "mouse"
# n.min = 10
# n.max = 500
# change.name = T
# only.end.terms = T
# FilterGOTerms(term, change.name = F)[1:5,1:5]
# parent = "R-MMU-1430728"
#
# FilterReactomeTerms(term, parent = parent, change.name = F)
