GetChilrenGO <- function(term, spe){
  Child <- unique(c(term, GO_Data[[spe]]$GO_ontology$children[term] %>% unlist))
  return(Child)
}

GetAllChilrenGO <- function(term, spe){
  Child <- GetChilrenGO(term, spe)
  if(length(Child) > length(term)) return(GetAllChilrenGO(Child, spe)) else return(term)
}

# GetChilrenGO(term, spe)
# GetAllChilrenGO(term, spe)
