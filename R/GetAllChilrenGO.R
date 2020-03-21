GetChilrenGO <- function(term, spe = getOption("spe")){
  check_spe(spe)
  Child <- union(term, GO_Data[[spe]]$GO_ontology$children[term] %>% unlist)
  return(Child)
}

GetAllChilrenGO <- function(term, spe = getOption("spe")){
  Child <- GetChilrenGO(term, spe)
  if(length(Child) > length(term)) return(GetAllChilrenGO(Child, spe)) else return(term)
}

# GetChilrenGO(term, spe)
# GetAllChilrenGO(term, spe)
