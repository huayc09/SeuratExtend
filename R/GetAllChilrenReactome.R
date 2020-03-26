GetChilrenReactome <- function(term, spe = getOption("spe")){
  check_spe(spe)
  Child <- union(term, Reactome_Data[[spe]]$Ontology$children[term] %>% unlist)
  return(Child)
}

GetAllChilrenReactome <- function(term, spe = getOption("spe")){
  Child <- GetChilrenReactome(term, spe)
  if(length(Child) > length(term)) return(GetAllChilrenReactome(Child, spe)) else return(term)
}

# options(spe = "mouse")
# term <- c("R-MMU-1430728", "R-MMU-162582")
# GetChilrenReactome(term)
# GetAllChilrenReactome(term)

