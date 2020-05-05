#' @title FUNCTION_TITLE
#' @description FUNCTION_DESCRIPTION
#' @param Uniprot PARAM_DESCRIPTION
#' @param spe PARAM_DESCRIPTION, Default: getOption("spe")
#' @return OUTPUT_DESCRIPTION
#' @details DETAILS
#' @examples
#' \dontrun{
#' if(interactive()){
#'  #EXAMPLE1
#'  }
#' }
#' @rdname UniprotToGenesymbol
#' @export

UniprotToGenesymbol <- function(Uniprot, spe = getOption("spe")) {
  check_spe(spe)
  if(any(!Uniprot %in% uniprot_db[[spe]]$Entry)) {
    Unknown <- setdiff(Uniprot, uniprot_db[[spe]]$Entry)
    message(paste0(length(Unknown), " item(s) not found in Uniprot ID list: \n",
                   paste0(head(Unknown, 10), collapse = "\n"),
                   ifelse(length(Unknown)>10, "\n...", "")))
  }
  return(uniprot_db[[spe]][Uniprot,"Gene names  (primary )"])
}

#' @title FUNCTION_TITLE
#' @description FUNCTION_DESCRIPTION
#' @param Symbol PARAM_DESCRIPTION
#' @param spe PARAM_DESCRIPTION, Default: getOption("spe")
#' @param loc PARAM_DESCRIPTION, Default: c("Secreted", "Membrane", "Cell membrane")
#' @return OUTPUT_DESCRIPTION
#' @details DETAILS
#' @examples
#' \dontrun{
#' if(interactive()){
#'  #EXAMPLE1
#'  }
#' }
#' @rdname check_sub_loc
#' @export

check_sub_loc <-
  function(Symbol, spe = getOption("spe"), loc = c("Secreted","Membrane","Cell membrane")) {
  check_spe(spe)
  if(!Symbol %in% uniprot_db[[spe]]$`Gene names  (primary )`) return(NA)
  cc <- uniprot_db[[spe]] %>%
    .[.$"Gene names  (primary )" %in% Symbol, "Subcellular location [CC]"] %>%
    strsplit(split = "[:;{}.]") %>%
    unlist() %>%
    trimws(which = "both")
  return(any(cc %in% loc))
}

# Reactome_interact_mouse_chemo_adh$Genesymbol.1 %>% sapply(check_sub_loc)
# Symbol <- "Trav19"
# filtered_genes %>% sapply(check_sub_loc)
