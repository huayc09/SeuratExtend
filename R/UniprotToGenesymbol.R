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
  return(uniprot_db[[spe]][Uniprot,"Gene symbol"])
}
