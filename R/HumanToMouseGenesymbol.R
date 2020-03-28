#' @title FUNCTION_TITLE
#' @description FUNCTION_DESCRIPTION
#' @param human_genes PARAM_DESCRIPTION
#' @param mirror PARAM_DESCRIPTION, Default: NULL
#' @return OUTPUT_DESCRIPTION
#' @details DETAILS
#' @examples
#' \dontrun{
#' if(interactive()){
#'  #EXAMPLE1
#'  }
#' }
#' @rdname HumanToMouseGenesymbol
#' @export

HumanToMouseGenesymbol <- function(human_genes, mirror = NULL){
  human_genes <- unique(human_genes)
  library("biomaRt")
  message("Posible mirrors: 'www', 'uswest', 'useast', 'asia'.")
  human <- useDataset("hsapiens_gene_ensembl", useEnsembl(biomart = "ensembl", mirror = mirror))
  mouse <- useDataset("mmusculus_gene_ensembl", useEnsembl(biomart = "ensembl", mirror = mirror))
  genelist = getLDS(attributes = c("hgnc_symbol"), filters = "hgnc_symbol", values = human_genes, mart = human,
                    attributesL = c("mgi_symbol"), martL = mouse, uniqueRows=T)
  return(genelist)
}
