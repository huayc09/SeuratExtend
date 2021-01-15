#' @title FUNCTION_TITLE
#' @description FUNCTION_DESCRIPTION
#' @param human_genes PARAM_DESCRIPTION
#' @param mirror PARAM_DESCRIPTION, Default: NULL
#' @param local.mode PARAM_DESCRIPTION, Default: T
#' @param keep.seq PARAM_DESCRIPTION, Default: F
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

HumanToMouseGenesymbol <- function(human_genes, mirror = NULL, local.mode = T, keep.seq = F){
  if(is.vector(human_genes)){
    if(local.mode){
      genelist <- mouse_human_genesymbols[mouse_human_genesymbols$HGNC.symbol %in% human_genes, 2:1]
    } else {
      human_genes <- unique(human_genes)
      import("biomaRt")
      message("Posible mirrors: 'www', 'uswest', 'useast', 'asia'.")
      human <- useDataset("hsapiens_gene_ensembl", useEnsembl(biomart = "ensembl", mirror = mirror))
      mouse <- useDataset("mmusculus_gene_ensembl", useEnsembl(biomart = "ensembl", mirror = mirror))
      genelist = getLDS(attributes = c("hgnc_symbol"), filters = "hgnc_symbol", values = human_genes, mart = human,
                        attributesL = c("mgi_symbol"), martL = mouse, uniqueRows=T)
    }
    if(nrow(genelist) == 0) message("No homologous genes found")
    if(keep.seq) {
      genelist2 <- vector()
      for (i in human_genes) {
        genelist2[i] <- genelist$MGI.symbol[genelist$HGNC.symbol %in% i][1]
        genelist <- genelist[!genelist$HGNC.symbol %in% i & !genelist$MGI.symbol %in% genelist2[i],]
      }
      return(genelist2)
    }
    return(genelist)
  }else{
    library(magrittr)
    library(dplyr)
    human_genes <- human_genes %>%
      set_rownames(HumanToMouseGenesymbol(rownames(.), local.mode = T, keep.seq = T)) %>%
      .[!is.na(rownames(.)),]
    return(human_genes)
  }
}

#' @title FUNCTION_TITLE
#' @description FUNCTION_DESCRIPTION
#' @param mouse_genes PARAM_DESCRIPTION
#' @param mirror PARAM_DESCRIPTION, Default: NULL
#' @param local.mode PARAM_DESCRIPTION, Default: T
#' @param keep.seq PARAM_DESCRIPTION, Default: F
#' @return OUTPUT_DESCRIPTION
#' @details DETAILS
#' @examples
#' \dontrun{
#' if(interactive()){
#'  #EXAMPLE1
#'  }
#' }
#' @rdname MouseToHumanGenesymbol
#' @export

MouseToHumanGenesymbol <- function(mouse_genes, mirror = NULL, local.mode = T, keep.seq = F){
  if(is.vector(mouse_genes)){
    if(local.mode){
      genelist <- mouse_human_genesymbols[mouse_human_genesymbols$MGI.symbol %in% mouse_genes, ]
    } else {
      mouse_genes <- unique(mouse_genes)
      import("biomaRt")
      message("Posible mirrors: 'www', 'uswest', 'useast', 'asia'.")
      human <- useDataset("hsapiens_gene_ensembl", useEnsembl(biomart = "ensembl", mirror = mirror))
      mouse <- useDataset("mmusculus_gene_ensembl", useEnsembl(biomart = "ensembl", mirror = mirror))
      genelist = getLDS(attributesL = c("hgnc_symbol"), filters = "mgi_symbol", values = mouse_genes, martL = human,
                        attributes = c("mgi_symbol"), mart = mouse, uniqueRows=T)
    }
    if(nrow(genelist) == 0) message("No homologous genes found")
    if(keep.seq) {
      genelist2 <- vector()
      for (i in mouse_genes) {
        genelist2[i] <- genelist$HGNC.symbol[genelist$MGI.symbol %in% i][1]
        genelist <- genelist[!genelist$MGI.symbol %in% i & !genelist$HGNC.symbol %in% genelist2[i],]
      }
      return(genelist2)
    }
    return(genelist)
  }else{
    library(magrittr)
    library(dplyr)
    mouse_genes <- mouse_genes %>%
      set_rownames(MouseToHumanGenesymbol(rownames(.), local.mode = T, keep.seq = T)) %>%
      .[!is.na(rownames(.)),]
    return(mouse_genes)
  }
}

