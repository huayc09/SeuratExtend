#' @title FUNCTION_TITLE
#' @description FUNCTION_DESCRIPTION
#' @param Ensembl PARAM_DESCRIPTION
#' @param spe PARAM_DESCRIPTION, Default: getOption("spe")
#' @param mirror PARAM_DESCRIPTION, Default: NULL
#' @return OUTPUT_DESCRIPTION
#' @details DETAILS
#' @examples
#' \dontrun{
#' if(interactive()){
#'  #EXAMPLE1
#'  }
#' }
#' @rdname EnsemblToGenesymbol
#' @export

EnsemblToGenesymbol<-function(Ensembl, spe = getOption("spe"), mirror = NULL){
  library("biomaRt")
  par <- list(mouse = c(Dataset = "mmusculus_gene_ensembl", symbol = "mgi_symbol"),
              human = c(Dataset = "hsapiens_gene_ensembl", symbol = "hgnc_symbol"))
  message("Posible mirrors: 'www', 'uswest', 'useast', 'asia'.")
  mart <- useDataset(par[[spe]]["Dataset"], useEnsembl(biomart = "ensembl", mirror = mirror))
  EnsemblGenes <- getBM(attributes = c(par[[spe]]["symbol"], "ensembl_gene_id"),
                        filters = "ensembl_gene_id",
                        values = Ensembl,
                        mart = mart)
  return(EnsemblGenes)
}

# options(spe = "human")
# par <- list(mouse = c(name = "Mus musculus", title = "MMU"),
#             human = c(name = "Homo sapiens", title = "HSA"))
# Ensembl2Reactome <- Ensembl2Reactome_PE_All_Levels %>% .[.$V8 == par[[spe]]["name"], ]
# Ensembl2Reactome$V1 %>% unique() %>% EnsemblToGenesymbol(mirror = "useast")

