#' @title FUNCTION_TITLE
#' @description FUNCTION_DESCRIPTION
#' @param Reactome_interactions_filtered PARAM_DESCRIPTION
#' @param path PARAM_DESCRIPTION, Default: getwd()
#' @return OUTPUT_DESCRIPTION
#' @details DETAILS
#' @examples
#' \dontrun{
#' if(interactive()){
#'  #EXAMPLE1
#'  }
#' }
#' @rdname CellphoneDB_GenerateCustomDB
#' @export

CellphoneDB_GenerateCustomDB <- function(Reactome_interactions_filtered, path = getwd()) {
  library(dplyr)
  dir.create("input")
  gene_input <-
    data.frame("gene_name" = c(Reactome_interactions_filtered$Genesymbol.1,
                               Reactome_interactions_filtered$Genesymbol.2),
               "uniprot" = c(Reactome_interactions_filtered$Uniprot.1,
                             Reactome_interactions_filtered$Uniprot.2)) %>%
    .[!duplicated(.$uniprot),] %>%
    mutate("hgnc_symbol" = .$gene_name,
           "ensembl" = .$gene_name)
  write.table(gene_input, file = "input/gene_input.csv", sep = ",", row.names = F, quote = F)

  protein_input <-
    data.frame("uniprot" = gene_input$uniprot,
               "protein_name" = gene_input$gene_name)
  write.table(protein_input, file = "input/protein_input.csv", sep = ",", row.names = F, quote = F)

  interaction_input <-
    data.frame("partner_a" = Reactome_interactions_filtered$Uniprot.1,
               "partner_b" = Reactome_interactions_filtered$Uniprot.2,
               "annotation_strategy" = "curated",
               "source" = Reactome_interactions_filtered$Pubmed.references)
  write.table(interaction_input, file = "input/interaction_input.csv", sep = ",", row.names = F, quote = F)
  system(paste0("cellphonedb database generate --user-interactions input/interaction_input.csv ",
                "--user-protein input/protein_input.csv --user-gene input/gene_input.csv"))
}


# library(Seurat)
# library(SeuratExtend)
# library(purrr)
#
# setwd("~/R documents/cellphonedb test")
# path = getwd()
# CellphoneDB_GenerateCustomDB(Reactome_interactions_filtered = Reactome_interactions_filtered)
