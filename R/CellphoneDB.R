#' @title FUNCTION_TITLE
#' @description FUNCTION_DESCRIPTION
#' @param Reactome_interactions_filtered PARAM_DESCRIPTION
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

CellphoneDB_GenerateCustomDB <- function(Reactome_interactions_filtered) {
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

#' @title FUNCTION_TITLE
#' @description FUNCTION_DESCRIPTION
#' @param seu PARAM_DESCRIPTION
#' @param group.by PARAM_DESCRIPTION
#' @param database PARAM_DESCRIPTION, Default: 'cellphonedb_mouse_cytokines'
#' @return OUTPUT_DESCRIPTION
#' @details DETAILS
#' @examples
#' \dontrun{
#' if(interactive()){
#'  #EXAMPLE1
#'  }
#' }
#' @rdname RunCellphoneDB
#' @export

RunCellphoneDB <- function(seu, group.by, database = "cellphonedb_mouse_cytokines") {
  message("Please make sure CellphoneDB is installed and able to run in terminal")
  message("https://github.com/Teichlab/cellphonedb")

  library(Seurat)
  library(dplyr)

  dir.create("input")
  if(is.null(database)) {
    command <- "cellphonedb method statistical_analysis input/meta_data.txt input/counts.txt"
  }else if(database == "cellphonedb_mouse_cytokines"){
    file.copy(system.file("extdata/CellphoneDB", "cellphonedb_mouse_cytokines.db", package = "SeuratExtend"),
              "input/custom_db.db")
    command <- paste0("cellphonedb method statistical_analysis input/meta_data.txt input/counts.txt ",
                      "--database input/custom_db.db")
  }else{
    file.copy(database, "input/custom_db.db")
    command <- paste0("cellphonedb method statistical_analysis input/meta_data.txt input/counts.txt ",
                      "--database input/custom_db.db")
  }

  data <- GetAssayData(seu)
  meta <- seu@meta.data
  data_export <- cbind(data.frame("Gene"=rownames(data)), data)
  meta_export <- data.frame("Cell" = rownames(meta), "cell_type" = meta[,group.by])
  write.table(meta_export, file = "input/meta_data.txt", row.names = F, sep = "\t", quote = F)
  write.table(data_export, file = "input/counts.txt", row.names = F, sep = "\t", quote = F)
  system(command)
}

# library(Seurat)
# library(SeuratExtend)
# library(purrr)
#
# setwd("~/R documents/cellphonedb test")
# path = getwd()
# CellphoneDB_GenerateCustomDB(Reactome_interactions_filtered = Reactome_interactions_filtered)
# seu <- readRDS("~/R documents/2020-2-10 EC PyMT and E0771/rds/PyMTEC_old.rds")
# group.by <- "cluster"
# RunCellphoneDB(seu, group.by)
