setwd("~/R documents/SeuratExtend_databases/2020-3-28 GSEA")

gmtlist <- list.files()
names(gmtlist) <-
  c("positional gene sets",
    "all curated gene sets",
    "chemical and genetic perturbations",
    "BioCarta gene sets",
    "KEGG gene sets",
    "PID gene sets",
    "all canonical pathways",
    "all motif gene sets",
    "transcription factor targets",
    "all computational gene sets",
    "all immunologic signatures gene sets",
    "hallmark gene sets"
  )

library(GSEABase)
library(dplyr)
library(readxl)
Gmt2list <- function(Gmt){
  GenesetList <- Gmt %>%
    lapply(geneIds) %>%
    setNames(names(Gmt))
  return(GenesetList)
}
Genesets <- list()
Genesets[["human"]][["GSEA"]] <-
  gmtlist %>%
  lapply(function(x) getGmt(x) %>% Gmt2list)

transformed_genes <-
  Genesets[["human"]] %>%
  unlist() %>%
  unique() %>%
  HumanToMouseGenesymbol(mirror = "uswest")

Genesets[["mouse"]][["GSEA_mouse_gene_transformed"]] <-
  Genesets[["human"]][["GSEA"]] %>%
  lapply(function(x){
    lapply(x, function(y){
      transformed_genes$MGI.symbol[transformed_genes$HGNC.symbol %in% y] %>%
        unique()
    })
  })

Genesets_data <- Genesets
usethis::use_data(Genesets_data, overwrite = TRUE) # this data is moved to SeuratExtendData package

hall50 <- list()
hall50[["human"]] <- Genesets_data$human$GSEA$`hallmark gene sets`
hall50[["mouse"]] <- Genesets_data$mouse$GSEA_mouse_gene_transformed$`hallmark gene sets`
usethis::use_data(hall50, overwrite = TRUE)
