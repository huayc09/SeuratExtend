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

HEV_markers_Girards <- readRDS("~/R documents/SeuratExtend_databases/Other gene sets/HEV_markers_Girards.rds")
Genesets[["mouse"]][["HEV_markers_Girards"]] <- HEV_markers_Girards %>% lapply(rownames)
Genesets[["mouse"]][["HEV_markers_Girards_top30"]] <- HEV_markers_Girards %>% lapply(function(x) rownames(x)[1:30])

EC_3T_markers <- read_excel("~/R documents/SeuratExtend_databases/Other gene sets/EC_3T_markers.xlsx")
EC_markers_Junbin <- as.list(EC_3T_markers)
Genesets[["human"]][["EC_markers_Junbin"]] <- EC_markers_Junbin

Genesets_data <- Genesets
usethis::use_data(Genesets_data, overwrite = TRUE)


