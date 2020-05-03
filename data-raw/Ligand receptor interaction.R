# UniProt database

setwd("~/R documents/SeuratExtend_databases/2020-4-30 UniProt")
library(readxl)
library(tidyr)
options(max.print = 50)

uniprot_db <- list()
uniprot_db[["mouse"]] <-
  read_excel("uniprot-organism__Mus+musculus+(Mouse)+[10090]_.xlsx") %>%
  as.data.frame()
uniprot_db[["human"]] <-
  read_excel("uniprot-filtered-organism__Homo+sapiens+(Human)+[9606]_.xlsx") %>%
  as.data.frame()
uniprot_db[["mouse"]]$"Gene symbol" <-
  uniprot_db[["mouse"]]$`Gene names` %>%
  strsplit(" ") %>%
  sapply(function(x) x[1])
uniprot_db[["human"]]$"Gene symbol" <-
  uniprot_db[["human"]]$`Gene names` %>%
  strsplit(" ") %>%
  sapply(function(x) x[1])
rownames(uniprot_db[["mouse"]]) <- uniprot_db[["mouse"]]$Entry
rownames(uniprot_db[["human"]]) <- uniprot_db[["human"]]$Entry

setwd("~/R documents/SeuratExtend")
usethis::use_data(uniprot_db, overwrite = TRUE)

# Reactome ligand-receptor database
setwd("~/R documents/SeuratExtend_databases/2020-4-29 Reactome interaction")

reactome.all_species.interactions <-
  read.csv("reactome.all_species.interactions.tab-delimited.txt",
           sep = "\t", header = T)
Reactome_interactions <- list()
Reactome_interactions[["mouse"]] <-
  reactome.all_species.interactions %>%
  .[grepl("MMU", .$Interaction.context),]
Reactome_interactions[["human"]] <-
  reactome.all_species.interactions %>%
  .[grepl("HSA", .$Interaction.context),]

Reactome_interactions[["mouse"]]$Genesymbol.1 <-
  sub("uniprotkb:", "", Reactome_interactions[["mouse"]]$X..Interactor.1.uniprot.id) %>%
  UniprotToGenesymbol()
Reactome_interactions[["mouse"]]$Genesymbol.2 <-
  sub("uniprotkb:", "", Reactome_interactions[["mouse"]]$Interactor.2.uniprot.id) %>%
  UniprotToGenesymbol()
Reactome_interactions[["human"]]$Genesymbol.1 <-
  sub("uniprotkb:", "", Reactome_interactions[["human"]]$X..Interactor.1.uniprot.id) %>%
  UniprotToGenesymbol(spe = "human")
Reactome_interactions[["human"]]$Genesymbol.2 <-
  sub("uniprotkb:", "", Reactome_interactions[["human"]]$Interactor.2.uniprot.id) %>%
  UniprotToGenesymbol(spe = "human")

setwd("~/R documents/SeuratExtend")
usethis::use_data(Reactome_interactions, overwrite = TRUE)

All_genes <-
  unique(c(Reactome_interactions[["mouse"]]$Genesymbol.1,
           Reactome_interactions[["mouse"]]$Genesymbol.2))
filtered_genes <-
  c("Ccl", "Ccr", "Cxc", "Xc", "Cx3c", "Sel", "cam", "Esam") %>%
  sapply(function(x){
    All_genes %>% grepl(x, .)
  }) %>%
  apply(1, any) %>%
  All_genes[.] %>%
  .[!. %in% c("Selenoi", "Ticam1", "Ticam2", "Fcamr", "Sel1l")]

# Reactome_interactions_filtered <-
Reactome_interactions_filtered <-
  Reactome_interactions[["mouse"]][,c("Genesymbol.1","Genesymbol.2")] %>%
  apply(2, function(x) x %in% filtered_genes) %>%
  apply(1, any) %>%
  Reactome_interactions[["mouse"]][.,] %>%
  drop_na() %>%
  mutate("pairs" = .[,c("Genesymbol.1","Genesymbol.2")] %>%
           apply(1, function(x) paste(sort(x), collapse = "_"))) %>%
  .[!duplicated(.$pairs),] %>%
  .[.$Interaction.type == "physical association", ] %>%
  mutate("Uniprot.1" = sub("uniprotkb:","", .$"X..Interactor.1.uniprot.id"),
         "Uniprot.2" = sub("uniprotkb:","", .$"Interactor.2.uniprot.id"),
         "pairs" = paste(.$Genesymbol.1, .$Genesymbol.2, sep = "_"))
setwd("~/R documents/SeuratExtend")
usethis::use_data(Reactome_interactions_filtered, overwrite = TRUE)

# setwd("~/R documents/cellphonedb test")
# saveRDS(Reactome_interactions_filtered, "Reactome_interactions_filtered.rds")


