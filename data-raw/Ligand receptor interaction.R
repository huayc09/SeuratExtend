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

# Reactome_interactions_filtered

library(MeSH.PCR.db)
library(MeSH.db)
library(MeSH.Mmu.eg.db)
library(biomaRt)
library(openxlsx)
setwd("~/MeSH 2020")
cls <- columns(MeSH.PCR.db)
# D015815: Cell Adhesion Molecules
# D018925: Chemokines
# children of these 2 categories
res <- select(MeSH.PCR.db, keys=c("D015815","D018925"), columns="CHILD", keytype="PARENT")[[1]]
# filter out non-gene terms (category)
res <- setdiff(res, select(MeSH.PCR.db, keys=res, columns="PARENT", keytype="CHILD")[[1]])
# get mesh term
res <- select(MeSH.db, keys=res, columns=c("MESHID","MESHTERM"), keytype="MESHID")
# match entrez id
res2 <- select(MeSH.Mmu.eg.db, keys=res$MESHID,
               columns=columns(MeSH.Mmu.eg.db), keytype="MESHID")
res2 <- res2[res2$SOURCEDB=="gene2pubmed",]
res2 <- sapply(res$MESHID, function(x) names(which.max(table(res2[res2$MESHID==x,"GENEID"]))))
res <- data.frame(res, "EntrezID" = res2)
# match entrez id to gene symbol
mart <- useDataset("mmusculus_gene_ensembl", useEnsembl(biomart = "ensembl"))
Genes <- getBM(attributes = c("entrezgene_id","mgi_symbol"),
               filters = "entrezgene_id",
               values = res$EntrezID,
               mart = mart)
rownames(Genes) <- Genes$entrezgene_id
res$Symbol <- Genes[res$EntrezID, "mgi_symbol"]
write.xlsx(res, "match.xlsx")
# manually correct mistake matches
match <- read_excel("match.xlsx")
match$Symbol[!is.na(match$Corrected)] <- match$Corrected[!is.na(match$Corrected)]
match <- data.frame(match[,c("MESHID","MESHTERM","Symbol")])
filtered_genes <- match$Symbol

library(dplyr)
library(SeuratExtend)
Reactome_interact_mouse_chemo_adh <-
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
         "Uniprot.2" = sub("uniprotkb:","", .$"Interactor.2.uniprot.id")) %>%
  .[,c("Genesymbol.1","Genesymbol.2","Uniprot.1","Uniprot.2","Pubmed.references")] %>%
  apply(1, function(x) {
    if(!x[["Genesymbol.1"]] %in% filtered_genes) x <- x[c(2,1,4,3,5)]
    return(x)
  }) %>% t() %>% as.data.frame(stringsAsFactors = F) %>%
  `colnames<-`(c("Genesymbol.1","Genesymbol.2","Uniprot.1","Uniprot.2","Pubmed.references"))

setwd("~/R documents/SeuratExtend")
usethis::use_data(Reactome_interact_mouse_chemo_adh, overwrite = TRUE)

# setwd("~/R documents/cellphonedb test")
# saveRDS(Reactome_interactions_filtered, "Reactome_interactions_filtered.rds")


