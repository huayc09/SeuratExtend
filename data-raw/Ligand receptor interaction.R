# UniProt database

setwd("~/R documents/SeuratExtend_databases/2020-5-5 UniProt")
library(readxl)
library(tidyr)
options(max.print = 50)

uniprot_db <- list()
uniprot_db[["mouse"]] <- read.csv("uniprot-filtered-organism__Mus+musculus+(Mouse)+[10090]_.tab",
                                  sep = "\t", stringsAsFactors = F, check.names = F)
uniprot_db[["human"]] <- read.csv("uniprot-filtered-organism__Homo+sapiens+(Human)+[9606]_.tab",
                                  sep = "\t", stringsAsFactors = F, check.names = F)
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
for (i in c("mouse","human")) {
  Reactome_interactions[[i]]$Uniprot.1 <-
    sub("uniprotkb:", "", Reactome_interactions[[i]]$X..Interactor.1.uniprot.id)
  Reactome_interactions[[i]]$Uniprot.2 <-
    sub("uniprotkb:", "", Reactome_interactions[[i]]$Interactor.2.uniprot.id)
  Reactome_interactions[[i]]$Genesymbol.1 <-
    Reactome_interactions[[i]]$Uniprot.1 %>%
    UniprotToGenesymbol(spe = i)
  Reactome_interactions[[i]]$Genesymbol.2 <-
    Reactome_interactions[[i]]$Uniprot.2 %>%
    UniprotToGenesymbol(spe = i)
}

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

# filter interaction pairs in Reactome
filter_interact_pairs <- function(filtered_genes, spe = getOption("spe")){
  check_spe(spe)
  library(dplyr)
  library(SeuratExtend)
  library(tidyr)
  Reactome_interact_filtered <-
    Reactome_interactions[["mouse"]][,c("Genesymbol.1","Genesymbol.2")] %>%
    apply(2, function(x) x %in% filtered_genes) %>%
    apply(1, any) %>%
    Reactome_interactions[["mouse"]][.,] %>%
    drop_na() %>%
    mutate("pairs" = .[,c("Genesymbol.1","Genesymbol.2")] %>%
             apply(1, function(x) paste(sort(x), collapse = "_"))) %>%
    .[!duplicated(.$pairs),] %>%
    mutate("if.mebrane.1" = .$Genesymbol.1 %>% sapply(check_sub_loc),
           "if.mebrane.2" = .$Genesymbol.2 %>% sapply(check_sub_loc)) %>%
    .[.$if.mebrane.1 & .$if.mebrane.2, ] %>%
    .[.$Interaction.type == "physical association", ]
  return(Reactome_interact_filtered)
}
Reactome_interact_mouse_chemo_adh <- filter_interact_pairs(filtered_genes)

setwd("~/R documents/SeuratExtend")
usethis::use_data(Reactome_interact_mouse_chemo_adh, overwrite = TRUE)

# setwd("~/R documents/cellphonedb test")
# saveRDS(Reactome_interactions_filtered, "Reactome_interactions_filtered.rds")


