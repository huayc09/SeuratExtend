options(max.print = 50)
library(mgsa)
library(ontologyIndex)
library(dplyr)
library(openxlsx)
library(readxl)
# original databases are moved to "SeuratExtend_databases"
GO_Annot_Mgi <- readGAF("data-raw/mgi.gaf")
GO_Annot_human <- readGAF("data-raw/goa_human.gaf")
GO_ontology <- get_OBO("data-raw/go-basic.obo")

usethis::use_data(GO_ontology, overwrite = TRUE)

GO_Data <- list()
GO_Annot <- GO_Annot_Mgi
spe <- "mouse"
# GO_Annot <- GO_Annot_human
# spe <- "human"

filtered_sets <- intersect(GO_Annot@sets %>% names, GO_ontology$id[!GO_ontology$obsolete])
check_parent_sets <- function(sets){
  if(any(!sets %in% filtered_sets)){
    sets_new <- c(intersect(sets, filtered_sets), GO_ontology$parents[setdiff(sets, filtered_sets)]) %>%
      unlist() %>%
      unique()
    return(check_parent_sets(sets_new))
  }else{
    return(sets)
  }
}
parents <- GO_ontology$parents[filtered_sets] %>% lapply(check_parent_sets)
GO_ontology_new <- ontology_index(parents, name = GO_ontology$name[filtered_sets])
sets_annot_new <- GO_Annot@sets[filtered_sets]
for (i in filtered_sets) {
  extra_genes <- lapply(sets_annot_new[GO_ontology_new$parents[[i]]],
                        function(x) setdiff(sets_annot_new[[i]], x)) %>% unlist
  if(length(extra_genes)>0){
    for (j in GO_ontology_new$ancestors[[i]]) {
      sets_annot_new[[j]] <- union(sets_annot_new[[j]], sets_annot_new[[i]])
    }
  }
}
sets_annot_new <- lapply(sets_annot_new, function(x) GO_Annot@itemAnnotations$symbol[x] %>% as.character)

GO_Data[[spe]] <- list(GO_ontology = GO_ontology_new,
                       GO_Annot = GO_Annot,
                       GO2Gene = sets_annot_new)
usethis::use_data(GO_Data, overwrite = TRUE)


# Create general terms
# Generate table for defining general terms
roots <- c("GO:0008150", "GO:0005575", "GO:0003674")
for (i in roots) {
  tmp <- unlist2(GO_ontology_new$children[GO_ontology_new$children[[i]]]) %>%
    data.frame("P" = names(.), "C" = .)
  tmp$P_name <- GO_ontology_new$name[as.vector(tmp$P)]
  tmp$C_name <- GO_ontology_new$name[as.vector(tmp$C)]
  write.xlsx(tmp, paste0("data-raw/define_general_terms_ref_",i,".xlsx"))
  tmp <- unlist2(GO_ontology_new$children[i]) %>%
    data.frame("P" = names(.), "C" = .)
  tmp$P_name <- GO_ontology_new$name[as.vector(tmp$P)]
  tmp$C_name <- GO_ontology_new$name[as.vector(tmp$C)]
  write.xlsx(tmp, paste0("data-raw/define_general_terms_",i,".xlsx"))
}

# edit the excels

GeneralTermsGO_mgi <- list()
GetChildOfLevel <- function(GO_term, level, start = 1){
  if(start < level) {
    start = start + 1
    children <- GO_ontology_new$children[GO_term] %>% unlist()
    return(GetChildOfLevel(children, level, start))}
  else
    return(GO_term)
}
for (i in roots) {
  test_2 <- read_excel(paste0("data-raw/define_general_terms_",i,".xlsx")) %>% as.data.frame()
  for (j in rownames(test_2)) {
    GeneralTermsGO_mgi[[i]] <- c(GeneralTermsGO_mgi[[i]], GetChildOfLevel(test_2[j,"C"], test_2[j,"level"]))
  }
  GeneralTermsGO_mgi[[i]] <- unique(GeneralTermsGO_mgi[[i]])
}
usethis::use_data(GeneralTermsGO_mgi, overwrite = TRUE)




################## not processed

library(mgsa)
library(ontologyIndex)
library(ontologyPlot)
library(dplyr)
library(openxlsx)
library(ggplot2)


GeneralSetsList <- list()
itemlist <- as.character(c(8,12))
for(i in itemlist){
  GeneralSetsList[[i]] <- general_terms[lapply(sets_annot_new[general_terms], function(x) i %in% x) %>% unlist()]
}
GeneralSetsList <-
  unlist2(GeneralSetsList) %>%
  data.frame("Genes" = names(.), "Sets" = .) %>%
  split(.$Genes) %>%
  lapply(function(x){mutate(x, rank = 1:nrow(x),
                            Sets_name = paste(x$Sets, GO_ontology_new$name[as.vector(x$Sets)], sep = " "))}) %>%
  list.rbind() %>%
  mutate(Sets_name = factor(.$Sets_name))

ggplot(GeneralSetsList, aes(x = rank, y = Genes, label = Sets, fill = Sets_name)) +
  geom_label() +
  guides(
    fill = guide_legend(
      title = "GO terms",
      override.aes = aes(label = "")
    )
  ) +
  theme(panel.grid = element_blank(),
        axis.title = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks = element_blank(),
        panel.background = element_blank())

LowestLevelSetsList <- list()
Gene2GO <- unlist2(sets_annot_new) %>%
  data.frame("Sets" = names(.), "Genes" = .) %>%
  split(.$Genes) %>%
  lapply(function(x) x$Sets)
itemlist <- as.character(c(8,12))
for (i in itemlist) {
  LowestLevelSetsList[[i]] <-
    setdiff(Gene2GO[[i]],
            GO_ontology_new$ancestors[Gene2GO[[i]]] %>%
              lapply(function(x) x[-length(x)]) %>%
              unlist() %>% unique())

}
LowestLevelSetsList <-
  unlist2(LowestLevelSetsList) %>%
  data.frame("Genes" = names(.), "Sets" = .) %>%
  split(.$Genes) %>%
  lapply(function(x){mutate(x, rank = 1:nrow(x),
                            Sets_name = paste(x$Sets, GO_ontology_new$name[as.vector(x$Sets)], sep = " "))}) %>%
  list.rbind() %>%
  mutate(Sets_name = factor(.$Sets_name))

ggplot(LowestLevelSetsList, aes(x = rank, y = Genes, label = Sets, fill = Sets_name)) +
  geom_label() +
  guides(
    fill = guide_legend(
      title = "GO terms",
      override.aes = aes(label = "")
    )
  ) +
  theme(panel.grid = element_blank(),
        axis.title = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks = element_blank(),
        panel.background = element_blank())
