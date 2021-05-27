PanglaoDB <-
  read.csv("/Documents/R documents/SeuratExtend_databases/2020-3-27 PanglaoDB/PanglaoDB_markers_27_Mar_2020.tsv",
           sep = "\t")
PanglaoDB_data <- list()
Hs <- grepl("Hs",PanglaoDB$species)
PanglaoDB_data[["marker_list_human"]] <-
  split(PanglaoDB$official.gene.symbol[Hs], PanglaoDB$cell.type[Hs])
Mm <- grepl("Mm",PanglaoDB$species)
PanglaoDB_data[["marker_list_mouse"]] <-
  split(PanglaoDB$official.gene.symbol[Mm], PanglaoDB$cell.type[Mm]) %>%
  lapply(HumanToMouseGenesymbol, match = F)
PanglaoDB_data[["original_data"]] <- PanglaoDB

object.size(PanglaoDB_data)
