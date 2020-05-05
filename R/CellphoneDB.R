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
  message(paste(Sys.time(), "Start runing CellphoneDB"))
  system(command)

  deconvoluted <- read.table("out/deconvoluted.txt",sep = "\t", header = T, check.names=FALSE)
  means <- read.table("out/means.txt",sep = "\t", header = T, check.names=FALSE)
  pvalues <- read.table("out/pvalues.txt",sep = "\t", header = T, check.names=FALSE)
  significant_means <- read.table("out/significant_means.txt",sep = "\t", header = T, check.names=FALSE)
  seu@misc[["CellphoneDB"]][[group.by]][["deconvoluted"]] <- deconvoluted
  seu@misc[["CellphoneDB"]][[group.by]][["means"]] <- means
  seu@misc[["CellphoneDB"]][[group.by]][["pvalues"]] <- pvalues
  seu@misc[["CellphoneDB"]][[group.by]][["significant_means"]] <- significant_means
  message(paste0(Sys.time(), " Results saved in the slot: SeuratObject@misc$CellphoneDB$", group.by))
  return(seu)
}

#' @title FUNCTION_TITLE
#' @description FUNCTION_DESCRIPTION
#' @param seu PARAM_DESCRIPTION
#' @param sender PARAM_DESCRIPTION
#' @param receiver PARAM_DESCRIPTION
#' @param group.by PARAM_DESCRIPTION
#' @param top_n PARAM_DESCRIPTION, Default: 20
#' @return OUTPUT_DESCRIPTION
#' @details DETAILS
#' @examples
#' \dontrun{
#' if(interactive()){
#'  #EXAMPLE1
#'  }
#' }
#' @rdname CellphoneDB_Plots
#' @export

CellphoneDB_Plots <- function(seu, sender, receiver, group.by, top_n = 20){
  library(dplyr)
  library(SeuratExtend)
  library(Seurat)
  library(reshape2)

  clu_pairs <-
    data.frame("sender" = rep(sender, each = length(receiver)),
               "receiver" = rep(receiver, n = length(sender)),
               stringsAsFactors = F)
  clu_pairs$pairs <- paste(clu_pairs$sender, clu_pairs$receiver, sep = "|")
  rownames(clu_pairs) <- clu_pairs$pairs

  significant_means_trimmed <-
    seu@misc[["CellphoneDB"]][[group.by]][["significant_means"]] %>%
    `rownames<-`(.$interacting_pair) %>%
    .[,clu_pairs$pairs]  %>%
    .[apply(.,1,function(x) !all(is.na(x))),]
  significant_means_trimmed_rev <-
    seu@misc[["CellphoneDB"]][[group.by]][["significant_means"]] %>%
    `rownames<-`(.$interacting_pair) %>%
    .[,pair_rev(clu_pairs$pairs, sep = "|")]  %>%
    .[apply(.,1,function(x) !all(is.na(x))),] %>%
    `rownames<-`(pair_rev(rownames(.))) %>%
    `colnames<-`(pair_rev(colnames(.), sep = "|"))
  significant_means_trimmed <-
    rbind(significant_means_trimmed,
          significant_means_trimmed_rev[
            setdiff(rownames(significant_means_trimmed_rev),
                    rownames(significant_means_trimmed)),])

  significant_means_trimmed_top <-
    significant_means_trimmed[
      apply(significant_means_trimmed, 1, function(x) max(x, na.rm = T)) %>%
        sort(decreasing = T) %>%
        head(top_n) %>%
        names(),
      ]

  p_gene_clu <- Heatmap(significant_means_trimmed_top, color_scheme = "B")

  clu_clu_count <- data.frame()
  clu_clu_sum <- data.frame()
  for (i in rownames(clu_pairs)) {
    s <- clu_pairs[i,"sender"]
    r <- clu_pairs[i,"receiver"]
    sr <- clu_pairs[i,"pairs"]
    clu_clu_count[s,r] <- sum(!is.na(significant_means_trimmed[,sr]))
    clu_clu_sum[s,r] <- sum(significant_means_trimmed[,sr], na.rm = T)
  }

  p_clu_clu_sum <- Heatmap(clu_clu_sum, color_scheme = "E", lab_fill = "Accumulated\nmean value")
  p_clu_clu_count <- Heatmap(clu_clu_count, color_scheme = "E", lab_fill = "Number of \ninteractions")

  gene_pairs <- data.frame()
  for (i in rownames(significant_means_trimmed_top)) {
    g <- strsplit(i, split = "_")[[1]]
    gene_pairs[g[2], g[1]] <- 1
  }
  gene_pairs[is.na(gene_pairs)] <- 0
  p_gene_pairs <- Heatmap(gene_pairs, lab_fill = "Interaction pairs")

  gene_clu_sender <-
    seu@misc[["CellphoneDB"]][[group.by]][["deconvoluted"]] %>%
    .[!duplicated(.$gene_name),c("gene_name",sender)] %>%
    `rownames<-`(.$gene_name) %>%
    .[colnames(gene_pairs),] %>%
    melt(value.name = "expression", variable.name = "cluster") %>%
    mutate("percent" = lapply(sender, function(x){
      gene_percent(seu, feature = colnames(gene_pairs), group.by = group.by, cluster = x)
    }) %>% unlist)

  p_gene_clu_sender <-
    ggplot(gene_clu_sender, mapping = aes(x = gene_name, y = cluster, size = percent, color = expression)) +
    geom_point() +
    theme_classic() +
    scale_color_gradientn(colors = viridis_pal(option = "B")(20))

  gene_clu_receiver <-
    seu@misc[["CellphoneDB"]][[group.by]][["deconvoluted"]] %>%
    .[!duplicated(.$gene_name),c("gene_name",receiver)] %>%
    `rownames<-`(.$gene_name) %>%
    .[rownames(gene_pairs),] %>%
    melt(value.name = "expression", variable.name = "cluster") %>%
    mutate("percent" = lapply(receiver, function(x){
      gene_percent(seu, feature = rownames(gene_pairs), group.by = group.by, cluster = x)
    }) %>% unlist)

  p_gene_clu_receiver <-
    ggplot(gene_clu_receiver, mapping = aes(x = cluster, y = gene_name, size = percent, color = expression)) +
    geom_point() +
    theme_classic() +
    scale_color_gradientn(colors = viridis_pal(option = "B")(20))

  CellphoneDB_Plots <-
    list(clu_pairs = clu_pairs,
         significant_means_trimmed = significant_means_trimmed,
         significant_means_trimmed_top = significant_means_trimmed_top,
         p_gene_clu = p_gene_clu,
         p_clu_clu_sum = p_clu_clu_sum,
         p_clu_clu_count = p_clu_clu_count,
         gene_pairs = gene_pairs,
         p_gene_pairs = p_gene_pairs,
         gene_clu_sender = gene_clu_sender,
         p_gene_clu_sender = p_gene_clu_sender,
         gene_clu_receiver = gene_clu_receiver,
         p_gene_clu_receiver = p_gene_clu_receiver)
  return(CellphoneDB_Plots)
}

# library(Seurat)
# library(SeuratExtend)
# library(purrr)
#
# setwd("~/R documents/cellphonedb test")
# CellphoneDB_GenerateCustomDB(Reactome_interactions_filtered = Reactome_interactions_filtered)
# seu <- readRDS("~/R documents/2020-2-10 EC PyMT and E0771/rds/PyMTEC_old.rds")
# group.by <- "cluster"
# RunCellphoneDB(seu, group.by)

# seu <- readRDS("~/R documents/2020-2-10 EC PyMT and E0771/cellphonedb/seu2.rds")
# sender <- c("Arterial","Cap Cxcl12+","Cap Cxcl9+","HEV","Lymphatic","Tip cell","Venous")
# receiver <- c("Bcell","DC","Macrophage","NK","pDC","Tcell","Tumor")
# group.by <- "merge_cluster"
# top_n <- 20
# CellphoneDB_Plots(seu, sender, receiver, group.by)
