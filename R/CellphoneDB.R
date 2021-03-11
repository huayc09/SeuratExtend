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
#' @param seu PARAM_DESCRIPTION, Default: NULL
#' @param dir PARAM_DESCRIPTION, Default: 'out'
#' @return OUTPUT_DESCRIPTION
#' @details DETAILS
#' @examples
#' \dontrun{
#' if(interactive()){
#'  #EXAMPLE1
#'  }
#' }
#' @rdname ImportCellphoneDBOutputs
#' @export

ImportCellphoneDBOutputs <- function(seu = NULL, dir = "out"){
  Output <- list()
  Output[["deconvoluted"]] <- read.table(file.path(dir,"deconvoluted.txt"), sep = "\t", header = T, check.names=FALSE)
  Output[["means"]] <- read.table(file.path(dir,"means.txt"), sep = "\t", header = T, check.names=FALSE)
  if("pvalues.txt" %in% list.files(dir)) Output[["pvalues"]] <-
    read.table(file.path(dir,"pvalues.txt"), sep = "\t", header = T, check.names=FALSE)
  Output[["significant_means"]] <- read.table(file.path(dir,"significant_means.txt"), sep = "\t", header = T, check.names=FALSE)
  if(is.null(seu)) return(Output)
  seu@misc[["CellphoneDB"]] <- Output
  return(seu)
}

#' @title FUNCTION_TITLE
#' @description FUNCTION_DESCRIPTION
#' @param seu PARAM_DESCRIPTION
#' @param group.by PARAM_DESCRIPTION
#' @param cpdb.dir PARAM_DESCRIPTION, Default: NULL
#' @param method PARAM_DESCRIPTION, Default: c("analysis", "statistical_analysis")
#' @param database PARAM_DESCRIPTION, Default: 'cellphonedb_mouse_cytokines'
#' @param export.to.Seurat PARAM_DESCRIPTION, Default: T
#' @param NormalizeData PARAM_DESCRIPTION, Default: T
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

RunCellphoneDB <- function(seu, group.by, cpdb.dir = NULL, method = c("analysis", "statistical_analysis"),
                           database = "cellphonedb_mouse_cytokines", export.to.Seurat = T, NormalizeData = T) {
  message("Check if CellphoneDB is installed and able to run in terminal")
  cpdb.path <- ifelse(is.null(cpdb.dir),
                      "cellphonedb",
                      file.path(cpdb.dir,"cellphonedb"))
  test <- system(cpdb.path)
  if(test != 0) {
    stop("For installation, please visit: https://github.com/Teichlab/cellphonedb\n",
         "If CellphoneDB is already installed, please make sure ",
         "you have set the correct direction in 'cpdb.dir' augument (e.g. cpdb-venv/bin/)")
  }

  library(Seurat)
  library(dplyr)

  dir.create("input")
  if(is.null(database)) {
    command <- paste0(cpdb.path, " method ", method[1]," input/meta_data.txt ",
                      "input/counts.txt --counts-data hgnc_symbol")
  }else if(database == "cellphonedb_mouse_cytokines"){
    file.copy(system.file("extdata/CellphoneDB", "cellphonedb_mouse_cytokines.db", package = "SeuratExtend"),
              "input/custom_db.db")
    command <- paste0(cpdb.path, " method ", method[1]," input/meta_data.txt input/counts.txt ",
                      "--database input/custom_db.db")
  }else{
    file.copy(database, "input/custom_db.db")
    command <- paste0(cpdb.path, " method ", method[1]," input/meta_data.txt input/counts.txt ",
                      "--database input/custom_db.db")
  }

  if(NormalizeData) {
    if(!any(grepl("NormalizeData",names(seu@commands)))) seu <- NormalizeData(seu)
    data <- GetAssayData(seu, slot = "data")
  } else {
    data <- GetAssayData(seu, slot = "counts")
  }
  meta <- seu@meta.data
  data_export <- cbind(data.frame("Gene"=rownames(data)), data)
  meta_export <- data.frame("Cell" = rownames(meta), "cell_type" = trimws(meta[,group.by]))
  write.table(meta_export, file = "input/meta_data.txt", row.names = F, sep = "\t", quote = F)
  write.table(data_export, file = "input/counts.txt", row.names = F, sep = "\t", quote = F)
  message(paste(Sys.time(), "Start runing CellphoneDB"))
  system(command)

  if(!export.to.Seurat) return(ImportCellphoneDBOutputs())
  seu <- ImportCellphoneDBOutputs(seu = seu)
  return(seu)
}

#' @title FUNCTION_TITLE
#' @description FUNCTION_DESCRIPTION
#' @param seu PARAM_DESCRIPTION, Default: NULL
#' @param output PARAM_DESCRIPTION, Default: NULL
#' @param sender PARAM_DESCRIPTION, Default: NULL
#' @param receiver PARAM_DESCRIPTION, Default: NULL
#' @param group.by PARAM_DESCRIPTION, Default: NULL
#' @param top_n PARAM_DESCRIPTION, Default: 20
#' @param Dotplot PARAM_DESCRIPTION, Default: F
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

CellphoneDB_Plots <- function(seu = NULL, output = NULL, sender = NULL, receiver = NULL, group.by = NULL,
                              top_n = 20, Dotplot = F){

  if(!is.null(output)) {
    if(is.character(output)) output <- ImportCellphoneDBOutputs(dir = output)
    if(!is.list(output)) stop("Please check 'output' argument")
  }else{
    if(!is.null(seu)) {
      if(!is.null(seu@misc[["CellphoneDB"]])) output <- seu@misc[["CellphoneDB"]] else
        stop("No CellphoneDB outputs found. Please check the arguments")
    } else stop("No Seurat object or CellphoneDB output defined")
  }

  library(dplyr)
  library(Seurat)
  library(reshape2)
  library(ggplot2)
  library(viridis)
  import("egg")
  library(rlang)
  library(rlist)
  library(tidyr)

  clusters <- colnames(output$deconvoluted)[-c(1:6)] %>% trimws()
  sender <- sender %||% clusters %>% trimws()
  receiver <- receiver %||% clusters %>% trimws()
  if(!all(sender %in% clusters) | !all(receiver %in% clusters)){
    message("Detected cell types according to CellphoneDB output: ", paste0(clusters, collapse = ", "))
    if(!all(sender %in% clusters)) message("Unknown cell type in 'sender' argument: ",
                                           paste0(setdiff(sender, clusters), collapse = ", "))
    if(!all(receiver %in% clusters)) message("Unknown cell type in 'receiver' argument: ",
                                             paste0(setdiff(receiver, clusters), collapse = ", "))
    stop("Please revise the sender/receiver argument")
  }
  #seu@meta.data[,group.by] <- trimws(seu@meta.data[,group.by])

  clu_pairs <-
    data.frame("sender" = rep(sender, each = length(receiver)),
               "receiver" = rep(receiver, n = length(sender)),
               stringsAsFactors = F)
  clu_pairs$pairs <- paste(clu_pairs$sender, clu_pairs$receiver, sep = "|")
  rownames(clu_pairs) <- clu_pairs$pairs

  sig <- output[["significant_means"]]
  sig <- sig[!duplicated(sig$interacting_pair),]
  check.diff <- setdiff(c(clu_pairs$pairs, pair_rev(clu_pairs$pairs, sep = "|")) %>% unique(), colnames(sig))
  if(length(check.diff) > 0) sig[,check.diff] <- NA
  significant_means_trimmed <-
    sig %>%
    `rownames<-`(.$interacting_pair) %>%
    .[,clu_pairs$pairs, drop = F]  %>%
    .[apply(.,1,function(x) !all(is.na(x))),, drop = F]
  significant_means_trimmed_rev <-
    sig %>%
    `rownames<-`(.$interacting_pair) %>%
    .[,pair_rev(clu_pairs$pairs, sep = "|"), drop = F]  %>%
    .[apply(.,1,function(x) !all(is.na(x))),, drop = F] %>%
    `rownames<-`(pair_rev(rownames(.))) %>%
    `colnames<-`(pair_rev(colnames(.), sep = "|"))
  significant_means_trimmed <-
    rbind(significant_means_trimmed,
          significant_means_trimmed_rev[
            setdiff(rownames(significant_means_trimmed_rev),
                    rownames(significant_means_trimmed)),, drop = F])

  significant_means_trimmed_top <-
    significant_means_trimmed[
      apply(significant_means_trimmed, 1, function(x) max(x, na.rm = T)) %>%
        sort(decreasing = T) %>%
        head(top_n) %>%
        names(),, drop = F
    ]

  significant_means_trimmed_top$gene_pair <- rownames(significant_means_trimmed_top)
  gene_clu <- melt(significant_means_trimmed_top,
                   id.vars = "gene_pair",
                   variable.name = "cluster",
                   value.name = "means")
  gene_clu[,c("sender","receiver")] <- clu_pairs[gene_clu$cluster,c("sender","receiver")]
  gene_clu$means[is.na(gene_clu$means)] <- 0
  gene_clu$gene_pair <- factor(gene_clu$gene_pair, levels = rev(rownames(significant_means_trimmed_top)))
  gene_clu$cluster <- factor(gene_clu$cluster, levels = colnames(significant_means_trimmed_top) %>% .[-length(.)])

  p_gene_clu <-
    ggplot(gene_clu, aes(x = receiver, y = gene_pair, fill = means)) +
    facet_grid(cols = vars(sender), scales = "free") +
    geom_tile(colour = "white") +
    scale_fill_viridis() +
    theme_classic()+
    labs(x = "", y = "")+
    scale_y_discrete(position = "right") +
    theme(axis.text.x=element_text(angle = 90, hjust = 1, vjust = 0.5))

  clu_clu_count <- data.frame()
  clu_clu_sum <- data.frame()
  for (i in rownames(clu_pairs)) {
    s <- clu_pairs[i,"sender"]
    r <- clu_pairs[i,"receiver"]
    sr <- clu_pairs[i,"pairs"]
    clu_clu_count[s,r] <- sum(!is.na(significant_means_trimmed[,sr]))
    clu_clu_sum[s,r] <- sum(significant_means_trimmed[,sr], na.rm = T)
  }

  gene_pairs <- data.frame()
  for (i in rownames(significant_means_trimmed_top)) {
    g <- strsplit(i, split = "_")[[1]]
    gene_pairs[g[2], g[1]] <- 1
  }
  for(i in rownames(significant_means_trimmed)) {
    a <- strsplit(i, split = "_")[[1]][1]
    b <- strsplit(i, split = "_")[[1]][2]
    if((a %in% rownames(gene_pairs)) & (b %in% colnames(gene_pairs))) gene_pairs[a,b] <- 1
  }
  gene_pairs[is.na(gene_pairs)] <- 0
  p_gene_pairs <- Heatmap(gene_pairs, lab_fill = "Interaction pairs")

  dec <- output[["deconvoluted"]]
  dec_sender <-
    dec[dec$gene_name %in% colnames(gene_pairs) |
          dec$complex_name %in% colnames(gene_pairs),
        c("gene_name","complex_name",sender)] %>%
    unique() %>%
    split(.$complex_name %>% replace_na(""))
  for (i in setdiff(names(dec_sender), "")) {
    dec_sender[[i]] <-
      cbind(data.frame(gene_name = i,
                       complex_name = i),
            apply(dec_sender[[i]][,-c(1:2), drop = F], 2, mean) %>% t)
  }
  dec_sender <- list.rbind(dec_sender)
  dec_receiver <-
    dec[dec$gene_name %in% rownames(gene_pairs) |
          dec$complex_name %in% rownames(gene_pairs),
        c("gene_name","complex_name",receiver)] %>%
    unique() %>%
    split(.$complex_name %>% replace_na(""))
  for (i in setdiff(names(dec_receiver), "")) {
    dec_receiver[[i]] <-
      cbind(data.frame(gene_name = i,
                       complex_name = i),
            apply(dec_receiver[[i]][,-c(1:2), drop = F], 2, mean) %>% t)
  }
  dec_receiver <- list.rbind(dec_receiver)

  if(Dotplot & !is.null(seu) & !is.null(group.by) &
     all(dec_sender$complex_name == "") &
     all(dec_receiver$complex_name == "")){
    gene_clu_sender <-
      dec_sender[,-2] %>%
      melt(value.name = "expression", variable.name = "cluster") %>%
      mutate("percent" = lapply(sender, function(x){
        gene_percent(seu, feature = colnames(gene_pairs), ident = x, group.by = group.by)
      }) %>% unlist)
    gene_clu_sender$gene_name <- factor(gene_clu_sender$gene_name, levels = colnames(gene_pairs))
    gene_clu_sender$cluster <- factor(gene_clu_sender$cluster, levels = rev(sender))

    p_gene_clu_sender <-
      ggplot(gene_clu_sender, mapping = aes(x = gene_name, y = cluster, size = percent, color = expression)) +
      geom_point() +
      theme_classic() +
      scale_color_gradientn(colors = viridis_pal(option = "D")(20))

    gene_clu_receiver <-
      dec_receiver[,-2] %>%
      melt(value.name = "expression", variable.name = "cluster") %>%
      mutate("percent" = lapply(receiver, function(x){
        gene_percent(seu, feature = rownames(gene_pairs), ident = x, group.by = group.by)
      }) %>% unlist)
    gene_clu_receiver$gene_name <- factor(gene_clu_receiver$gene_name, levels = rev(rownames(gene_pairs)))

    p_gene_clu_receiver <-
      ggplot(gene_clu_receiver, mapping = aes(x = cluster, y = gene_name, size = percent, color = expression)) +
      geom_point() +
      theme_classic() +
      scale_color_gradientn(colors = viridis_pal(option = "D")(20))
    p_clu_clu_sum <- Heatmap(clu_clu_sum, color_scheme = "D", lab_fill = "Accumulated\nmean value")
    p_clu_clu_count <- Heatmap(clu_clu_count, color_scheme = "D", lab_fill = "Number of \ninteractions")
  }else{
    gene_clu_sender <-
      dec_sender %>%
      `rownames<-`(.$gene_name) %>%
      .[colnames(gene_pairs),sender, drop = F]
    p_gene_clu_sender <- Heatmap(t(gene_clu_sender), color_scheme = c("white",muted("red")))
    gene_clu_receiver <-
      dec_receiver %>%
      `rownames<-`(.$gene_name) %>%
      .[rownames(gene_pairs),receiver, drop = F]
    p_gene_clu_receiver <- Heatmap(gene_clu_receiver, color_scheme = c("white",muted("blue")))
    p_clu_clu_sum <- Heatmap(clu_clu_sum, color_scheme = c("white",muted("magenta")), lab_fill = "Accumulated\nmean value")
    p_clu_clu_count <- Heatmap(clu_clu_count, color_scheme = c("white",muted("magenta")), lab_fill = "Number of \ninteractions")
  }

  p_combined <-
    ggarrange(p_clu_clu_sum +
                scale_y_discrete(position = "left") +
                scale_x_discrete(position = "top") +
                theme(legend.position = "none",
                      axis.text.x = element_text(angle = 90, hjust = 0, vjust = 0.5)),
              p_gene_clu_sender +
                scale_y_discrete(position = "right") +
                scale_x_discrete(position = "top") +
                theme(legend.position = "none",
                      axis.text.x = element_text(angle = 90, hjust = 0, vjust = 0.5),
                      axis.title = element_blank(),
                      axis.line = element_blank(),
                      axis.ticks = element_blank()),
              p_gene_clu_receiver +
                scale_y_discrete(position = "left") +
                scale_x_discrete(position = "bottom") +
                theme(legend.position = "none",
                      axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
                      axis.title = element_blank(),
                      axis.line = element_blank(),
                      axis.ticks = element_blank()),
              p_gene_pairs +
                scale_fill_gradient(low = "white", high = "#238A8DFF") +
                theme(legend.position = "none",
                      axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)),
              widths = c(length(receiver), ncol(gene_pairs)),
              heights = c(length(sender), nrow(gene_pairs)))

  CellphoneDB_Plots <-
    list(sender = sender,
         receiver = receiver,
         clu_pairs = clu_pairs,
         significant_means_trimmed = significant_means_trimmed,
         significant_means_trimmed_top = significant_means_trimmed_top,
         gene_clu = gene_clu,
         p_gene_clu = p_gene_clu,
         p_clu_clu_sum = p_clu_clu_sum,
         p_clu_clu_count = p_clu_clu_count,
         gene_pairs = gene_pairs,
         p_gene_pairs = p_gene_pairs,
         gene_clu_sender = gene_clu_sender,
         p_gene_clu_sender = p_gene_clu_sender,
         gene_clu_receiver = gene_clu_receiver,
         p_gene_clu_receiver = p_gene_clu_receiver,
         p_combined = p_combined)
  return(CellphoneDB_Plots)
}

#' @title FUNCTION_TITLE
#' @description FUNCTION_DESCRIPTION
#' @param seu PARAM_DESCRIPTION, Default: NULL
#' @param output PARAM_DESCRIPTION, Default: NULL
#' @param sender PARAM_DESCRIPTION, Default: NULL
#' @param receiver PARAM_DESCRIPTION, Default: NULL
#' @param nLink PARAM_DESCRIPTION, Default: 50
#' @param ignore PARAM_DESCRIPTION, Default: NULL
#' @param alpha.low.link PARAM_DESCRIPTION, Default: 0.1
#' @return OUTPUT_DESCRIPTION
#' @details DETAILS
#' @examples
#' \dontrun{
#' if(interactive()){
#'  #EXAMPLE1
#'  }
#' }
#' @rdname CellphoneDB_Plots_Circlize
#' @export

CellphoneDB_Plots_Circlize <-
  function(seu = NULL, output = NULL, sender = NULL, receiver = NULL, nLink = 50,
           ignore = NULL, alpha.low.link = 0.1){
    if(!is.null(output)) {
      if(is.character(output)) output <- ImportCellphoneDBOutputs(dir = output)
      if(!is.list(output)) stop("Please check 'output' argument")
    }else{
      if(!is.null(seu)) {
        if(!is.null(seu@misc[["CellphoneDB"]])) output <- seu@misc[["CellphoneDB"]] else
          stop("No CellphoneDB outputs found. Please check the arguments")
      } else stop("No Seurat object or CellphoneDB output defined")
    }

    library(dplyr)
    library(Seurat)
    library(reshape2)
    library(ggplot2)
    library(viridis)
    library(rlang)
    library(rlist)
    import("circlize")
    library(scales)
    library(tidyr)
    import("ComplexHeatmap", detach = T)

    clusters <- colnames(output$deconvoluted)[-c(1:6)] %>% trimws()
    sender <- sender %||% clusters %>% trimws()
    receiver <- receiver %||% clusters %>% trimws()
    if(!all(sender %in% clusters) | !all(receiver %in% clusters)){
      message("Detected cell types according to CellphoneDB output: ", paste0(clusters, collapse = ", "))
      if(!all(sender %in% clusters)) message("Unknown cell type in 'sender' argument: ",
                                             paste0(setdiff(sender, clusters), collapse = ", "))
      if(!all(receiver %in% clusters)) message("Unknown cell type in 'receiver' argument: ",
                                               paste0(setdiff(receiver, clusters), collapse = ", "))
      stop("Please revise the sender/receiver argument")
    }

    clu_pairs <-
      data.frame("sender" = rep(sender, each = length(receiver)),
                 "receiver" = rep(receiver, n = length(sender)),
                 stringsAsFactors = F)
    clu_pairs$pairs <- paste(clu_pairs$sender, clu_pairs$receiver, sep = "|")
    rownames(clu_pairs) <- clu_pairs$pairs

    sig <- output[["significant_means"]]
    sig <- sig[!duplicated(sig$interacting_pair),]
    if(!is.null(ignore)) sig <- sig[!sig$interacting_pair %in% ignore, ]
    significant_means_trimmed <-
      sig %>%
      `rownames<-`(.$interacting_pair) %>%
      .[,clu_pairs$pairs, drop = F]  %>%
      .[apply(.,1,function(x) !all(is.na(x))),, drop = F]
    significant_means_trimmed_rev <-
      sig %>%
      `rownames<-`(.$interacting_pair) %>%
      .[,pair_rev(clu_pairs$pairs, sep = "|"), drop = F]  %>%
      .[apply(.,1,function(x) !all(is.na(x))),, drop = F] %>%
      `rownames<-`(pair_rev(rownames(.))) %>%
      `colnames<-`(pair_rev(colnames(.), sep = "|"))
    significant_means_trimmed <-
      rbind(significant_means_trimmed,
            significant_means_trimmed_rev[
              setdiff(rownames(significant_means_trimmed_rev),
                      rownames(significant_means_trimmed)),, drop = F])

    thr <- significant_means_trimmed %>% as.matrix %>% .[order(.,decreasing = T)[nLink]]
    significant_means_trimmed_top <-
      significant_means_trimmed %>%
      .[apply(., 1, function(x) max(x, na.rm = T) >= thr),
        apply(., 2, function(x) max(x, na.rm = T) >= thr), drop = F]
    LinkDf <-
      melt(cbind(significant_means_trimmed_top,
                 Interaction = rownames(significant_means_trimmed_top)),
           id.vars = "Interaction") %>%
      .[!is.na(.$value),] %>%
      separate(col = "Interaction", into = c("gene_a","gene_b"), sep = "_", remove = F) %>%
      separate(col = "variable", into = c("CellType1","CellType2"), sep = "[|]")

    dec <- output[["deconvoluted"]]
    dec <-
      dec[, c("gene_name","complex_name",union(sender, receiver))] %>%
      unique() %>%
      split(.$complex_name)
    for (i in setdiff(names(dec), "")) {
      dec[[i]] <-
        cbind(data.frame(gene_name = i,
                         complex_name = i),
              apply(dec[[i]][,-c(1:2), drop = F], 2, mean) %>% t)
    }
    dec <- list.rbind(dec)

    FactorDf <- data.frame()
    for (i in rownames(LinkDf)) {
      FactorDf[paste0(i,"a"),"Celltype"]<-LinkDf[i,"CellType1"]
      FactorDf[paste0(i,"a"),"Gene"]<-LinkDf[i,"gene_a"]
      FactorDf[paste0(i,"a"),"id"]<-paste0(LinkDf[i,"CellType1"],"_", LinkDf[i,"gene_a"])
      FactorDf[paste0(i,"b"),"Celltype"]<-LinkDf[i,"CellType2"]
      FactorDf[paste0(i,"b"),"Gene"]<-LinkDf[i,"gene_b"]
      FactorDf[paste0(i,"b"),"id"]<-paste0(LinkDf[i,"CellType2"],"_", LinkDf[i,"gene_b"])
    }
    FactorDf <- FactorDf %>% .[!duplicated(.$id),]
    for (i in rownames(FactorDf)) {
      FactorDf[i,"Level"] <- dec[dec$gene_name==FactorDf[i,"Gene"], FactorDf[i,"Celltype"]]
    }
    FactorDf$Celltype <- factor(FactorDf$Celltype, levels = unique(c(sender,receiver)) %>% intersect(FactorDf$Celltype))
    FactorDf <- FactorDf %>%
      .[order(.$Level, decreasing = T),] %>%
      .[order(.$Celltype),]
    for (i in levels(FactorDf$Celltype)) {
      Index<-0
      for (j in rownames(FactorDf)) {
        if(FactorDf[j,"Celltype"]==i){
          Index<-Index+1
          FactorDf[j,"Index"]<-Index
        }
      }
    }
    FactorDf$col<-viridis_pal(option = "D")(20)[as.numeric(cut(FactorDf$Level,breaks = 20))]

    xlimDf<-data.frame(row.names = levels(FactorDf$Celltype))
    for (i in rownames(xlimDf)) {
      xlimDf[i,"min"]<-0
      xlimDf[i,"max"]<-max(FactorDf$Index[FactorDf$Celltype==i])
    }
    circos.initialize(factors = FactorDf$Celltype, x = FactorDf$Index, xlim = xlimDf)
    circos.par(points.overflow.warning = F)
    circos.par(track.margin=c(0.01,0.2))
    circos.track(ylim = c(0, 1), track.height = 0.1, panel.fun = function(x, y) {
      xlim = CELL_META$xlim
      ylim = CELL_META$ylim
      tmp<-FactorDf[FactorDf$Celltype==CELL_META$sector.index,]
      circos.rect(tmp$Index-1, rep(0, nrow(tmp)),
                  tmp$Index, rep(1, nrow(tmp)),
                  col = tmp$col, border = NA)
      circos.text(tmp$Index-0.5, rep(1.4, nrow(tmp)),
                  FactorDf$Gene[FactorDf$Celltype==CELL_META$sector.index],
                  facing = "clockwise", niceFacing = TRUE,
                  adj = c(0, 0), cex = 0.6)
    })
    col <- gg_color_hue(nlevels(FactorDf$Celltype))
    CelltypeCol <- data.frame(Celltype=unique(FactorDf$Celltype))
    CelltypeCol$col <- sample(col)

    circos.par(track.margin=c(0.01,0.01))
    circos.track(factors = CelltypeCol$Celltype, ylim=c(0,1),
                 bg.col = CelltypeCol$col,
                 track.height = 0.1,
                 panel.fun = function(x, y) {
                   circos.text(CELL_META$xcenter, 0.5,
                               CELL_META$sector.index, cex = 0.8, col = "white")
                 })


    for (i in rownames(LinkDf)) {
      LinkDf[i,"col"]<-alpha(CelltypeCol[CelltypeCol$Celltype==LinkDf[i,"CellType1"],"col"],
                             alpha = (LinkDf[i, "value"]-min(LinkDf$value))/(max(LinkDf$value)-min(LinkDf$value))*
                               (1-alpha.low.link) + alpha.low.link)
    }
    for (i in rownames(LinkDf)) {
      circos.link(LinkDf[i,"CellType1"],
                  FactorDf[FactorDf$Celltype==LinkDf[i,"CellType1"] &
                             FactorDf$Gene==LinkDf[i,"gene_a"],
                           "Index"] -0.5,
                  LinkDf[i,"CellType2"],
                  FactorDf[FactorDf$Celltype==LinkDf[i,"CellType2"] &
                             FactorDf$Gene==LinkDf[i,"gene_b"],
                           "Index"] -0.5,
                  col = LinkDf[i,"col"],
                  lwd = 2
      )
    }
    labs <- levels(cut(FactorDf$Level,breaks = 20))
    labs <- cbind(lower = as.numeric( sub("\\((.+),.*", "\\1", labs) ),
                upper = as.numeric( sub("[^,]*,([^]]*)\\]", "\\1", labs) ))
    labs <- apply(labs, 1, mean)
    col_fun <- colorRamp2(labs, viridis_pal(option = "D")(20))
    lgd_gene <- ComplexHeatmap::Legend(col_fun = col_fun, title = "Gene Expression",
                                       direction = "horizontal", legend_width = unit(20, "mm"))
    pushViewport(viewport(x = unit(1, "npc"), y = unit(0, "npc"), width = 0.25, height = 0.2, just = c("right", "bottom")))
    ComplexHeatmap::draw(lgd_gene)
  }

# library(Seurat)
# library(SeuratExtend)
# library(purrr)

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
# p <- CellphoneDB_Plots(seu, sender, receiver, group.by)
