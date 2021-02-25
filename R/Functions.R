#' @title FUNCTION_TITLE
#' @description FUNCTION_DESCRIPTION
#' @param pairs PARAM_DESCRIPTION
#' @param sep PARAM_DESCRIPTION, Default: '_'
#' @return OUTPUT_DESCRIPTION
#' @details DETAILS
#' @examples
#' \dontrun{
#' if(interactive()){
#'  #EXAMPLE1
#'  }
#' }
#' @rdname pair_rev
#' @export

pair_rev <- function(pairs, sep = "_") {
  library(dplyr)
  pairs_rev <-
    strsplit(pairs, split = sep, fixed = T) %>%
    sapply(function(x) paste(rev(x), collapse = sep))
  return(pairs_rev)
}

#' @title FUNCTION_TITLE
#' @description FUNCTION_DESCRIPTION
#' @param seu PARAM_DESCRIPTION
#' @param feature PARAM_DESCRIPTION
#' @param ident PARAM_DESCRIPTION
#' @param group.by PARAM_DESCRIPTION, Default: NULL
#' @param DefaultAssay PARAM_DESCRIPTION, Default: 'RNA'
#' @return OUTPUT_DESCRIPTION
#' @details DETAILS
#' @examples
#' \dontrun{
#' if(interactive()){
#'  #EXAMPLE1
#'  }
#' }
#' @rdname gene_percent
#' @export

gene_percent <- function(seu, feature, ident, group.by = NULL, DefaultAssay = "RNA") {
  library(Seurat)
  DefaultAssay(seu) <- DefaultAssay
  if(is.null(group.by)) {
    cells <- names(Idents(seu))[Idents(seu) == ident]
  }else{
    cells <- colnames(seu)[seu[[group.by]] == ident]
  }
  m <- FetchData(seu, vars = feature, cells = cells)
  return(apply(m, 2, function(x) sum(x>0)) / nrow(m))
}

#' @title FUNCTION_TITLE
#' @description FUNCTION_DESCRIPTION
#' @param ident PARAM_DESCRIPTION
#' @param seu PARAM_DESCRIPTION
#' @param feature PARAM_DESCRIPTION, Default: rownames(seu)
#' @param group.by PARAM_DESCRIPTION, Default: NULL
#' @param pct PARAM_DESCRIPTION, Default: 0.1
#' @param DefaultAssay PARAM_DESCRIPTION, Default: 'RNA'
#' @return OUTPUT_DESCRIPTION
#' @details DETAILS
#' @examples
#' \dontrun{
#' if(interactive()){
#'  #EXAMPLE1
#'  }
#' }
#' @rdname gene_expressed
#' @export

gene_expressed <- function(ident, seu, feature = rownames(seu), group.by = NULL, pct = 0.1,
                           DefaultAssay = "RNA") {
  gene_percent <- gene_percent(seu = seu, feature = feature, ident = ident, group.by = group.by,
                               DefaultAssay = DefaultAssay)
  return(names(gene_percent)[gene_percent >= pct])
}

#' @title FUNCTION_TITLE
#' @description FUNCTION_DESCRIPTION
#' @param fun PARAM_DESCRIPTION
#' @param next.line PARAM_DESCRIPTION, Default: T
#' @param ... PARAM_DESCRIPTION
#' @return OUTPUT_DESCRIPTION
#' @details DETAILS
#' @examples
#' \dontrun{
#' if(interactive()){
#'  #EXAMPLE1
#'  }
#' }
#' @rdname write.py.fun
#' @export

write.py.fun <- function(fun, next.line = T, ...) {
  par <- c(...)
  par <- paste(names(par), par, sep = " = ", collapse = ", ")
  text <- paste0(fun, "(", par, ")")
  if(next.line) text <- paste0(text, "\n")
  return(text)
}


{
RelationPlot<-function(nodes, relation){
  require(Rgraphviz)
  rEG <- new("graphNEL", nodes=nodes, edgemode="directed")
  relation_subset<-relation[relation$V1 %in% nodes,]
  relation_subset$V1<-as.vector(relation_subset$V1)
  relation_subset$V2<-as.vector(relation_subset$V2)
  for (i in rownames(relation_subset)) {
    rEG <- addEdge(relation_subset[i,"V1"], relation_subset[i,"V2"], rEG, 1)
  }
  return(rEG)
}
PlotHierachyWithColor<-function(nodes, relations, top_node, scored_matr, r, g, b){
  if(top_node=="Im"){
    TN<-"R-MMU-168256"
  }else if(top_node=="Sg"){
    TN<-"R-MMU-162582"
  }else if(top_node=="Mt"){
    TN<-"R-MMU-1430728"
  }else{
    TN<-top_node
  }
  all_nodes<-GetAllChild(TN, relations)
  Relations_sub<-relations[relations$V1 %in% all_nodes,]
  plot_nodes<-GetAllParent(nodes, Relations_sub)
  rEG<-RelationPlot(all_nodes, relations)
  matr<-scored_matr
  matr[plot_nodes[!plot_nodes %in% rownames(matr)],]<-0
  matr<-matr[plot_nodes,]
  m<-max(matr)
  transform_for_rgb<-function(x){
    if(x<0) n<-0
    else n <- x/m
  }
  for (i in rownames(matr)) {
    matr[i,"r"]<-transform_for_rgb(matr[i,r])
    matr[i,"g"]<-transform_for_rgb(matr[i,g])
    matr[i,"b"]<-transform_for_rgb(matr[i,b])
    matr[i,"rgb"]<-rgb(1-matr[i,"g"]/2-matr[i,"b"]/2,
                       1-matr[i,"r"]/2-matr[i,"b"]/2,
                       1-matr[i,"r"]/2-matr[i,"g"]/2)
  }
  nAttrs<-list()
  nAttrs$fillcolor <- matr$rgb
  names(nAttrs$fillcolor) <- rownames(matr)
  nAttrs$label<-ChangeName(rownames(matr),Ensembl2ReactomeAll, "V2","V4")
  names(nAttrs$label)<-rownames(matr)
  nAttrs$fontsize<-rep(80,nrow(matr))
  names(nAttrs$fontsize)<-rownames(matr)
  rEG<-subGraph(plot_nodes,rEG)
  return(plot(rEG, nodeAttrs=nAttrs))
}
check_spe <- function(spe){
  if(is.null(spe)) stop("species undefined: options(spe = c(\"mouse\", \"human\"))")
}
}
{
  # ScoreAndOrderDown<-function(matr, meta, group, score, order_method="none", n=0){
  #   tmp<-as.matrix(CalcScoreGeneral(matr, meta, group, score))
  #   if(order_method=="none") return(tmp)
  #   tmp_min<-vector()
  #   for (i in rownames(tmp)) {
  #     tmp_min[i]<-colnames(tmp)[tmp[i,]==min(tmp[i,])]
  #   }
  #   tmp_cluster<-colnames(tmp) %>% .[. %in% tmp_min]
  #   if(order_method=="value"){
  #     tmp2<-list()
  #     for (i in tmp_cluster) {
  #       tmp2[[i]]<-tmp[tmp_min==i,] %>% .[order(.[,i]),]
  #     }
  #   }
  #   if(order_method=="sd"){
  #     tmp2<-list()
  #     for (i in tmp_cluster) {
  #       tmp2[[i]]<-tmp[tmp_min==i,]
  #       tmp_sd<-apply(tmp2[[i]],1,sd)
  #       tmp2[[i]]<-tmp2[[i]] %>% .[order(tmp_sd),]
  #     }
  #   }
  #   if(order_method=="p"){
  #     tmp2<-list()
  #     tmp_p<-CalcScoreGeneral(matr, meta, group, "p")
  #     for (i in tmp_cluster) {
  #       tmp2[[i]]<-as.data.frame(tmp)[tmp_min==i,]
  #       tmp_p_i<-tmp_p[rownames(tmp2[[i]]),i]
  #       tmp2[[i]]<-tmp2[[i]] %>% .[order(tmp_p_i),]
  #     }
  #   }
  #   if(n==0){
  #     tmp3<-tmp2[[names(tmp2)[1]]]
  #     for (i in names(tmp2)[2:length(tmp2)]) {
  #       tmp3<-rbind(tmp3, tmp2[[i]])
  #     }
  #   }else{
  #     tmp3<-head(tmp2[[names(tmp2)[1]]],n)
  #     for (i in names(tmp2)[2:length(tmp2)]) {
  #       tmp3<-rbind(tmp3, head(tmp2[[i]],n))
  #     }
  #   }
  #   return(tmp3)
  # }
  # EnsemblToGenesymbol<-function(genelist_list){
  #   GeneListEnsembl<-unique(unlist(genelist_list))
  #   require("biomaRt")
  #   mart <- useDataset("mmusculus_gene_ensembl", useMart("ensembl"))
  #   EnsemblGenes<-getBM(attributes = c("mgi_symbol", "ensembl_gene_id"),
  #                       filters = "ensembl_gene_id",
  #                       values = GeneListEnsembl,
  #                       mart = mart)
  #   GeneListGenesymbol<-list()
  #   for (i in names(genelist_list)) {
  #     GeneListGenesymbol[[i]]<-unique(EnsemblGenes$mgi_symbol[EnsemblGenes$ensembl_gene_id %in% genelist_list[[i]]])
  #     if (length(GeneListGenesymbol[[i]])==0){
  #       GeneListGenesymbol[[i]] <- NULL
  #     }
  #   }
  #   return(GeneListGenesymbol)
  # }
  # GetChild<-function(MMU,relation){
  #   child<-unique(as.vector(relation$V2)[relation$V1 %in% MMU])
  #   return(child)
  # }
  # GetAllChild<-function(MMU,relation){
  #   tmp<-MMU
  #   child<-GetChild(MMU,relation)
  #   tmp<-unique(c(tmp,child))
  #   if(length(tmp)>length(MMU)){
  #     return(GetAllChild(tmp,relation))
  #   }else{
  #     return(tmp)
  #   }
  # }
  # ChangeName<-function(char, ref, from, to){
  #   tmp<-char
  #   for (i in 1:length(char)) {
  #     tmp[i]<-as.character(ref[ref[,from]==char[i],to][1])
  #   }
  #   return(tmp)
  # }
  # GetParent<-function(MMU,relation){
  #   parent<-unique(as.vector(relation$V1)[relation$V2 %in% MMU])
  #   return(parent)
  # }
  # GetAllParent<-function(MMU,relation){
  #   tmp<-MMU
  #   parent<-GetParent(MMU,relation)
  #   tmp<-unique(c(tmp,parent))
  #   if(length(tmp)>length(MMU)){
  #     return(GetAllParent(tmp,relation))
  #   }else{
  #     return(tmp)
  #   }
  # }
  # stacked_violin<-function(matr, meta, group, genelist){
  #   require(ggplot2)
  #   tmp_meta<-meta
  #   for (i in genelist) {
  #     tmp_meta[,i]<-matr[i,]
  #   }
  #   tmp<-reshape2::melt(tmp_meta[,c(group,genelist)])
  #   p<-ggplot(tmp, aes(x=tmp[,group],y=value, fill=tmp[,group]))+
  #     geom_violin(scale = "width")+
  #     facet_wrap( ~variable, ncol = 1, strip.position="left")+
  #     ylab(NULL) +
  #     xlab(NULL) +
  #     theme(strip.background = element_blank(),
  #           strip.placement = "outside",
  #           legend.position = "none",
  #           axis.text.x=element_text(angle = 45,hjust = 1))
  #   return(p)
  # }
  # AUC_violin<-function(Obj, genelist_list){
  #   tmp_meta<-Obj@meta.data
  #   tmp_data<-Obj@data
  #   require(AUCell)
  #   cells_rankings <- AUCell_buildRankings(tmp_data, plotStats=TRUE)
  #   for (i in names(genelist_list)){
  #     AUC <- getAUC(AUCell_calcAUC(genelist_list[[i]], cells_rankings))
  #     tmp_meta[,i]<-AUC[1,]
  #   }
  #   tmp<-reshape2::melt(tmp_meta[,c("cluster",names(genelist_list))])
  #   p<-ggplot(tmp, aes(x=tmp[,"cluster"],y=value, fill=cluster))+
  #     geom_violin(scale = "width")+
  #     facet_wrap( ~variable, ncol = 1, strip.position="left")+
  #     ylab(NULL) +
  #     xlab(NULL) +
  #     theme(strip.background = element_blank(),
  #           strip.placement = "outside",
  #           legend.position = "none",
  #           axis.text.x=element_text(angle = 45,hjust = 1))
  #   return(p)
  # }
  # AUC_matrix<-function(matr, genelist_list, ratio=0.2){
  #   require(AUCell)
  #   cells_rankings <- AUCell_buildRankings(matr, plotStats=TRUE)
  #   tmp<-genelist_list
  #   for (i in names(tmp)){
  #     if(sum(tmp[[i]] %in% rownames(matr))/length(tmp[[i]]) < (1-ratio)){
  #       tmp[[i]]<-NULL
  #     }
  #   }
  #   AUC <- getAUC(AUCell_calcAUC(tmp, cells_rankings, nCores = 2))
  #   return(AUC)
  # }
  # stacked_violin_v3<-function(Seu, group = "seurat_clusters", genelist, cell = NULL, ncol = 1){
  #   require(ggplot2)
  #   require(rlang)
  #   require(dplyr)
  #   cell <- cell %||% colnames(Seu)
  #   genelist <- genelist %>% .[. %in% rownames(Seu)]
  #   matr <- GetAssayData(Seu)[genelist, cell]
  #   tmp_meta<-Seu@meta.data[cell,]
  #   for (i in genelist) {
  #     tmp_meta[,i]<-matr[i,]
  #   }
  #   tmp<-reshape2::melt(tmp_meta[,c(group,genelist)])
  #   p<-ggplot(tmp, aes(x=tmp[,group],y=value, fill=tmp[,group]))+
  #     geom_violin(scale = "width")+
  #     facet_wrap( ~variable, ncol = ncol, strip.position="left")+
  #     ylab(NULL) +
  #     xlab(NULL) +
  #     theme(strip.background = element_blank(),
  #           strip.placement = "outside",
  #           legend.position = "none",
  #           axis.text.x=element_text(angle = 45,hjust = 1)) +
  #     theme_classic() +
  #     labs(fill=group)
  #   return(p)
  # }
  # HumanToMouseGeneList<-function(genelist_list){
  #   GeneListHuman<-unique(unlist(genelist_list))
  #   require("biomaRt")
  #   human = useMart("ensembl", dataset = "hsapiens_gene_ensembl")
  #   mouse = useMart("ensembl", dataset = "mmusculus_gene_ensembl")
  #   genelists = getLDS(attributes = c("hgnc_symbol"), filters = "hgnc_symbol", values = GeneListHuman , mart = human,
  #                      attributesL = c("mgi_symbol"), martL = mouse, uniqueRows=T)
  #   GeneListMouse<-list()
  #   for (i in names(genelist_list)) {
  #     GeneListMouse[[i]]<-unique(genelists$MGI.symbol[genelists$HGNC.symbol %in% genelist_list[[i]]])
  #     if (length(GeneListMouse[[i]])==0){
  #       GeneListMouse[[i]] <- NULL
  #     }
  #   }
  #   return(GeneListMouse)
  # }
  # MouseToHumanGeneExpression<-function(Matr){
  #   require("biomaRt")
  #   human = useMart("ensembl", dataset = "hsapiens_gene_ensembl")
  #   mouse = useMart("ensembl", dataset = "mmusculus_gene_ensembl")
  #   genelists = getLDS(attributes = c("mgi_symbol"), filters = "mgi_symbol", values = rownames(Matr) , mart = mouse,
  #                      attributesL = c("hgnc_symbol"), martL = human, uniqueRows=T)
  #   genes<-genelists$HGNC.symbol %>% .[!duplicated(.)]
  #   cells<-colnames(Matr)
  #   new_matrix<-matrix(nrow = length(genes), ncol = length(cells))
  #   rownames(new_matrix)<-genes
  #   colnames(new_matrix)<-cells
  #   mouse_matrix<-as.matrix(Matr)
  #   for (i in genes){
  #     gene_mouse<-genelists$MGI.symbol[genelists$HGNC.symbol==i][1]
  #     new_matrix[i,]<-mouse_matrix[gene_mouse,]
  #   }
  #   return(new_matrix)
  # }
  # GetGeneSet<-function(genelist_list,char,n=0){
  #   hits<-sum(grepl(char,names(genelist_list)))
  #   if(hits==0) return("no hits")
  #   if(hits>1&n==0) return(names(genelist_list)[grep(char,names(genelist_list))])
  #   if(hits==1) return(genelist_list[[names(genelist_list)[grep(char,names(genelist_list))]]])
  #   if(hits>1&n>0) return(genelist_list[[names(genelist_list)[grep(char,names(genelist_list))][n]]])
  # }
  # GenesetOfGenes<-function(geneset_list, genes, type="onlygene"){
  #   tmp_list<-list()
  #   if(type=="all"){
  #     for (i in names(geneset_list)) {
  #       if(sum(geneset_list[[i]] %in% genes)>0) tmp_list[[i]]<-geneset_list[[i]]
  #     }
  #   }else{
  #     for (i in names(geneset_list)) {
  #       if(sum(geneset_list[[i]] %in% genes)>0) tmp_list[[i]]<-geneset_list[[i]] %>% .[. %in% genes]
  #     }
  #   }
  #   return(tmp_list)
  # }
  #
  # AddMetaData<-function(Df, MetadataDf, ID, col){
  #   MetaIDs<-as.data.frame(MetadataDf)[,ID]
  #   for (i in rownames(Df)) {
  #     if(Df[i,ID] %in% MetaIDs){
  #       for (j in col) {
  #         Df[i,j]<-MetadataDf[MetadataDf[,ID]==Df[i,ID],j][1]
  #       }
  #     }
  #   }
  #   return(Df)
  # }
  # GeneSetAnalysis_v3<-function(seu = NULL, database = "GSEA", dataset = NULL, ratio = 0.2){
  #   DatabaseList<-list("GSEA"=c("Hall50","BP","MF","CC","Biocarta"),
  #                      "Reactome"=c("Immune",
  #                                   "Metabolism",
  #                                   "Signal_transduction"),
  #                      "GO"=c("BP","MF","CC",
  #                             "immune_system_process",
  #                             "response_to_stimulus",
  #                             "signaling",
  #                             "metabolic_process",
  #                             "regulation_of_vasculature_development",
  #                             "signal_transduction"))
  #   require(Seurat)
  #   if(is.null(seu)){
  #     return(DatabaseList)
  #   }
  #   require(dplyr)
  #   if(database == "GO"){
  #     if(is.null(dataset)){
  #       print(DatabaseList[["GO"]])
  #       return(seu)
  #     }else if(sum(grepl(dataset, DatabaseList[["GO"]]))==0){
  #       print(DatabaseList[["GO"]])
  #       return(seu)
  #     }else if(sum(grepl(dataset, DatabaseList[["GO"]]))>1){
  #       print(DatabaseList[["GO"]] %>% .[grepl(dataset, .)])
  #       return(seu)
  #     }else{
  #       dataset <- DatabaseList[["GO"]] %>% .[grepl(dataset, .)]
  #       if(dataset == "BP"){
  #         genelist_list <- readRDS("~/Documents/scRNA/GO/BP_GO2gene_high_quality.rds")
  #       }else if(dataset == "MF"){
  #         genelist_list <- readRDS("~/Documents/scRNA/GO/MF_GO2gene_high_quality.rds")
  #       }else if(dataset == "CC"){
  #         genelist_list <- readRDS("~/Documents/scRNA/GO/CC_GO2gene_high_quality.rds")
  #       }else if(dataset == "immune_system_process"){
  #         genelist_list <- readRDS("~/Documents/scRNA/GO/GOTerms_BP/GO:0002376_GO2Gene_all.rds")
  #       }else if(dataset == "response_to_stimulus"){
  #         genelist_list <- readRDS("~/Documents/scRNA/GO/GOTerms_BP/GO:0050896_GO2Gene_all.rds")
  #       }else if(dataset == "signaling"){
  #         genelist_list <- readRDS("~/Documents/scRNA/GO/GOTerms_BP/GO:0023052_GO2Gene_all.rds")
  #       }else if(dataset == "metabolic_process"){
  #         genelist_list <- readRDS("~/Documents/scRNA/GO/GOTerms_BP/GO:0008152_GO2Gene_all.rds")
  #       }else if(dataset == "regulation_of_vasculature_development"){
  #         genelist_list <- readRDS("~/Documents/scRNA/GO/GOTerms_BP/GO:1901342_GO2Gene_all.rds")
  #       }else if(dataset == "signal_transduction"){
  #         genelist_list <- readRDS("~/Documents/scRNA/GO/GOTerms_BP/GO:0007165_GO2Gene_all.rds")
  #       }
  #     }
  #   }else if(database == "Reactome"){
  #     if(is.null(dataset)){
  #       print(DatabaseList[["Reactome"]])
  #       return(seu)
  #     }else if(sum(grepl(dataset, DatabaseList[["Reactome"]]))==0){
  #       print(DatabaseList[["Reactome"]])
  #       return(seu)
  #     }else if(sum(grepl(dataset, DatabaseList[["Reactome"]]))>1){
  #       print(DatabaseList[["Reactome"]] %>% .[grepl(dataset, .)])
  #       return(seu)
  #     }else{
  #       dataset <- DatabaseList[["Reactome"]] %>% .[grepl(dataset, .)]
  #       Path2GenesymbolAll <- readRDS("~/Documents/scRNA/Reactome/Path2GenesymbolAll.rds")
  #       if(dataset == "Immune"){
  #         ImmuneSystem <- readRDS("~/Documents/scRNA/Reactome/ImmuneSystem.rds")
  #         genelist_list <- Path2GenesymbolAll[names(Path2GenesymbolAll)[names(Path2GenesymbolAll) %in% ImmuneSystem]]
  #       }else if(dataset == "Metabolism"){
  #         Metabolism <- readRDS("~/Documents/scRNA/Reactome/Metabolism.rds")
  #         genelist_list <- Path2GenesymbolAll[names(Path2GenesymbolAll)[names(Path2GenesymbolAll) %in% Metabolism]]
  #       }else if(dataset == "Signal_transduction"){
  #         SignalTransduction <- readRDS("~/Documents/scRNA/Reactome/SignalTransduction.rds")
  #         genelist_list <- Path2GenesymbolAll[names(Path2GenesymbolAll)[names(Path2GenesymbolAll) %in% SignalTransduction]]
  #       }
  #     }
  #   }else if(database == "GSEA"){
  #     if(is.null(dataset)){
  #       print(DatabaseList[["GSEA"]])
  #       return(seu)
  #     }else if(sum(grepl(dataset, DatabaseList[["GSEA"]]))==0){
  #       print(DatabaseList[["GSEA"]])
  #       return(seu)
  #     }else if(sum(grepl(dataset, DatabaseList[["GSEA"]]))>1){
  #       print(DatabaseList[["GSEA"]] %>% .[grepl(dataset, .)])
  #       return(seu)
  #     }else{
  #       dataset <- DatabaseList[["GSEA"]] %>% .[grepl(dataset, .)]
  #       if(dataset == "Hall50"){
  #         genelist_list <- readRDS("~/Documents/scRNA/HEV project/GSEA/h.all.v6.2.symbols.gmt.rds")
  #       }else if(dataset == "BP"){
  #         genelist_list <- readRDS("~/Documents/scRNA/HEV project/GSEA/c5.bp.v6.2.symbols.gmt.rds")
  #       }else if(dataset == "MF"){
  #         genelist_list <- readRDS("~/Documents/scRNA/HEV project/GSEA/c5.mf.v6.2.symbols.gmt.rds")
  #       }else if(dataset == "CC"){
  #         genelist_list <- readRDS("~/Documents/scRNA/HEV project/GSEA/c5.cc.v6.2.symbols.gmt.rds")
  #       }else if(dataset == "Biocarta"){
  #         genelist_list <- readRDS("~/Documents/scRNA/HEV project/GSEA/c2.cp.biocarta.v6.2.symbols.gmt.txt.rds")
  #       }
  #     }
  #   }
  #   matr <- GetAssayData(seu)
  #   require(AUCell)
  #   if(is.null(seu@misc$AUCell[["cells_rankings"]])){
  #     seu@misc$AUCell<-list()
  #     seu@misc$AUCell[["cells_rankings"]] <- AUCell_buildRankings(matr, plotStats=TRUE)
  #   }
  #   tmp<-genelist_list
  #   for (i in names(tmp)){
  #     if(sum(tmp[[i]] %in% rownames(matr))/length(tmp[[i]]) < (1-ratio)){
  #       tmp[[i]]<-NULL
  #     }
  #   }
  #   seu@misc[["AUCell"]][[database]][[dataset]] <-
  #     getAUC(AUCell_calcAUC(tmp, seu@misc$AUCell[["cells_rankings"]], nCores = 2)) %>%
  #     .[apply(., 1, sum)>0, ]
  #   return(seu)
  # }
  # HeatmapGeneSet_v3<-function(seu = NULL, group = NULL, order_by_group = T,
  #                             score = "tscore", order_method = NULL, n = NULL,
  #                             database = "none", dataset = NULL, change_name = NULL, only_end_term = NULL,
  #                             color_scheme="D"){
  #   DatabaseList<-list("GSEA"=c("Hall50","BP","MF","CC","Biocarta"),
  #                      "Reactome"=c("Immune",
  #                                   "Metabolism",
  #                                   "Signal_transduction"),
  #                      "GO"=c("BP","MF","CC",
  #                             "immune_system_process",
  #                             "response_to_stimulus",
  #                             "signaling",
  #                             "metabolic_process",
  #                             "regulation_of_vasculature_development",
  #                             "signal_transduction"))
  #   require(Seurat)
  #   if(is.null(seu)){
  #     return(DatabaseList)
  #   }
  #   require(dplyr)
  #   for (i in names(DatabaseList)) {
  #     if(database == i){
  #       if(is.null(dataset)){
  #         print(DatabaseList[[i]])
  #         return(seu)
  #       }else if(sum(grepl(dataset, DatabaseList[[i]]))==0){
  #         print(DatabaseList[[i]])
  #         return(seu)
  #       }else if(sum(grepl(dataset, DatabaseList[[i]]))>1){
  #         print(DatabaseList[[i]] %>% .[grepl(dataset, .)])
  #         return(seu)
  #       }else{
  #         dataset <- DatabaseList[[i]] %>% .[grepl(dataset, .)]
  #       }
  #     }
  #   }
  #   require(rlang)
  #   AUCmatr_new <- seu@misc[["AUCell"]][[database]][[dataset]]
  #   meta <- seu@meta.data
  #   if(database %in% c("GO","Reactome")){
  #     change_name <- change_name %||% T
  #     only_end_term <- only_end_term %||% T
  #   }else{
  #     change_name <- change_name %||% F
  #     only_end_term <- only_end_term %||% F
  #   }
  #   if(only_end_term){
  #     if(database == "GO"){
  #       require(ontologyIndex)
  #       GOdata <- readRDS("~/Documents/scRNA/GO/GOdata_obo.rds")
  #       tmp<-vector()
  #       for (i in rownames(AUCmatr_new)) {
  #         if(sum(rownames(AUCmatr_new) %in% GOdata$children[[i]]) == 0) tmp <- c(tmp, i)
  #       }
  #       AUCmatr_new <- AUCmatr_new[tmp,]
  #     }else if(database == "Reactome"){
  #       Path2Genesymbol <- readRDS("~/Documents/scRNA/Reactome/Path2Genesymbol.rds")
  #       AUCmatr_new <- AUCmatr_new[rownames(AUCmatr_new) %in% names(Path2Genesymbol),]
  #     }
  #   }
  #   CalcScoreGeneral<-function(matr, meta, group, method){
  #     scores<-data.frame(row.names = rownames(matr))
  #     if(is.factor(meta[,group])){
  #       cluster<-levels(meta[,group])
  #     }else{
  #       cluster<-unique(meta[,group])
  #     }
  #     if(method=="tscore"){
  #       for (i in cluster) {
  #         for (j in rownames(matr)) {
  #           x<-t.test(matr[j,meta[,group]==i],
  #                     matr[j,meta[,group]!=i])
  #           scores[j,i]<-x$statistic
  #         }
  #       }
  #     }
  #     if(method=="p"){
  #       for (i in cluster) {
  #         for (j in rownames(matr)) {
  #           x<-t.test(matr[j,meta[,group]==i],
  #                     matr[j,meta[,group]!=i])
  #           scores[j,i]<- -log10(x$p.value)
  #         }
  #       }
  #     }
  #     if(method=="mean"){
  #       for (i in cluster) {
  #         for (j in rownames(matr)) {
  #           scores[j,i]<-mean(matr[j,meta[,group]==i])
  #         }
  #       }
  #     }
  #     if(method=="median"){
  #       for (i in cluster) {
  #         for (j in rownames(matr)) {
  #           scores[j,i]<-median(matr[j,meta[,group]==i])
  #         }
  #       }
  #     }
  #     if(method=="zscore"){
  #       clustermean<-as.matrix(CalcScoreGeneral(matr, meta, group, "mean"))
  #       for (i in cluster) {
  #         for (j in rownames(matr)) {
  #           scores[j,i]<-(clustermean[j,i]-mean(clustermean[j,]))/sd(clustermean[j,])
  #         }
  #       }
  #     }
  #     return(scores)
  #   }
  #   if(order_by_group){
  #     ScoreAndOrder<-function(matr, meta, group, score, order_method="none", n=0){
  #       tmp<-as.matrix(CalcScoreGeneral(matr, meta, group, score))
  #       if(order_method=="none") return(tmp)
  #       tmp_max<-vector()
  #       for (i in rownames(tmp)) {
  #         tmp_max[i]<-colnames(tmp)[tmp[i,]==max(tmp[i,])]
  #       }
  #       tmp_cluster<-colnames(tmp) %>% .[. %in% tmp_max]
  #       if(order_method=="value"){
  #         tmp2<-list()
  #         for (i in tmp_cluster) {
  #           tmp2[[i]]<-tmp[tmp_max==i,] %>% .[order(.[,i],decreasing = T),]
  #         }
  #       }
  #       if(order_method=="sd"){
  #         tmp2<-list()
  #         for (i in tmp_cluster) {
  #           tmp2[[i]]<-tmp[tmp_max==i,]
  #           tmp_sd<-apply(tmp2[[i]],1,sd)
  #           tmp2[[i]]<-tmp2[[i]] %>% .[order(tmp_sd,decreasing = T),]
  #         }
  #       }
  #       if(order_method=="p"){
  #         tmp2<-list()
  #         tmp_p<-CalcScoreGeneral(matr, meta, group, "p")
  #         for (i in tmp_cluster) {
  #           tmp2[[i]]<-as.data.frame(tmp)[tmp_max==i,]
  #           tmp2[[i]][,"p"]<-tmp_p[rownames(tmp2[[i]]),i]
  #           tmp2[[i]]<-tmp2[[i]] %>% .[order(.$p,decreasing = T),] %>% .[.$p>2,] %>% .[,-ncol(.)]
  #         }
  #       }
  #       if(n==0){
  #         tmp3<-tmp2[[names(tmp2)[1]]]
  #         for (i in names(tmp2)[2:length(tmp2)]) {
  #           tmp3<-rbind(tmp3, tmp2[[i]])
  #         }
  #       }else{
  #         tmp3<-head(tmp2[[names(tmp2)[1]]],n)
  #         for (i in names(tmp2)[2:length(tmp2)]) {
  #           tmp3<-rbind(tmp3, head(tmp2[[i]],n))
  #         }
  #       }
  #       return(tmp3)
  #     }
  #     n <- n %||% 6
  #     ToPlot <- ScoreAndOrder(AUCmatr_new, meta, group, score, order_method, n)
  #   }else{
  #     order_method <- order_method %||% "sd"
  #     n <- n %||% nrow(AUCmatr_new)
  #     ToPlot<-CalcScoreGeneral(AUCmatr_new, meta, group, score)
  #     if(order_method == "sd"){
  #       ToPlot <- ToPlot %>%
  #         .[order(apply(., 1, sd), decreasing = T)[c(1:n)],]
  #     }
  #   }
  #   if(change_name){
  #     ChangeName<-function(char, ref, from, to){
  #       tmp<-char
  #       for (i in 1:length(char)) {
  #         tmp[i]<-as.character(ref[ref[,from]==char[i],to][1])
  #       }
  #       return(tmp)
  #     }
  #     if(database == "GO"){
  #       GO_Names <- readRDS("~/Documents/scRNA/GO/GO_Names.rds")
  #       rownames(ToPlot) <- ChangeName(rownames(ToPlot), GO_Names, "V1", "V2")
  #     }else if(database == "Reactome"){
  #       Ensembl2ReactomeAll <- readRDS("~/Documents/scRNA/Reactome/Ensembl2ReactomeAll.rds")
  #       rownames(ToPlot) <- ChangeName(rownames(ToPlot), Ensembl2ReactomeAll, "V2","V4")
  #     }
  #   }
  #   heatmap_marker<-function(markers_score, color_scheme="D"){
  #     require(ggplot2)
  #     tmp<-as.data.frame(markers_score)
  #     tmp$id<-rownames(tmp)
  #     tmp$id<-factor(tmp$id, levels=unique(rev(tmp$id)))
  #     tmp<-reshape2::melt(tmp)
  #     color_set<-c("#08306b","#2171b5","#6baed6","#c6dbef","#f7fbff",
  #                  "#fff5f0","#fcbba1","#fb6a4a","#cb181d","#67000d")
  #     if(color_scheme=="b-r-light"){
  #       color_set<-c("#2171b5","#6baed6","#c6dbef","#f7fbff",
  #                    "#fff5f0","#fcbba1","#fb6a4a","#cb181d")
  #     }
  #     if(color_scheme=="g"){
  #       color_set<-c("#f7fcfd","#e5f5f9","#ccece6","#99d8c9",
  #                    "#66c2a4","#41ae76","#238b45","#006d2c","#00441b")
  #     }
  #     if(color_scheme=="g-light"){
  #       color_set<-c("#edf8fb","#b2e2e2","#66c2a4","#2ca25f","#006d2c")
  #     }
  #     if(color_scheme %in% c("A","B","C","D","E")){
  #       require(viridis)
  #       color_set<-viridis_pal(option = color_scheme)(20)
  #     }
  #     p <- ggplot(tmp, aes(variable, id)) +
  #       geom_tile(aes(fill = value), colour = "white") +
  #       scale_fill_gradientn(colors = color_set)+
  #       theme_classic()+
  #       labs(x = "", y = "")+
  #       scale_y_discrete(position = "right")+
  #       theme(axis.text.x=element_text(angle = 45,hjust = 1))
  #     return(p)
  #   }
  #   p<-heatmap_marker(ToPlot, color_scheme) + labs(fill=score)
  #   return(p)
  # }
  # run_GSVA <- function(Obj, Gmtfilename){
  #   require("biomaRt")
  #   human = useMart("ensembl", dataset = "hsapiens_gene_ensembl")
  #   mouse = useMart("ensembl", dataset = "mmusculus_gene_ensembl")
  #   genelists = getLDS(attributes = c("mgi_symbol"), filters = "mgi_symbol", values = Obj@var.genes , mart = mouse,
  #                      attributesL = c("hgnc_symbol"), martL = human, uniqueRows=T)
  #   genes<-genelists$HGNC.symbol %>% .[!duplicated(.)]
  #   cells<-rownames(Obj@meta.data)
  #   new_matrix<-matrix(nrow = length(genes), ncol = length(cells))
  #   rownames(new_matrix)<-genes
  #   colnames(new_matrix)<-cells
  #   mouse_matrix<-as.matrix(Obj@data)
  #   for (i in genes){
  #     gene_mouse<-genelists$MGI.symbol[genelists$HGNC.symbol==i][1]
  #     new_matrix[i,]<-mouse_matrix[gene_mouse,]
  #   }
  #   require("GSEABase")
  #   Gmtfile<-getGmt(paste0("GSEA/",Gmtfilename))
  #   require("GSVA")
  #   es <- gsva(new_matrix, Gmtfile, min.sz=10, max.sz=500, verbose=TRUE, parallel.sz=1)
  #   return(es)
  # }
}
