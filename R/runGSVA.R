#' @title FUNCTION_TITLE
#' @description FUNCTION_DESCRIPTION
#' @param matr PARAM_DESCRIPTION
#' @param genesets PARAM_DESCRIPTION
#' @param min.sz PARAM_DESCRIPTION, Default: 0
#' @param max.sz PARAM_DESCRIPTION, Default: Inf
#' @return OUTPUT_DESCRIPTION
#' @details DETAILS
#' @examples
#' \dontrun{
#' if(interactive()){
#'  #EXAMPLE1
#'  }
#' }
#' @rdname RunGSVA
#' @export

RunGSVA <- function(matr, genesets, min.sz = 0, max.sz = Inf) {
  library("GSVA")
  es <- gsva(matr, genesets, min.sz = min.sz, max.sz = max.sz, verbose=TRUE)
  return(es)
}

#' @title FUNCTION_TITLE
#' @description FUNCTION_DESCRIPTION
#' @param seu PARAM_DESCRIPTION
#' @param group.by PARAM_DESCRIPTION
#' @param genesets PARAM_DESCRIPTION
#' @param min.sz PARAM_DESCRIPTION, Default: 0
#' @param max.sz PARAM_DESCRIPTION, Default: Inf
#' @param title PARAM_DESCRIPTION, Default: 'genesets'
#' @return OUTPUT_DESCRIPTION
#' @details DETAILS
#' @examples
#' \dontrun{
#' if(interactive()){
#'  #EXAMPLE1
#'  }
#' }
#' @rdname RunGSVA_v3
#' @export

RunGSVA_v3 <- function(seu, group.by, genesets, min.sz = 0, max.sz = Inf, title = "genesets") {
  library(rlist)
  library(Seurat)
  f <- factor(seu@meta.data[, group.by])
  matr <-
    GetAssayData(seu, assay = "RNA") %>%
    expm1 %>%
    t() %>%
    as.data.frame() %>%
    split(f) %>%
    lapply(function(x) apply(x,2,mean)) %>%
    list.cbind() %>%
    log1p()
  es <- RunGSVA(matr, genesets, min.sz = min.sz, max.sz = max.sz)
  seu@misc[["GSVA"]][[title]] <- es
  return(seu)
}

#' @title FUNCTION_TITLE
#' @description FUNCTION_DESCRIPTION
#' @param seu PARAM_DESCRIPTION
#' @param group.by PARAM_DESCRIPTION
#' @param geneset PARAM_DESCRIPTION
#' @param ident.1 PARAM_DESCRIPTION, Default: NULL
#' @param ident.2 PARAM_DESCRIPTION, Default: NULL
#' @param logFC.cutoff PARAM_DESCRIPTION, Default: c(-1, 1)
#' @param title PARAM_DESCRIPTION, Default: deparse(substitute(geneset))
#' @return OUTPUT_DESCRIPTION
#' @details DETAILS
#' @examples
#' \dontrun{
#' if(interactive()){
#'  #EXAMPLE1
#'  }
#' }
#' @rdname GSEAplot
#' @export

GSEAplot <- function(seu, group.by, geneset, ident.1 = NULL, ident.2 = NULL,
                     logFC.cutoff = c(-1,1), title = deparse(substitute(geneset))) {
  library(ggplot2)
  library(cowplot)
  library(Seurat)
  library(rlang)
  library(scales)
  library(dplyr)

  f <- factor(seu@meta.data[, group.by])
  ident.1 <- ident.1 %||% levels(f)[1]
  cell.1 <- (f == ident.1)
  cell.2 <- if(is.null(ident.2)) f != ident.1 else f == ident.2
  title <- paste0("Enrichment Plot: ", title)
  title2 <- paste0("Rank in gene list (", ident.1, " vs. ", ident.2 %||% paste0("non ", ident.1), ")")

  matr <-
    GetAssayData(seu) %>%
    expm1() %>%
    .[apply(., 1, function(x) sum(x) > 0),]
  geneset <- intersect(geneset, rownames(matr))
  scores <-
    GetAssayData(seu) %>%
    apply(1, function(x) log(mean(x[cell.1]+1)/mean(x[cell.2]+1))) %>%
    .[order(., decreasing = T)] %>%
    as.data.frame() %>%
    `colnames<-`("LogFC") %>%
    mutate(genes = factor(rownames(.), levels = rownames(.)),
           geneset = as.numeric(rownames(.) %in% geneset),
           rank = 1:nrow(.),
           centered_rank = mean(which.min(abs(.$LogFC))) - (1:nrow(.)))

  step.T <- 1/sum(scores$geneset)
  step.F <- 1/sum(!scores$geneset)
  scores$`Enrichment Score (ES)` <-
    1:nrow(scores) %>%
    sapply(function(x) sum(scores$geneset[1:x]) * step.T - sum(!scores$geneset[1:x]) * step.F)

  p1 <- ggplot(scores, aes(x = rank, y = `Enrichment Score (ES)`)) +
    geom_step() +
    theme_bw() +
    theme(axis.text.x = element_blank(),
          axis.ticks.x = element_blank(),
          axis.title.x = element_blank(),
          plot.margin = margin(1,1,0,1)) +
    scale_x_continuous(expand = c(0, 0), limits = c(0, nrow(scores))) +
    ggtitle(title)

  p2 <- ggplot(scores) +
    geom_bar(mapping = aes(x = rank, y = 1, fill = centered_rank, alpha = 0.8),
             stat = "identity", width = 1) +
    geom_bar(mapping = aes(x = rank, y = geneset), stat = "identity", width = 10) +
    theme_bw() +
    theme(axis.text = element_blank(),
          axis.ticks = element_blank(),
          panel.grid = element_blank(),
          axis.title = element_blank(),
          legend.position = "none",
          plot.margin = margin(0,1,0,1)) +
    scale_fill_gradient2(low = muted("blue"), high = muted("red")) +
    scale_x_continuous(expand = c(0, 0), limits = c(0, nrow(scores)))

  p3 <- ggplot(scores, aes(x = rank, y = LogFC)) +
    geom_bar(stat = "identity", width = 1) +
    theme_bw() +
    theme(axis.text.x=element_blank(),
          axis.ticks.x=element_blank(),
          panel.grid.major.x = element_blank(),
          plot.margin = margin(0,1,1,1)) +
    coord_cartesian(ylim = logFC.cutoff) +
    scale_x_continuous(expand = c(0, 0), limits = c(0, nrow(scores))) +
    xlab(title2)

  p <- plot_grid(p1, p2, p3, ncol = 1, align = "v", rel_heights = c(6,1,3))

  return(p)
}



# seu <- readRDS("~/R documents/2020-2-10 EC PyMT and E0771/rds/PyMTEC_old.rds")
# genesets <- Genesets_data$mouse$GSEA_mouse_gene_transformed$`hallmark gene sets`
# geneset <- genesets$HALLMARK_INTERFERON_GAMMA_RESPONSE
# seu <- RunGSVA_v3(seu, group.by = "cluster", genesets = genesets)
# seu@misc$GSVA$genesets %>% Heatmap()

# GSEAplot(seu, "cluster", geneset = geneset, ident.1 = "TU_HEV", ident.2 = "LN_HEV")
