#' @title FUNCTION_TITLE
#' @description FUNCTION_DESCRIPTION
#' @param matr PARAM_DESCRIPTION
#' @param f PARAM_DESCRIPTION
#' @param f2 PARAM_DESCRIPTION, Default: NULL
#' @param features PARAM_DESCRIPTION, Default: NULL
#' @param ncol PARAM_DESCRIPTION, Default: 1
#' @param lab_fill PARAM_DESCRIPTION, Default: 'group'
#' @param scales PARAM_DESCRIPTION, Default: 'free_y'
#' @param violin PARAM_DESCRIPTION, Default: T
#' @param box PARAM_DESCRIPTION, Default: T
#' @param width PARAM_DESCRIPTION, Default: 0.9
#' @param pt PARAM_DESCRIPTION, Default: F
#' @param pt.style PARAM_DESCRIPTION, Default: c("quasirandom", "jitter")
#' @param pt.size PARAM_DESCRIPTION, Default: 1.5
#' @param pt.alpha PARAM_DESCRIPTION, Default: 0.35
#' @return OUTPUT_DESCRIPTION
#' @details DETAILS
#' @examples
#' \dontrun{
#' if(interactive()){
#'  #EXAMPLE1
#'  }
#' }
#' @rdname StackedViolin
#' @export

StackedViolin <- function(matr, f, f2 = NULL, features = NULL, ncol = 1, lab_fill = "group",
                          scales = "free_y", violin = T, box = T, width = 0.9,
                          pt = F, pt.style = c("quasirandom", "jitter"), pt.size = 1.5, pt.alpha = 0.35){
  library(ggplot2)
  library(rlang)
  library(dplyr)
  library(reshape2)
  features <- features %||% rownames(matr)
  if(!is_empty(setdiff(features, rownames(matr)))){
    message(paste0(setdiff(features, rownames(matr)), collapse = ", "), " not found")
    features <- intersect(features, rownames(matr))
  }
  f2 <- f2 %||% data.frame(row.names = colnames(matr))
  ToPlot <-
    cbind(f, f2, as.data.frame(t(matr[features,,drop = F]))) %>%
    melt(measure.vars = features)
  if(!violin & !box & !pt) stop("No plot type defined (violin/boxplot/point)")
  x <- ifelse(is_empty(f2), "f", "f2")
  p <- ggplot(ToPlot, aes_string(x = x, y = "value", fill = x))
  if(violin) {
    p <- p +
      geom_violin(scale = "width", width = width)
  }
  if(box & !violin) {
    outlier.size <- ifelse(pt, 0, pt.size)
    p <- p +
      geom_boxplot(outlier.size = outlier.size, width = width)
  }
  if(pt) {
    pt.style <- pt.style[1]
    if(!pt.style %in% c("quasirandom", "jitter")) stop('"pt.style" shoule be "quasirandom" or "jitter"')
    if(pt.style == "jitter") p <- p + geom_jitter(width = width/2, size = pt.size, alpha= pt.alpha)
    if(pt.style == "quasirandom") {
      import("ggbeeswarm")
      p <- p + geom_quasirandom(size = pt.size, width = width/2, alpha= pt.alpha)
    }
  }
  if(box & violin) {
    outlier.size <- ifelse(pt, 0, pt.size)
    p <- p +
      geom_boxplot(fill = "white", outlier.size = outlier.size, width = 0.1, outlier.alpha = pt.alpha)
  }

  if(is_empty(f2)){
    p <- p +
      facet_wrap( ~variable, ncol = ncol, strip.position="left", scales = scales)+
      ylab(NULL) +
      xlab(NULL) +
      theme_classic() +
      theme(strip.background = element_blank(),
            strip.placement = "outside",
            legend.position = "none",
            axis.text.x=element_text(angle = 45,hjust = 1)) +
      labs(fill = lab_fill) +
      scale_y_continuous(expand = expansion(mult = c(0,0.05)))
  }else{
    p <- p +
      facet_grid(vars(variable), vars(f), switch = c("both"), scales = scales)+
      ylab(NULL) +
      xlab(NULL) +
      theme_classic() +
      theme(strip.background = element_blank(),
            strip.placement = "outside",
            axis.text.x = element_blank()) +
      labs(fill = lab_fill) +
      scale_y_continuous(expand = expansion(mult = c(0,0.05)))
  }

  return(p)
}

#' @title FUNCTION_TITLE
#' @description FUNCTION_DESCRIPTION
#' @param seu PARAM_DESCRIPTION
#' @param features PARAM_DESCRIPTION
#' @param group.by PARAM_DESCRIPTION, Default: 'seurat_clusters'
#' @param split.by PARAM_DESCRIPTION, Default: NULL
#' @param cell PARAM_DESCRIPTION, Default: NULL
#' @param ... PARAM_DESCRIPTION
#' @return OUTPUT_DESCRIPTION
#' @details DETAILS
#' @examples
#' \dontrun{
#' if(interactive()){
#'  #EXAMPLE1
#'  }
#' }
#' @rdname StackedViolin_v3
#' @export

StackedViolin_v3 <-
  function(seu, features, group.by = "seurat_clusters", split.by = NULL, cell = NULL, ...){
  require(rlang)
  cell <- cell %||% colnames(seu)
  matr <- t(FetchData(seu, vars = features, cells = cell))
  f <- factor(seu[[group.by]][cell,])
  f2 <- seu[[split.by]][cell,]
  p <- StackedViolin(matr, f, f2, features, ...)
  return(p)
}

# seu <- readRDS("~/R documents/2020-2-10 EC PyMT and E0771/rds/PyMTEC_old.rds")
# matr <- as.matrix(GetAssayData(seu))[1:20,]
# f <- seu@meta.data$cluster
# f2 <- NULL
# features = c("Selp","Sele","Vwf","Glycam1")
# cell = NULL
# split.by = "orig.ident"
# group.by = "seurat_clusters"
# ncol = 1
# lab_fill = "group"
# StackedViolin(matr, f, f2, ncol = 2)
# StackedViolin_v3(seu = seu, features = features, split.by = split.by)
# StackedViolin_v3(pbmc, features = c("CD3D"), violin = T, box = T, pt = F)
