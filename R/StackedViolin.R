#' @title FUNCTION_TITLE
#' @description FUNCTION_DESCRIPTION
#' @param matr PARAM_DESCRIPTION
#' @param f PARAM_DESCRIPTION
#' @param f2 PARAM_DESCRIPTION, Default: NULL
#' @param features PARAM_DESCRIPTION, Default: NULL
#' @param ncol PARAM_DESCRIPTION, Default: 1
#' @param lab_fill PARAM_DESCRIPTION, Default: 'group'
#' @param scales PARAM_DESCRIPTION, Default: 'free_y'
#' @param type PARAM_DESCRIPTION, Default: c("violin", "boxplot")
#' @param outlier.size PARAM_DESCRIPTION, Default: 1
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
                          scales = "free_y", type = c("violin","boxplot"), outlier.size = 1){
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
    cbind(f, f2, as.data.frame(t(matr[features,]))) %>%
    melt(measure.vars = features)
  type <- type[1]
  if(type == "violin") {
    p <- ggplot(ToPlot, aes(x=f,y=value, fill=f)) +
      geom_violin(scale = "width")
  } else if(type == "boxplot") {
    p <- ggplot(ToPlot, aes(x=f,y=value, fill=f)) +
      geom_boxplot(outlier.size = outlier.size)
  } else {
    stop('"type" argument should be "violin" or "boxplot"')
  }
  if(is_empty(f2)){
    p <- p +
      facet_wrap( ~variable, ncol = ncol, strip.position="left", scales = scales)+
      ylab(NULL) +
      xlab(NULL) +
      theme(strip.background = element_blank(),
            strip.placement = "outside",
            legend.position = "none",
            axis.text.x=element_text(angle = 45,hjust = 1)) +
      theme_classic() +
      labs(fill = lab_fill)
  }else{
    p <- p +
      facet_grid(vars(variable), vars(f2), switch = c("both"), scales = scales)+
      ylab(NULL) +
      xlab(NULL) +
      theme(strip.background = element_blank(),
            strip.placement = "outside",
            legend.position = "none",
            axis.text.x=element_text(angle = 45,hjust = 1)) +
      theme_classic() +
      labs(fill = lab_fill)
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
#' @param ncol PARAM_DESCRIPTION, Default: 1
#' @param scales PARAM_DESCRIPTION, Default: 'free_y'
#' @param type PARAM_DESCRIPTION, Default: c("violin", "boxplot")
#' @param outlier.size PARAM_DESCRIPTION, Default: 1
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
  function(seu, features, group.by = "seurat_clusters", split.by = NULL, cell = NULL,
           ncol = 1, scales = "free_y", type = c("violin","boxplot"), outlier.size = 1){
  require(rlang)
  cell <- cell %||% colnames(seu)
  matr <- t(FetchData(seu, vars = features, cells = cell))
  f <- factor(seu[[group.by]][cell,])
  f2 <- seu[[split.by]][cell,]
  p <- StackedViolin(matr, f, f2, features, ncol, lab_fill = group.by,
                     scales = scales, type = type, outlier.size = outlier.size)
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
