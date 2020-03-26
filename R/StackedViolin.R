#' @title FUNCTION_TITLE
#' @description FUNCTION_DESCRIPTION
#' @param matr PARAM_DESCRIPTION
#' @param f PARAM_DESCRIPTION
#' @param features PARAM_DESCRIPTION, Default: NULL
#' @param ncol PARAM_DESCRIPTION, Default: 1
#' @param lab_fill PARAM_DESCRIPTION, Default: 'group'
#' @return OUTPUT_DESCRIPTION
#' @details DETAILS
#' @examples
#' \dontrun{
#' if(interactive()){
#'  #EXAMPLE1
#'  }
#' }
#' @seealso
#'  \code{\link[reshape2]{melt}}
#' @rdname StackedViolin
#' @export
#' @importFrom reshape2 melt

StackedViolin <- function(matr, f, features = NULL, ncol = 1, lab_fill = "group"){
  require(ggplot2)
  require(rlang)
  require(dplyr)
  features <- features %||% rownames(matr)
  if(!is_empty(setdiff(features, rownames(matr)))){
    message(paste0(setdiff(features, rownames(matr)), collapse = ", "), " not found")
    features <- intersect(features, rownames(matr))
  }
  ToPlot <-
    cbind(f, as.data.frame(t(matr[features,]))) %>%
    reshape2::melt()
  p<-ggplot(ToPlot, aes(x=f,y=value, fill=f))+
    geom_violin(scale = "width")+
    facet_wrap( ~variable, ncol = ncol, strip.position="left")+
    ylab(NULL) +
    xlab(NULL) +
    theme(strip.background = element_blank(),
          strip.placement = "outside",
          legend.position = "none",
          axis.text.x=element_text(angle = 45,hjust = 1)) +
    theme_classic() +
    labs(fill = lab_fill)
  return(p)
}

#' @title FUNCTION_TITLE
#' @description FUNCTION_DESCRIPTION
#' @param seu PARAM_DESCRIPTION
#' @param features PARAM_DESCRIPTION
#' @param group.by PARAM_DESCRIPTION, Default: 'seurat_clusters'
#' @param cell PARAM_DESCRIPTION, Default: NULL
#' @param ncol PARAM_DESCRIPTION, Default: 1
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

StackedViolin_v3 <- function(seu, features, group.by = "seurat_clusters", cell = NULL, ncol = 1){
  require(rlang)
  cell <- cell %||% colnames(seu)
  matr <- GetAssayData(seu)[features, cell]
  f <- factor(seu[[group.by]][cell,])
  p <- StackedViolin(matr, f, features, ncol, lab_fill = group.by)
  return(p)
}

# matr <- as.matrix(GetAssayData(seu))[1:20,]
# f <- seu@meta.data$cluster
# features = NULL
# StackedViolin(matr, f, ncol = 2)
# StackedViolin_v3(seu = seu, features = features)
