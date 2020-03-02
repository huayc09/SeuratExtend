#' @title FUNCTION_TITLE
#' @description FUNCTION_DESCRIPTION
#' @param origin PARAM_DESCRIPTION
#' @param cluster PARAM_DESCRIPTION
#' @param rev PARAM_DESCRIPTION, Default: F
#' @param normalize PARAM_DESCRIPTION, Default: rev
#' @param percent PARAM_DESCRIPTION, Default: T
#' @param plot PARAM_DESCRIPTION, Default: T
#' @param flip PARAM_DESCRIPTION, Default: T
#' @return OUTPUT_DESCRIPTION
#' @details DETAILS
#' @examples
#' \dontrun{
#' if(interactive()){
#'  #EXAMPLE1
#'  }
#' }
#' @rdname ClusterDistrBar
#' @export
ClusterDistrBar <- function(origin, cluster, rev = F, normalize = rev, percent = T,
                            plot = T, flip = T){
  library(dplyr)
  library(rlist)
  library(ggplot2)
  library(reshape2)

  origin <- factor(origin)
  cluster <- factor(cluster)

  ToPlot <-
    split(cluster, origin) %>%
    lapply(table) %>%
    list.cbind()
  if(normalize) ToPlot <- apply(ToPlot, 2, function(x) x * 100/sum(x))
  if(rev) ToPlot <- t(ToPlot)
  if(percent) ToPlot <- apply(ToPlot, 2, function(x) x * 100/sum(x))
  if(!plot) return(ToPlot)

  ToPlot <- melt(ToPlot)
  x.label <- ifelse(rev, "Cluster", "Origin")
  fill.label <- ifelse(rev, "Origin", "Cluster")
  y.label <- ifelse(normalize, "Normalized cell counts", "Cell counts")
  y.label <- ifelse(percent, paste("Percentage of", tolower(y.label)), y.label)

  p <-
    ggplot(ToPlot, aes(x = Var2, y = value, fill = Var1)) +
    geom_bar(stat="identity", position = position_stack(reverse = TRUE)) +
    theme_classic() +
    labs(x = x.label, y = y.label, fill = fill.label) +
    if(flip) coord_flip()
  return(p)
}

BarOfCluster <- function(meta,cluster,origin,method="value",do.heatmap=T){
  warning("BarOfCluster() is deprecated. Please use ClusterDistrBar() instead.")
  ClusterDistrBar(origin = meta[,origin], cluster = meta[,cluster], percent = method=="percent", plot = do.heatmap)
}

# seuratObj <- readRDS("~/Documents/scRNA/HEV-seurat3/Robjects/seuratObj.rds")
# meta <- seuratObj@meta.data
# origin <- meta$orig.ident
# cluster <- meta$cluster
# ClusterDistrBar(origin, cluster, plot = F)
# BarOfCluster(meta, "cluster", "orig.ident", method = "percent", do.heatmap = T)
