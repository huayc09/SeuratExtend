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
#' ClusterDistrBar(origin = pbmc$orig.ident, cluster = Idents(pbmc))
#'
#' # absolute cell count
#' ClusterDistrBar(origin = pbmc$orig.ident, cluster = Idents(pbmc), percent = F)
#'
#' # reverse x and y axis, normalized by sample size
#' ClusterDistrBar(origin = pbmc$orig.ident, cluster = Idents(pbmc), rev = T, normalize = T)
#'
#' # reverse x and y axis, not normalized by sample size
#' ClusterDistrBar(origin = pbmc$orig.ident, cluster = Idents(pbmc), rev = T, normalize = F)
#'
#' # vertical bar plot
#' ClusterDistrBar(origin = pbmc$orig.ident, cluster = Idents(pbmc), flip = F)
#'
#' # export matrix
#' ClusterDistrBar(origin = pbmc$orig.ident, cluster = Idents(pbmc), plot = F)
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

  if(flip) ToPlot$Var2 <- ToPlot$Var2 %>% factor(levels = rev(levels(.)))
  p <-
    ggplot(ToPlot, aes(x = Var2, y = value, fill = Var1)) +
    geom_bar(stat="identity", position = position_stack(reverse = TRUE)) +
    theme_classic() +
    labs(x = x.label, y = y.label, fill = fill.label) +
    if(flip) coord_flip()
  return(p)
}

#' @title FUNCTION_TITLE
#' @description FUNCTION_DESCRIPTION
#' @param meta PARAM_DESCRIPTION
#' @param cluster PARAM_DESCRIPTION
#' @param origin PARAM_DESCRIPTION
#' @param method PARAM_DESCRIPTION, Default: 'value'
#' @param do.heatmap PARAM_DESCRIPTION, Default: T
#' @return OUTPUT_DESCRIPTION
#' @details DETAILS
#' @examples
#' \dontrun{
#' if(interactive()){
#'  #EXAMPLE1
#'  }
#' }
#' @rdname BarOfCluster
#' @export

BarOfCluster <- function(meta,cluster,origin,method="value",do.heatmap=T){
  warning("BarOfCluster() is deprecated. Please use ClusterDistrBar() instead.")
  ClusterDistrBar(origin = meta[,origin], cluster = meta[,cluster], percent = method=="percent", plot = do.heatmap)
}

# meta <- seu@meta.data
# origin <- meta$orig.ident
# cluster <- meta$cluster
# ClusterDistrBar(origin, cluster, rev = T)
# BarOfCluster(meta, "cluster", "orig.ident", method = "percent", do.heatmap = T)
