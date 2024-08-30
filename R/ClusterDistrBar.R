#' @title Cluster proportion bar plot
#' @description Plot the percentage/absolute cell count of each cluster in each sample
#' @param origin factor/vector of sample
#' @param cluster factor/vector of cluster
#' @param rev If TRUE, plot the proportion of sample in each cluster, Default: F
#' @param normalize Normalize sample size to 100, Default: rev
#' @param percent If FALSE, plot absolute cell number instead of percentage, Default: T
#' @param plot Generate plot (TRUE) or matrix (FALSE), Default: T
#' @param flip If plot bars horizontally, Default: T
#' @param reverse_order If TRUE, will reverse the default stacking order. Default: T
#' @param width Width of bars, Default: 0.9
#' @param stack If TRUE, plot stacked bars, Default: T
#' @param cols Colors to use for plotting. Default: "pro_default"
#' @param border Border color, Default: "white"
#' @return ggplot object or matrix
#' @details See example
#' @examples
#' ClusterDistrBar(origin = pbmc$orig.ident, cluster = pbmc$cluster)
#'
#' # absolute cell count
#' ClusterDistrBar(origin = pbmc$orig.ident, cluster = pbmc$cluster, percent = F)
#'
#' # reverse x and y axis, normalized by sample size
#' ClusterDistrBar(origin = pbmc$orig.ident, cluster = pbmc$cluster, rev = T, normalize = T)
#'
#' # reverse x and y axis, not normalized by sample size
#' ClusterDistrBar(origin = pbmc$orig.ident, cluster = pbmc$cluster, rev = T, normalize = F)
#'
#' # vertical bar plot
#' ClusterDistrBar(origin = pbmc$orig.ident, cluster = pbmc$cluster, flip = F)
#'
#' # not stacking bars
#' ClusterDistrBar(origin = pbmc$orig.ident, cluster = pbmc$cluster, flip = FALSE, stack = FALSE)
#'
#' # export matrix
#' ClusterDistrBar(origin = pbmc$orig.ident, cluster = pbmc$cluster, plot = F)
#' @rdname ClusterDistrBar
#' @export

ClusterDistrBar <- function(
    origin,
    cluster,
    rev = F,
    normalize = rev,
    percent = T,
    plot = T,
    flip = T,
    reverse_order = T,
    width = 0.9,
    stack = T,
    cols = "pro_default",
    border = "white"
    ){
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
  if(stack) position <- position_stack(reverse = reverse_order) else position <- "dodge"
  p <-
    ggplot(ToPlot, aes(x = Var2, y = value, fill = factor(Var1))) +
    geom_bar(stat="identity", position = position, width = width, color = border) +
    theme_classic() +
    labs(x = x.label, y = y.label, fill = fill.label) +
    scale_y_continuous(expand = c(0, 0)) +
    if(flip) coord_flip()
  p <- p + scale_fill_disc_auto(cols, nlevels(factor(ToPlot$Var1)))
  return(p)
}

#' @title Cluster proportion bar plot
#' @description deprecated function of \code{\link[SeuratExtend:ClusterDistrBar]{ClusterDistrBar()}}
#' @param meta meta data data.frame
#' @param cluster variable name of cluster
#' @param origin variable name of sample
#' @param method 'percent'or 'value'
#' @param do.heatmap Generate plot (TRUE) or matrix (FALSE), Default: T
#' @return ggplot object or matrix
#' @seealso \code{\link[SeuratExtend:ClusterDistrBar]{ClusterDistrBar()}}
#' @rdname BarOfCluster
#' @export

BarOfCluster <- function(meta,cluster,origin,method="value",do.heatmap=T){
  warning("BarOfCluster() is deprecated. Please use ClusterDistrBar() instead.")
  ClusterDistrBar(origin = meta[,origin], cluster = meta[,cluster], percent = method=="percent", plot = do.heatmap)
}

