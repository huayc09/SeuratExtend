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
#' @param cols Colors to use for plotting. Default: "light"
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
#' # vertical bar plot, and keep the default stacking order
#' ClusterDistrBar(origin = pbmc$orig.ident, cluster = pbmc$cluster, flip = F, reverse_order = F)
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
    cols = "light",
    border = "white"
    ){
  library(dplyr)
  library(rlist)
  library(ggplot2)
  library(reshape2)
  
  # Check if origin and cluster have the same length
  if(length(origin) != length(cluster)) {
    stop("'origin' and 'cluster' must have the same length. Please check if they are from the same Seurat object or data frame.")
  }

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

  levels_var1 <- rownames(ToPlot)
  levels_var2 <- colnames(ToPlot)
  ToPlot <- melt(ToPlot)
  x.label <- ifelse(rev, "Cluster", "Origin")
  fill.label <- ifelse(rev, "Origin", "Cluster")
  y.label <- ifelse(normalize, "Normalized cell counts", "Cell counts")
  y.label <- ifelse(percent, paste("Percentage of", tolower(y.label)), y.label)

  ToPlot$Var1 <- factor(ToPlot$Var1, levels = levels_var1)
  if(flip) ToPlot$Var2 <- factor(ToPlot$Var2, levels = rev(levels_var2))
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

#' @title Cluster proportion boxplot by condition
#' @description Plot the percentage/absolute cell count of each cluster grouped by condition
#' @param origin factor/vector of sample
#' @param cluster factor/vector of cluster
#' @param condition factor/vector of condition, should have the same number of unique values as origin. If provided, will create boxplots instead of bar plots, Default: NULL
#' @param rev If TRUE, plot the proportion of sample in each cluster, Default: F
#' @param normalize Normalize sample size to 100, Default: rev
#' @param percent If FALSE, plot absolute cell number instead of percentage, Default: T
#' @param plot Generate plot (TRUE) or matrix (FALSE), Default: T
#' @param flip If plot bars horizontally (only used when condition is NULL), Default: T
#' @param reverse_order If TRUE, will reverse the default stacking order (only used when condition is NULL), Default: T
#' @param stack If TRUE, plot stacked bars (only used when condition is NULL), Default: T
#' @param cols Colors to use for plotting. Default: "light"
#' @param border Border color for bar plots (only used when condition is NULL), Default: "white"
#' @param violin Indicates whether to generate a violin plot (only used when condition is not NULL), Default: FALSE
#' @param box Indicates whether to depict a box plot (only used when condition is not NULL), Default: TRUE
#' @param width Width parameter. When condition is NULL, controls bar width; when condition is provided, controls box plot width. Default: 0.9
#' @param pt Indicates if points should be plotted (only used when condition is not NULL), Default: TRUE
#' @param pt.style Position adjustment. Default choices: "quasirandom", "jitter" (only used when condition is not NULL)
#' @param pt.size Point size setting (only used when condition is not NULL), Default: 1
#' @param pt.alpha Adjusts the transparency of points (only used when condition is not NULL), Default: 1
#' @param style Plot style: "fill" or "outline" (only used when condition is not NULL). Default: "outline"
#' @param ncol Number of columns for multi-panel plots. Default: NULL
#' @param nrow Number of rows for multi-panel plots. Default: NULL
#' @param stat.method Determines if pairwise statistics are added to the plot. Either "wilcox.test" or "t.test", Default: "wilcox.test"
#' @param p.adjust.method Method for adjusting p-values, Default: "none"
#' @param label Specifies label type. Options include "p.signif", "p", "p.adj", "p.format", Default: "p.signif"
#' @param hide.ns If TRUE, the 'ns' symbol is concealed when displaying significance levels, Default: TRUE
#' @param ... Additional arguments passed to VlnPlot2
#' @return ggplot object or matrix
#' @details This function extends ClusterDistrBar to allow grouping by condition, creating boxplots instead of bar plots.
#' When condition is NULL, it behaves exactly like ClusterDistrBar, supporting all its parameters
#' and returning a stacked bar plot. When condition is provided, it transforms the distribution data
#' into a boxplot format using VlnPlot2, with each box representing the distribution of a cluster's
#' proportion across samples within the same condition. This is particularly useful for comparing cluster
#' distributions between experimental conditions.
#'
#' All styling parameters from VlnPlot2 are supported, including violin/box display options, point
#' styles, colors, and statistical testing. For detailed formatting options, refer to the VlnPlot2
#' documentation. Common customizations include:
#'
#' - Changing point display: pt=TRUE/FALSE, pt.size, pt.alpha, pt.style
#'
#' - Adjusting plot style: style="fill"/"outline", violin=TRUE/FALSE
#'
#' - Adding statistical tests: stat.method="wilcox.test"/"t.test", p.adjust.method
#'
#' - Controlling layout: ncol, nrow for multi-cluster displays
#'
#' @examples
#' # Example with condition
#' ClusterDistrPlot(
#'   origin = pbmc$sample_id,
#'   cluster = pbmc$cluster,
#'   condition = pbmc$condition
#' )
#'
#' @rdname ClusterDistrPlot
#' @export

ClusterDistrPlot <- function(
    origin,
    cluster,
    condition = NULL,
    rev = F,
    normalize = rev,
    percent = T,
    plot = T,
    flip = T,
    reverse_order = T,
    stack = T,
    cols = "light",
    border = "white",
    violin = F,
    box = T,
    width = 0.9,
    pt = T,
    pt.style = "quasirandom",
    pt.size = 1,
    pt.alpha = 1,
    style = "outline",
    ncol = NULL,
    nrow = NULL,
    stat.method = "wilcox.test",
    p.adjust.method = "none",
    label = "p.signif",
    hide.ns = TRUE,
    ...
    ){
  library(dplyr)
  library(rlist)
  library(ggplot2)
  library(reshape2)
  
  # Check if origin and cluster have the same length
  if(length(origin) != length(cluster)) {
    stop("'origin' and 'cluster' must have the same length. Please check if they are from the same Seurat object or data frame.")
  }

  # Check for invalid parameter combinations
  if(!is.null(condition) && rev) {
    stop("Cannot use rev=TRUE with condition parameter. When condition is provided, rev must be FALSE.")
  }

  # Get distribution matrix using standard ClusterDistrBar functionality
  distr_matrix <- ClusterDistrBar(
    origin = origin,
    cluster = cluster,
    rev = rev,
    normalize = normalize,
    percent = percent,
    plot = FALSE
  )

  # If condition is not provided or plot=FALSE, return the distribution matrix
  if(is.null(condition) || !plot) {
    if(!plot) return(distr_matrix)
    # If condition is NULL but plot=TRUE, use standard ClusterDistrBar plot
    return(ClusterDistrBar(
      origin = origin,
      cluster = cluster,
      rev = rev,
      normalize = normalize,
      percent = percent,
      plot = TRUE,
      flip = flip,
      reverse_order = reverse_order,
      width = width,
      stack = stack,
      cols = cols,
      border = border
    ))
  }

  # Process condition parameter
  if(length(condition) == length(origin)) {
    # Create condition mapping table from the full data
    cond_map <- data.frame(origin = origin, condition = condition)
    cond_map <- unique(cond_map)
  } else {
    # If condition is provided as a mapping table or doesn't match origin length
    stop("The condition parameter must be the same length as origin or a valid mapping table")
  }

  # Ensure each origin is mapped to exactly one condition
  if(any(table(cond_map$origin) > 1)) {
    stop("Each origin (sample) must be associated with exactly one condition")
  }

  # Create condition factor for the columns of distr_matrix
  col_origins <- colnames(distr_matrix)
  f <- rep(NA, length(col_origins))

  for(i in seq_along(col_origins)) {
    matches <- which(cond_map$origin == col_origins[i])
    if(length(matches) == 1) {
      f[i] <- as.character(cond_map$condition[matches])
    } else if(length(matches) == 0) {
      stop(paste("Origin", col_origins[i], "not found in condition mapping"))
    } else {
      stop(paste("Origin", col_origins[i], "has multiple condition mappings"))
    }
  }

  # Convert to factor, preserving original levels if condition was a factor
  if(is.factor(condition)) {
    f <- factor(f, levels = levels(condition))
  } else {
    f <- factor(f)
  }

  # Create boxplot using VlnPlot2
  y_label <- if(percent) {
    if(normalize) {
      "Percentage of normalized cell counts"
    } else {
      "Percentage of cell counts"
    }
  } else {
    if(normalize) {
      "Normalized cell counts"
    } else {
      "Cell counts"
    }
  }

  # Generate plot using VlnPlot2
  p <- VlnPlot2(
    distr_matrix,
    f = f,
    violin = violin,
    box = box,
    style = style,
    pt.style = pt.style,
    pt = pt,
    pt.size = pt.size,
    pt.alpha = pt.alpha,
    width = width,
    stat.method = stat.method,
    p.adjust.method = p.adjust.method,
    label = label,
    hide.ns = hide.ns,
    cols = cols,
    lab_fill = "Condition",
    ncol = ncol,
    nrow = nrow,
    ...
  )

  # Update y-axis label
  p <- p + labs(y = y_label)

  return(p)
}

