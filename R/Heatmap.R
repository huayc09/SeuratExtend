#' @title Heatmap
#' @description Generates a heatmap plot.
#' @param score A matrix for input, for instance, one generated using the CalcStats function.
#' @param color_scheme Specifies the color gradient for the heatmap visualization.
#'   This parameter accepts multiple input formats to provide flexibility in defining color schemes:
#'
#'     - Predefined color schemes from the `viridis` package ("A" to "H").
#'
#'     - Named vector with keys "low", "mid", and "high" for three-point gradients. Example: `c(low = muted("blue"), mid = "white", high = muted("red"))`.
#'
#'     - Two-point gradient with keys "low" and "high". Example: `c(low = "lightblue", high = "red")`.
#'
#'     - RColorBrewer sequential palettes: "Blues", "BuGn", "BuPu", "GnBu", "Greens", "Greys", "Oranges", "OrRd", "PuBu", "PuBuGn", "PuRd", "Purples", "RdPu", "Reds", "YlGn", "YlGnBu", "YlOrBr", "YlOrRd".
#'
#'     - RColorBrewer diverging palettes: "BrBG", "PiYG", "PRGn", "PuOr", "RdBu", "RdGy", "RdYlBu", "RdYlGn", "Spectral".
#'
#'     - Custom diverging palettes: "GnYlRd", "BuYlRd", "GyRd", "BuRd", "PuOr".
#'
#'     - Append "-rev" to any RColorBrewer palette name to reverse the color order. Example: "RdBu-rev".
#'
#'     - Custom color gradient using a vector of colors.
#'
#' @param center_color Logical or NULL. Determines whether the color scale should be centered at zero.
#'   If TRUE, the color scale will be centered at zero, with the midpoint color representing zero.
#'   If FALSE, the color scale will span the full range of the data without centering.
#'   If NULL (default), it will automatically determine based on the color scheme:
#'   TRUE for diverging color palettes, FALSE for sequential palettes or custom color schemes.
#'   This is particularly useful for visualizing data with both positive and negative values,
#'   such as z-scores or log fold changes.
#' @param border_color Color for the tile borders. Default: NULL.
#' @param lab_fill Label for the color. Default: 'score'.
#' @param angle Angle of the x-axis text. Passed to element_text(). Default: 45.
#' @param hjust Horizontal justification of x-axis text. Passed to element_text(). Default: 1.
#' @param vjust Vertical justification of x-axis text. Passed to element_text(). Default: 1.
#' @param legend_position Position of the figure legend. Default: 'right'.
#' @param y_text_position Position of the row name text. Default: 'right'
#' @param feature_text_subset A subset of feature names to be shown. Useful when there are many features. Default: NULL.
#' @param segment.width Adjust the width ratio of the line connecting heatmap rows and feature text. Only effective when `feature_text_subset` is set. Default: c(1, 2.5, 1).
#' @param segment.size Thickness of the line connecting heatmap rows and feature text. Only effective when `feature_text_subset` is set. Default: 0.2.
#' @param text.spacing Spacing between feature texts. Only effective when `feature_text_subset` is set. Default: 2.
#' @param text.size Font size of the feature text. Only effective when `feature_text_subset` is set. Default: 2.5.
#' @param hide_axis_line Whether to hide the axis line or not. Default: TRUE.
#' @param plot.margin Adjusts the space around the plot to prevent text labels from being cut off. This is particularly useful when axis names are too long and may extend beyond the plot boundaries. The margin is defined for top (t), right (r), bottom (b), and left (l) edges of the plot area.
#'   To increase the left margin and provide more space for the first name on the x-axis, increase the `l` value. Default: `margin(t = 5.5, r = 5.5, b = 5.5, l = 5.5)` in grid units.
#' @param expand_limits_x This parameter is softly deprecated in favor of using `plot.margin`.
#'   Previously used to increase the space on the left side of the plot by adjusting the x-axis limits. It is still available for backward compatibility but its use is not recommended. Default: NULL.
#' @param facet_col Vector or factor to split the heatmap columns. Default: NULL.
#' @param facet_row Vector or factor to split the heatmap rows. Default: NULL.
#' @param panel.spacing Spacing between panels when rows or columns are split. Default: unit(5, "pt").
#' @param strip.placement Placement of panel names when rows or columns are split. Default: 'outside'.
#' @param ncol Number of columns when either `facet_col` or `facet_row` is not NULL, wrapping the heatmap. Default: NULL.
#' @param nrow Number of rows when either `facet_col` or `facet_row` is not NULL, wrapping the heatmap. Default: NULL.
#' @param ... Additional parameters passed to the `theme` function of ggplot2.
#' @return A ggplot object.
#' @details For more detailed usage, see the examples provided.
#' @examples
#' # First, create a matrix using the CalcStats function.
#' genes <- VariableFeatures(pbmc)
#' toplot <- CalcStats(pbmc, features = genes, method = "zscore", order = "p", n = 5)
#'
#' # Generate a basic heatmap.
#' Heatmap(toplot, lab_fill = "zscore")
#'
#' # Modify the color theme to range from white to dark green.
#' Heatmap(toplot, lab_fill = "zscore", color_scheme = c("white", muted("green")))
#'
#' # Use a color theme that transitions from dark blue to light yellow (centered at 0) to dark red.
#' Heatmap(toplot, lab_fill = "zscore", color_scheme = c(
#'   low = muted("blue"),
#'   mid = "lightyellow",
#'   high = muted("red"))
#' )
#'
#' # Employ the viridis color theme, options include A, B, C, D, or E.
#' Heatmap(toplot, lab_fill = "zscore", color_scheme = "A")
#'
#' # Adjust the left margin to ensure the x-axis labels fit within the plot boundaries
#' Heatmap(toplot, lab_fill = "zscore", plot.margin = margin(l = 30))
#'
#' # Construct a dense matrix with data for 500 genes.
#' toplot2 <- CalcStats(pbmc, features = genes[1:500], method = "zscore", order = "p")
#'
#' # Display only a subset of gene names.
#' Heatmap(toplot2, lab_fill = "zscore", feature_text_subset = genes[1:20], expand_limits_x = c(-0.5, 11))
#'
#' # Divide the heatmap into rows based on gene groups.
#' gene_groups <- sample(c("group1", "group2", "group3"), nrow(toplot2), replace = TRUE)
#' Heatmap(toplot2, lab_fill = "zscore", facet_row = gene_groups) +
#'   theme(axis.text.y = element_blank())
#'
#' @rdname Heatmap
#' @export

Heatmap <- function(
    score,
    color_scheme = "BuRd",
    center_color = NULL,
    border_color = NULL,
    lab_fill = "score",
    angle = 45,
    hjust = 1,
    vjust = 1,
    legend_position = "right",
    y_text_position = "right",
    feature_text_subset = NULL,
    segment.width = c(1,2.5,1),
    segment.size = 0.2,
    text.spacing = 2,
    text.size = 2.5,
    hide_axis_line = TRUE,
    plot.margin = margin(t = 5.5, r = 5.5, b = 5.5, l = 5.5),
    expand_limits_x = NULL,
    facet_col = NULL,
    facet_row = NULL,
    panel.spacing = unit(5, "pt"),
    strip.placement = "outside",
    ncol = NULL,
    nrow = NULL,
    ...
) {
  library(ggplot2)
  library(reshape2)
  library(scales)
  library(rlang)
  library(dplyr)

  ToPlot <-
    data.frame(score, id = factor(rownames(score), levels=unique(rev(rownames(score))))) %>%
    melt()
  value_range <- range(ToPlot$value)
  ToPlot$variable <- colnames(score)[ToPlot$variable] %>% factor(levels = unique(.))
  if(!is.null(facet_col)){
    if(length(facet_col)==ncol(score)){
      if(!is.null(feature_text_subset)) facet_col <- factor(facet_col) %>% `levels<-`(c(levels(.),""))
      ToPlot$facet_col <- rep(facet_col, each = nrow(score))
    }else{
      stop('"facet_col" must be the same length as the number of input matrix columns')
    }
  }
  if(!is.null(facet_row)){
    if(!is.null(feature_text_subset)) stop('"facet_row" and "feature_text_subset" cannot be set at the same time')
    if(length(facet_row)==nrow(score)){
      ToPlot$facet_row <- rep(facet_row, times = ncol(score))
    }else{
      stop('"facet_row" must be the same length as the number of input matrix rows')
    }
  }
  border_color <- border_color %||% ifelse(ncol(score) > 80 | nrow(score) > 80, NA, "white")

  p <- ggplot(ToPlot, aes(variable, id)) +
    geom_tile(aes(fill = value), colour = border_color) +
    theme_classic()+
    labs(x = "", y = "", fill = lab_fill)+
    scale_y_discrete(position = y_text_position, expand = c(0, 0)) +
    scale_x_discrete(expand = c(0, 0))+
    theme(axis.text.x=element_text(angle = angle, hjust = hjust, vjust = vjust),
          legend.position = legend_position,
          plot.margin = plot.margin) +
    theme(...)

  p <- p + scale_fill_cont_auto(color_scheme, center_color = center_color, value_range = value_range)

  if(hide_axis_line) {
    p <- p + theme(axis.line = element_blank(),
                   axis.ticks = element_blank())
  }

  if(!is.null(expand_limits_x)) {
    p <- p + expand_limits(x = expand_limits_x)
  }

  if(!is.null(feature_text_subset)) {
    text.spacing = nrow(score) * text.spacing / 100
    text_table <-
      data.frame(text = levels(ToPlot$id) %>% .[. %in% feature_text_subset],
                 y.orig = which(levels(ToPlot$id) %in% feature_text_subset))
    text_table[nrow(text_table),"y1"] <- text_table[nrow(text_table),"y.orig"]
    if(text_table[nrow(text_table),"y1"] > nrow(score) - text.spacing / 3) {
      text_table[nrow(text_table),"y1"] <- nrow(score) - text.spacing / 3
    }
    for (i in (nrow(text_table)-1):1) {
      if(text_table[i,"y.orig"] > text_table[(i+1),"y1"] - text.spacing) {
        text_table[i,"y1"] <- text_table[(i+1),"y1"] - text.spacing
      } else {
        text_table[i,"y1"] <- text_table[i,"y.orig"]
      }
    }
    if(text_table[1, "y1"] < 1 + text.spacing / 3) {
      i <- 1
      text_table[1, "y1"] <- 1 + text.spacing / 3
      while (text_table[(i+1), "y1"] < text_table[i, "y1"] + text.spacing) {
        text_table[(i+1), "y1"] <- text_table[i, "y1"] + text.spacing
        i <- i + 1
      }
    }
    if(!is.null(facet_col)) {
      last.level <- factor("", levels = levels(facet_col))
      text_table$facet_col <- last.level
      x0 <- 0
    } else x0 <- ncol(score) + 0.5 + ncol(score)*0.01

    text_table$x1 <- x0
    text_table$x2 <- text_table$x1 + ncol(score) * segment.width[1] / 100
    text_table$x3 <- text_table$x2 + ncol(score) * segment.width[2] / 100
    text_table$x4 <- text_table$x3 + ncol(score) * segment.width[3] / 100

    p <- p +
      theme(axis.line.y = element_blank(),
            axis.text.y = element_blank(),
            axis.ticks.y = element_blank()) +
      geom_segment(data = text_table, aes(x = x1, y = y.orig, xend = x2, yend = y.orig), lineend = "round", size = segment.size) +
      geom_segment(data = text_table, aes(x = x2, y = y.orig, xend = x3, yend = y1), lineend = "round", size = segment.size) +
      geom_segment(data = text_table, aes(x = x3, y = y1, xend = x4, yend = y1), lineend = "round", size = segment.size) +
      geom_text(data = text_table, aes(label = text, x = x4 + ncol(score)*0.01, y = y1), hjust = 0, size = text.size)
  }

  if(!is.null(facet_col) | !is.null(facet_row)){
    if(!is.null(ncol) | !is.null(nrow)) {
      p <- p +
        facet_rep_wrap(
          facets = vars(facet_col), ncol = ncol, nrow = nrow,
          scales = "fixed")
    } else {
      p <- p +
        facet_grid(
          rows = if(is.null(facet_row)) NULL else vars(facet_row),
          cols = if(is.null(facet_col)) NULL else vars(facet_col),
          scales = "free", space = "free",
          switch = "y")
    }
    p <- p +
      theme(panel.grid = element_blank(),
            strip.background = element_rect(size = 0),
            panel.spacing = panel.spacing,
            strip.placement = strip.placement)
    if(!is.null(feature_text_subset)){
      g <- ggplot_gtable(ggplot_build(p))
      stripr <- which(grepl('strip-r', g$layout$name) | grepl('strip-t', g$layout$name)) %>% tail(1)
      j <- which(grepl('rect', g$grobs[[stripr]]$grobs[[1]]$childrenOrder))
      g$grobs[[stripr]]$grobs[[1]]$children[[j]]$gp$col <- NA
      p <- grid::grid.draw(g)
    }
  }
  return(p)
}


