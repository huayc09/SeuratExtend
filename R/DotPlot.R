#' @title Enhanced Dot Plot for Single-Cell Data Visualization (Under Development)
#' @description Creates an enhanced dot plot for visualizing gene expression across different cell types or clusters in single-cell data. Note: This function is still under development and may change in future versions.
#' @param seu A Seurat object containing the single-cell data.
#' @param features A vector of gene names or a list of named vectors for grouped features.
#' @param group.by Column name in seu@meta.data for grouping cells. Default: NULL (uses current Idents).
#' @param color_scheme Color scheme for the plot. Default: 'A'.
#' @param center_color Center color for diverging color schemes. Default: NULL.
#' @param angle Angle of x-axis labels. Default: NULL (auto-determined).
#' @param hjust Horizontal justification of x-axis labels. Default: NULL (auto-determined).
#' @param vjust Vertical justification of x-axis labels. Default: NULL (auto-determined).
#' @param legend_position Position of the legend. Default: 'right'.
#' @param plot.margin Margins around the plot. Default: margin(t = 5.5, r = 5.5, b = 5.5, l = 5.5).
#' @param panel.spacing Spacing between facet panels. Default: unit(5, "pt").
#' @param strip.placement Placement of facet labels. Default: 'outside'.
#' @param border Whether to draw borders around points. Default: TRUE.
#' @param flip Whether to flip the coordinates of the plot. Default: FALSE.
#' @param free_space Whether to allow free space in facets. Default: TRUE.
#' @param show_grid Whether to show grid lines. Default: TRUE
#' @param ... Additional arguments passed to theme().
#' @return A ggplot object representing the dot plot.
#' @details
#' This function creates a dot plot where the size of each dot represents the percentage of cells
#' expressing the gene, and the color represents the average expression level. It supports
#' grouped features, coordinate flipping, and various customization options.
#' @examples
#' # Basic usage
#' genes <- VariableFeatures(pbmc)[1:10]
#' DotPlot2(pbmc, features = genes)
#'
#' # Grouped features
#' DotPlot2(pbmc, features = list(group1 = genes[1:3], group2 = genes[4:10]))
#'
#' # Flipped coordinates
#' DotPlot2(pbmc, features = genes, flip = TRUE)
#'
#' # Custom color scheme and grid
#' DotPlot2(pbmc, features = genes, color_scheme = "D", show_grid = FALSE)
#' @rdname DotPlot2
#' @export

DotPlot2 <- function(
    seu,
    features,
    group.by = NULL,
    color_scheme = "A",
    center_color = NULL,
    angle = NULL,
    hjust = NULL,
    vjust = NULL,
    legend_position = "right",
    plot.margin = margin(t = 5.5, r = 5.5, b = 5.5, l = 5.5),
    panel.spacing = unit(5, "pt"),
    strip.placement = "outside",
    border = TRUE,
    flip = FALSE,
    free_space = TRUE,
    show_grid = TRUE,
    ...
) {
  library(ggplot2)
  library(reshape2)
  library(dplyr)

  # Check if features is a list
  if (is.list(features)) {
    feature_groups <- unlist(lapply(names(features), function(group) {
      rep(group, length(features[[group]]))
    }))
    tp <- unlist(features)
  } else {
    tp <- features
    feature_groups <- NULL
  }

  pct <- feature_percent(seu, tp, group.by = group.by)
  pct.m <- melt(pct, value.name = "pct")
  z <- CalcStats(seu, tp, group.by = group.by) %>% as.matrix %>% melt(value.name = "z")
  ToPlot <- inner_join(pct.m, z, by = c("Var1","Var2"))

  if (!is.null(feature_groups)) {
    ToPlot$FeatureGroup <- rep(feature_groups, times = ncol(pct))
  }

  if (flip) {
    ToPlot$Var2 <- factor(ToPlot$Var2, levels = rev(unique(ToPlot$Var2)))
    ToPlot$Var1 <- factor(ToPlot$Var1, levels = unique(ToPlot$Var1))
  } else {
    ToPlot$Var2 <- factor(ToPlot$Var2, levels = unique(ToPlot$Var2))
    ToPlot$Var1 <- factor(ToPlot$Var1, levels = rev(unique(ToPlot$Var1)))
  }

  value_range <- range(ToPlot$z)

  # Determine default angle based on label lengths
  if (is.null(angle)) {
    if (flip) {
      max_label_length <- max(nchar(levels(ToPlot$Var1)))
    } else {
      max_label_length <- max(nchar(levels(ToPlot$Var2)))
    }
    angle <- if (max_label_length <= 2) 0 else 45
  }

  # Check if angle is within the recommended range
  if (abs(angle) > 90) {
    warning("Angle should be between -90 and 90 degrees for optimal readability.")
  }

  # Determine hjust based on angle
  if (is.null(hjust)) {
    if (angle > 0) {
      hjust <- 1  # Right align
    } else if (angle < 0) {
      hjust <- 0  # Left align
    } else {
      hjust <- 0.5  # Center align
    }
  }

  # Determine vjust based on angle
  if (is.null(vjust)) {
    if (abs(angle) == 90) {
      vjust <- 0.5
    } else {
      vjust <- 1
    }
  }

  # Modify the geom_point and color scaling based on the border parameter
  if (border) {
    p <- ggplot(ToPlot, aes(x = Var2, y = Var1, size = pct, fill = z)) +
      geom_point(shape = 21, color = "black", stroke = 0.5)
    color_scale <- scale_fill_cont_auto(color_scheme, center_color = center_color, value_range = value_range)
    color_lab <- "zscore"
  } else {
    p <- ggplot(ToPlot, aes(x = Var2, y = Var1, size = pct, color = z)) +
      geom_point()
    color_scale <- scale_color_cont_auto(color_scheme)
    color_lab <- "zscore"
  }

  p <- p +
    theme_bw() +
    theme(
      axis.title = element_blank(),
      axis.text.x = element_text(angle = angle, hjust = hjust, vjust = vjust),
      strip.background = element_rect(fill = NA, size = 0),
      panel.spacing = panel.spacing,
      strip.placement = strip.placement,
      legend.position = legend_position
    ) +
    labs(size = "Percent\nexpressed", color = color_lab, fill = color_lab) +
    theme(...) +
    color_scale

  if (!is.null(feature_groups)) {
    facet_scales <- ifelse(flip, "free_x", "free_y")
    facet_space <- ifelse(free_space, "free", "fixed")

    if (flip) {
      p <- p + facet_grid(cols = vars(FeatureGroup),
                          scales = facet_scales,
                          space = facet_space)
    } else {
      p <- p + facet_grid(rows = vars(FeatureGroup),
                          scales = facet_scales,
                          space = facet_space)
    }
  }

  if(!show_grid) {
    p <- p + theme(panel.grid = element_blank())
  }

  if (flip) {
    p <- p + coord_flip()
  }

  return(p)
}
