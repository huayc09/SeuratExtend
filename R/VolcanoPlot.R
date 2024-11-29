#' @include generics.R
#'
NULL

#' @param seu A Seurat object. Only applicable when using the Seurat method.
#' @param features Features to be plotted, which can include gene expression or any other data that can be retrieved using the `FetchData()` function. Only applicable for the Seurat method.
#' @param group.by A variable from `meta.data` for grouping, or a character vector of the same length as the number of cells. Only applicable for the Seurat method.
#' @param cell Cell identifiers to be used in the plot. Defaults to all cells. Only applicable for the Seurat method.
#' @param cells Alternative parameter name for cell identifiers. Same functionality as 'cell'. Defaults to all cells.
#' @param slot Slot from which to retrieve feature data. Only applicable for the Seurat method.
#' @param assay Name of the assay to use. If not specified, the active assay will be used. Only applicable for the Seurat method.
#' @param priority If set to "expr", the function will fetch data from the expression matrix rather than `meta.data`. Only applicable for the Seurat method.
#' @rdname VolcanoPlot
#' @export

VolcanoPlot.Seurat <- function(
    seu,
    features = rownames(seu),
    group.by = NULL,
    cell = NULL,
    cells = NULL,
    slot = "data",
    assay = NULL,
    priority = c("expr","none"),
    ident.1 = NULL,
    ident.2 = NULL,
    exp.transform = TRUE,
    x = "logFC",
    y = "p",
    x.threshold = NULL,
    y.threshold = NULL,
    x.quantile = 0.99,
    y.quantile = 0.99,
    top.n = 10,
    color = c("grey20", "grey50", "red3"),
    title = NULL
) {

  # Use cells if provided, otherwise fall back to cell
  cell_subset <- cells %||% cell

  Std.matr <- Seu2Matr(
    seu = seu,
    features = features,
    group.by = group.by,
    cells = cell_subset,
    slot = slot,
    assay = assay,
    priority = priority
  )

  p <- VolcanoPlot.default(
    matr = t(Std.matr$matr),
    f = Std.matr$f,
    ident.1 = ident.1,
    ident.2 = ident.2,
    exp.transform = exp.transform,
    x = x,
    y = y,
    x.threshold = x.threshold,
    y.threshold = y.threshold,
    x.quantile = x.quantile,
    y.quantile = y.quantile,
    top.n = top.n,
    color = color,
    title = title
  )

  return(p)
}

#' @param matr A matrix or data frame where rows represent features and columns represent cells.
#' @param f A factor or vector indicating the identity of each cell. The length should match the number of columns in `matr`.
#' @param ident.1 The primary identity class. If not specified, the first class will be used.
#' @param ident.2 An optional secondary identity class for comparison. If NULL, comparisons will be made against all other cells.
#' @param exp.transform Indicates whether to transform data using `expm1`. This is particularly useful when calculating the log fold change of normalized gene counts. Defaults to TRUE for Seurat object input and FALSE for matrix input.
#' @param x A character specifying the statistic for x-axis. Must be one of c("logFC", "tscore"). Defaults to "logFC".
#' @param y A character specifying the statistic for y-axis. Must be one of c("p", "tscore"). Defaults to "p".
#' @param x.threshold Threshold for the x-axis statistic. Features beyond this threshold will be labeled. Defaults to NULL (auto-detected).
#' @param y.threshold Threshold for the y-axis statistic. Features beyond this threshold will be labeled. Defaults to NULL (auto-detected).
#' @param x.quantile Quantile threshold for x-axis (between 0 and 1) used when auto-calculating x.threshold. Defaults to 0.99.
#' @param y.quantile Quantile threshold for y-axis (between 0 and 1) used when auto-calculating y.threshold. Defaults to 0.99.
#' @param top.n Number of top features to label. Defaults to 10.
#' @param color A vector of three colors specifying the color scheme: c(center_color, non_significant_color, significant_color).
#'   Defaults to c("grey20", "grey50", "red3"). The first color is used for points near zero (below thresholds),
#'   the second for points passing single threshold, and the third for points passing both thresholds.
#' @param title Title of the plot. Defaults to NULL.
#' @rdname VolcanoPlot
#' @export

VolcanoPlot.default <- function(
    matr,
    f,
    ident.1 = NULL,
    ident.2 = NULL,
    exp.transform = FALSE,
    x = "logFC",
    y = c("p", "tscore"),
    x.threshold = NULL,
    y.threshold = NULL,
    x.quantile = 0.99,
    y.quantile = 0.99,
    top.n = 10,
    color = c("grey20", "grey50", "red3"),
    title = NULL
){
  # Match y argument
  y <- match.arg(y)

  scores <- WaterfallPlot_Calc(
    matr = matr,
    f = f,
    ident.1 = ident.1,
    ident.2 = ident.2,
    exp.transform = exp.transform,
    order = FALSE,
    length = x,
    color = y,
    len.threshold = 0,
    col.threshold = 0,
    top.n = NULL
  )

  scores$color <- abs(scores$color)
  scores[[y]] <- abs(scores[[y]])

  titles <- WaterfallPlot_Title(
    f = f,
    ident.1 = ident.1,
    ident.2 = ident.2,
    title = title,
    length_label = x,
    y.label = NULL,
    flip = TRUE)

  p <- VolcanoPlot_Plot(
    scores = scores,
    x = x,
    y = y,
    x.threshold = x.threshold,
    y.threshold = y.threshold,
    x.quantile = x.quantile,
    y.quantile = y.quantile,
    top.n = top.n,
    color = color,
    title = titles[[1]],
    x.label = titles[[2]]
  )

  return(p)
}

# Internal ----------------------------------------------------------------

VolcanoPlot_Plot <- function(
    scores,
    x,
    y,
    x.threshold,
    y.threshold,
    x.quantile,
    y.quantile,
    top.n,
    color,
    title,
    x.label
) {
  library(ggplot2)
  library(ggrepel)
  library(dplyr)

  # Auto-calculate thresholds if not provided
  if(is.null(x.threshold)) {
    x.threshold <- get_auto_threshold(scores$length, quantile = x.quantile)
  }
  if(is.null(y.threshold)) {
    y.threshold <- get_auto_threshold(scores$color, quantile = y.quantile)
  }

  # Determine which points to label
  scores$significant <- abs(scores[[x]]) > x.threshold & scores[[y]] > y.threshold
  scores$label <- ifelse(scores$significant, rownames(scores), NA)

  # Determine point colors based on thresholds
  scores$point_color <- case_when(
    abs(scores$length) <= x.threshold & scores$color <= y.threshold ~ color[1],  # Near zero points
    abs(scores$length) > x.threshold & scores$color > y.threshold ~ color[3],    # Significant points
    TRUE ~ color[2]                                                              # Other points
  )

  # Select top features for each direction
  top_features_up <- scores %>%
    filter(significant & length > 0) %>%
    arrange(desc(color)) %>%
    head(top.n) %>%
    pull(label)

  top_features_down <- scores %>%
    filter(significant & length < 0) %>%
    arrange(desc(color)) %>%
    head(top.n) %>%
    pull(label)

  scores$label <- ifelse(scores$label %in% c(top_features_up, top_features_down), scores$label, NA)

  # Create the plot
  p <- ggplot(scores, aes(x = length, y = color, color = point_color, label = label)) +
    geom_point() +
    geom_text_repel(na.rm = TRUE) +
    theme_bw() +
    geom_hline(yintercept = y.threshold, linetype = "dashed") +
    geom_vline(xintercept = c(-x.threshold, x.threshold), linetype = "dashed") +
    labs(x = x.label, y = if(y == "p") "-log10(p)" else y, title = title) +
    scale_color_identity() +
    scale_y_continuous(expand = expansion(mult = c(0, 0.05))) +
    theme(legend.position = "none") +
    CenterTitle() + theme(plot.title = element_text(face = "bold"))

  return(p)
}

# Internal function to get appropriate scale for the data
get_data_scale <- function(x) {
  max_abs <- max(abs(x), na.rm = TRUE)
  10^floor(log10(max_abs))
}

# Get nice numbers at appropriate scale
get_nice_numbers <- function(scale) {
  base_numbers <- c(0.1, 0.2, 0.25, 0.5, 1, 1.5, 2, 2.5, 3, 4, 5, 6, 8, 10)
  sort(unique(c(base_numbers * scale, base_numbers * (scale/10))))
}

get_auto_threshold <- function(x, quantile = 0.99) {
  # Handle edge cases
  if (length(x) == 0 || all(is.na(x))) return(1)
  if (all(x == 0, na.rm = TRUE)) return(0.1)

  # Get data scale and nice numbers
  scale <- get_data_scale(x)
  nice_numbers <- get_nice_numbers(scale)

  # Calculate threshold
  x_abs <- abs(x)
  quant_val <- quantile(x_abs, quantile, na.rm = TRUE)

  # Find closest nice number
  nice_threshold <- nice_numbers[which.min(abs(nice_numbers - quant_val))]

  return(nice_threshold)
}
