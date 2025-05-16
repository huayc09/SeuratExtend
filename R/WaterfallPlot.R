#' @include generics.R
#'
NULL

#' @param seu A Seurat object. Only applicable when using the Seurat method.
#' @param features Features to be plotted, which can include gene expression, metrics, PC scores, or any other data that can be retrieved using the `FetchData()` function. Only applicable for the Seurat method.
#' @param group.by A variable from `meta.data` for grouping, or a character vector of the same length as the number of cells. Only applicable for the Seurat method.
#' @param cell Cell identifiers to be used in the plot. Defaults to all cells. Only applicable for the Seurat method.
#' @param cells Alternative parameter name for cell identifiers. Same functionality as 'cell'. Defaults to all cells.
#' @param slot Slot from which to retrieve feature data. Only applicable for the Seurat method.
#' @param assay Name of the assay to use. If not specified, the active assay will be used. Only applicable for the Seurat method.
#' @param priority If set to "expr", the function will fetch data from the expression matrix rather than `meta.data`. Only applicable for the Seurat method.
#' @rdname WaterfallPlot
#' @export

WaterfallPlot.Seurat <- function(
    seu,
    features,
    group.by = NULL,
    cell = NULL,
    cells = NULL,
    slot = "data",
    assay = NULL,
    priority = c("expr","none"),
    ident.1 = NULL,
    ident.2 = NULL,
    exp.transform = (length == "logFC"),
    order = TRUE,
    length = "logFC",
    color = "p",
    len.threshold = 0,
    col.threshold = 0,
    color_theme = "BuRd",
    center_color = NULL,
    top.n = NULL,
    flip = FALSE,
    y.label = NULL,
    angle = NULL,
    hjust = NULL,
    vjust = NULL,
    title = NULL,
    style = c("bar", "segment"),
    border = NA,
    log.base = "e",
    pseudocount = NULL
) {
  style <- match.arg(style)

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

  p <- WaterfallPlot.default(
    matr = t(Std.matr$matr),
    f = Std.matr$f,
    ident.1 = ident.1,
    ident.2 = ident.2,
    exp.transform = exp.transform,
    order = order,
    length = length,
    color = color,
    len.threshold = len.threshold,
    col.threshold = col.threshold,
    color_theme = color_theme,
    center_color = center_color,
    top.n = top.n,
    flip = flip,
    y.label = y.label,
    angle = angle,
    hjust = hjust,
    vjust = vjust,
    title = title,
    style = style,
    border = border,
    log.base = log.base,
    pseudocount = pseudocount
  )

  return(p)
}

#' @param matr A matrix or data frame where rows represent features and columns represent cells.
#' @param f A factor or vector indicating the identity of each cell. The length should match the number of columns in `matr`.
#' @param ident.1 The primary identity class. If not specified, the first class will be used.
#' @param ident.2 An optional secondary identity class for comparison. If NULL, comparisons will be made against all other cells.
#' @param exp.transform Indicates whether to transform data using `expm1`. This is particularly useful when calculating the log fold change of normalized gene counts. Defaults to FALSE for matrix input and TRUE for Seurat object input.
#' @param order Determines whether the features should be ordered. Defaults to TRUE.
#' @param length Specifies the statistic to determine the length of the bar. Possible values are "tscore" (default for matrix input), "p", or "logFC" (default for Seurat object input).
#' @param color Specifies the statistic to determine the color of the bar. Possible values are "tscore" (default), "p" (default for Seurat object input), or "logFC".
#' @param len.threshold Excludes features with a value for the `length` parameter below this threshold. Defaults to 0.
#' @param col.threshold Excludes features with a value for the `color` parameter below this threshold. Defaults to 0.
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
#' @param top.n Retains only the top `n` bars in both positive and negative directions. If `length(top.n)` is 1, the function retains `top.n` bars for both positive and negative directions. If `length(top.n)` is 2, it retains `top.n[1]` positive bars and `top.n[2]` negative bars. Defaults to NULL.
#' @param flip Determines whether the plot should be flipped. Defaults to TRUE for matrix input and FALSE for Seurat object input.
#' @param y.label Label for the y-axis. Defaults to "length".
#' @param angle Angle of the x-axis labels. This argument is passed to `element_text()`.
#' @param hjust Horizontal justification for the x-axis labels. This argument is passed to `element_text()`.
#' @param vjust Vertical justification for the x-axis labels. This argument is passed to `element_text()`.
#' @param title Title of the plot. Defaults to NULL.
#' @param style Character string specifying the plot style. Either "bar" (default) for traditional bar plot style, or
#'   "segment" for thin segments with end points.
#' @param border Color for the border of bars when style="bar". Use NA for no border (default), or specify a color (e.g., "black").
#' @param log.base The base for logarithmic calculations when using logFC. Can be "e" (natural logarithm, default), "2" (log2), or "10" (log10).
#' @rdname WaterfallPlot
#' @export

WaterfallPlot.default <- function(
    matr,
    f,
    ident.1 = NULL,
    ident.2 = NULL,
    exp.transform = FALSE,
    order = TRUE,
    length = "tscore",
    color = "tscore",
    len.threshold = 0,
    col.threshold = 0,
    color_theme = "BuRd",
    center_color = NULL,
    top.n = NULL,
    flip = TRUE,
    y.label = NULL,
    angle = NULL,
    hjust = NULL,
    vjust = NULL,
    title = NULL,
    style = c("bar", "segment"),
    border = NA,
    log.base = "e",
    pseudocount = NULL
) {
  style <- match.arg(style)

  scores <- WaterfallPlot_Calc(
    matr = matr,
    f = f,
    ident.1 = ident.1,
    ident.2 = ident.2,
    exp.transform = exp.transform,
    order = order,
    length = length,
    color = color,
    len.threshold = len.threshold,
    col.threshold = col.threshold,
    top.n = top.n,
    log.base = log.base,
    pseudocount = pseudocount
  )

  titles <- WaterfallPlot_Title(
    f = f,
    ident.1 = ident.1,
    ident.2 = ident.2,
    title = title,
    length_label = length,
    y.label = y.label,
    flip = flip,
    log.base = log.base)

  p <- WaterfallPlot_Plot(
    scores = scores,
    color = color,
    color_theme = color_theme,
    center_color = center_color,
    flip = flip,
    y.label = titles[[2]],
    angle = angle,
    hjust = hjust,
    vjust = vjust,
    title = titles[[1]],
    style = style,
    border = border
  )

  return(p)
}

#' @title WaterfallPlot_v3
#' @description Alias of \code{\link[SeuratExtend:WaterfallPlot]{WaterfallPlot()}}, Seurat version
#' @seealso \code{\link[SeuratExtend:WaterfallPlot]{WaterfallPlot()}}
#' @rdname WaterfallPlot_v3
#' @export

WaterfallPlot_v3 <- WaterfallPlot.Seurat

# Internal ----------------------------------------------------------------

WaterfallPlot_Calc <- function(
    matr,
    f,
    ident.1,
    ident.2,
    exp.transform,
    order,
    length,
    color,
    len.threshold,
    col.threshold,
    top.n,
    log.base = "e",
    pseudocount = NULL
) {
  library(dplyr)
  library(rlist)
  library(tidyr)

  # Validate log.base at the beginning of the function
  valid_bases <- c("e", "2", "10")
  use_natural_log <- FALSE
  
  # Convert log.base to numeric if it's a numeric string
  if (is.character(log.base) && !(log.base %in% valid_bases)) {
    # Try to convert to numeric
    numeric_base <- suppressWarnings(as.numeric(log.base))
    if (is.na(numeric_base)) {
      warning("Invalid log.base '", log.base, "' provided, using natural log (base e) instead")
      log.base <- "e"
      use_natural_log <- TRUE
    }
  }
  
  # Validate if it's a non-character, non-numeric value
  if (!is.character(log.base) && !is.numeric(log.base)) {
    warning("Invalid log.base type provided, using natural log (base e)")
    log.base <- "e"
    use_natural_log <- TRUE
  }

  f <- factor(f)
  ident.1 <- ident.1 %||% levels(f)[1]
  cell.1 <- (f == ident.1)
  cell.2 <- if(is.null(ident.2)) f != ident.1 else f == ident.2
  
  # Check if logFC calculation will be used
  will_use_logfc <- "logFC" %in% c(length, color)
  
  # Prepare matrix for analysis
  if(exp.transform) {
    original_matr <- matr
    matr <- expm1(matr)
  }
  
  # Automatically determine pseudocount if not provided and logFC will be used
  if (is.null(pseudocount) && will_use_logfc) {
    # Check the range on the transformed matrix if applicable
    check_matr <- if(exp.transform) matr else matr
    data_range <- range(check_matr, na.rm = TRUE)
    
    if (data_range[1] >= 0 && data_range[2] <= 1) {
      # 0-1 range data (like AUCell)
      pseudocount <- 0.01
      message("Data range detected as 0-1. Using pseudocount = 0.01 for logFC calculation.")
    } else {
      # Default for count data or other types
      pseudocount <- 1
      message("Using pseudocount = 1 for logFC calculation.")
    }
    
    # Check for negative values and warn if found
    if (data_range[1] < 0) {
      warning("Negative values detected in data. LogFC calculation may be affected.")
    }
  } else if (is.null(pseudocount)) {
    # Set default pseudocount silently if logFC is not used
    pseudocount <- 1
  }

  # Helper function to handle t.test with zero variance
  safe_ttest <- function(x, idx1, idx2) {
    # Handle all-zero or constant value cases
    if (all(x == 0) || length(unique(x)) == 1) {
      return(list(statistic = 0, p.value = 1))
    }
    # Perform t.test
    tryCatch({
      t.test(x[idx1], x[idx2])
    }, error = function(e) {
      list(statistic = 0, p.value = 1)
    })
  }

  scores <- list()
  if("tscore" %in% c(length, color)){
    scores[["tscore"]] <-
      apply(matr, 1, function(x) safe_ttest(x, cell.1, cell.2)[["statistic"]], simplify = TRUE)
  }
  if("p" %in% c(length, color)){
    scores[["p"]] <-
      apply(matr, 1, function(x){
        p_value <- safe_ttest(x, cell.1, cell.2)[["p.value"]]
        log10p <- pmin(-log10(p_value), 325) * ifelse(mean(x[cell.1]) > mean(x[cell.2]), 1, -1)
        return(log10p)
      }, simplify = TRUE)
  }
  if("logFC" %in% c(length, color)){
    # Calculate logFC with the specified base and pseudocount
    scores[["logFC"]] <- apply(matr, 1, function(x) {
      ratio <- mean(x[cell.1] + pseudocount) / mean(x[cell.2] + pseudocount)
      if (log.base == "e" || use_natural_log) {
        return(log(ratio))
      } else if (log.base == "2") {
        return(log2(ratio))
      } else if (log.base == "10") {
        return(log10(ratio))
      } else {
        # Must be a valid numeric at this point
        base <- as.numeric(log.base)
        return(log(ratio, base))
      }
    }, simplify = TRUE)
  }

  scores <- as.data.frame(list.cbind(scores))
  scores <- scores %>%
    mutate(rank = rownames(.),
           length = .[[length]],
           color = .[[color]])
  scores <- scores[
    abs(scores$length) > len.threshold &
      abs(scores$color) > col.threshold,
  ]

  if(order) {
    scores <- arrange(scores, desc(length))
  }
  if(!is.null(top.n)) {
    if(length(top.n) == 1) top.n <- c(top.n, top.n)
    scores2 <- split(scores, ifelse(scores$length > 0, "pos","neg"))
    scores <- rbind( head(scores2[["pos"]], top.n[1]), tail(scores2[["neg"]], top.n[2]) )
  }
  return(scores)
}

WaterfallPlot_Plot <- function(
    scores,
    color,
    color_theme,
    center_color,
    flip,
    y.label,
    angle,
    hjust,
    vjust,
    title,
    style,
    border = NA
) {
  library(ggplot2)
  library(scales)
  
  # Score preprocessing
  scores$rank <- factor(scores$rank, levels = unique(scores$rank))
  if(flip){
    scores <- scores[nrow(scores):1, ]
  }
  
  # Auto-determine angle, hjust, and vjust if not provided
  if (is.null(angle)) {
    max_label_length <- max(nchar(as.character(scores$rank)))
    if (flip) {
      angle <- 0  # Horizontal labels for flipped plot
    } else {
      angle <- if (max_label_length <= 2) 0 else -90  # Vertical for longer labels
    }
  }
  
  if (abs(angle) > 90) {
    warning("Angle should be between -90 and 90 degrees for optimal readability.")
  }
  
  if (is.null(hjust)) {
    if (angle > 0) {
      hjust <- 1  # Right align
    } else if (angle < 0) {
      hjust <- 0  # Left align
    } else {
      hjust <- 0.5  # Center align
    }
  }
  
  if (is.null(vjust)) {
    if (abs(angle) == 90) {
      vjust <- 0.5
    } else {
      vjust <- 1
    }
  }
  
  # Set fill label
  if(color == "p") lab_fill <- "-log10(p)" else lab_fill <- color

  # Base plot with common elements
  p <- ggplot(scores, aes(x = rank, y = length, colour = color)) +
    labs(fill = lab_fill, color = lab_fill, x = element_blank(), y = y.label, title = title)

  # Add style-specific elements
  if (style == "bar") {
    p <- p + geom_bar(stat = "identity", aes(fill = color), colour = border) +
      theme_classic()
  } else {
    p <- p +
      geom_segment(aes(xend = rank, y = 0, yend = length), size = 0.8) +
      geom_point(aes(y = length), size = 3) +
      theme_minimal()
  }
  p <- p +
    theme(axis.text.x = element_text(angle = angle, hjust = hjust, vjust = vjust),
          plot.title = element_text(hjust = 0.5, face = "bold"))

  value_range <- range(scores$color)
  if (style == "bar") {
    p <- p + scale_fill_cont_auto(color_theme, center_color = center_color, value_range = value_range)
  } else {
    p <- p + scale_color_cont_auto(color_theme, center_color = center_color, value_range = value_range)
  }

  if(flip) p <- p + scale_x_discrete(position = "top") + coord_flip()
  return(p)
}

WaterfallPlot_Title <- function(f, ident.1, ident.2, title, length_label, y.label, flip, log.base = "e") {
  f <- factor(f)
  ident.1 <- ident.1 %||% levels(f)[1]
  ident.2 <- ident.2 %||% paste0("non-", ident.1)
  title <- title %||% paste0(ident.1, " vs. ", ident.2)
  
  # Validate log.base for display
  valid_bases <- c("e", "2", "10")
  if (is.character(log.base) && !(log.base %in% valid_bases)) {
    # Check if it can be converted to numeric
    numeric_base <- suppressWarnings(as.numeric(log.base))
    if (is.na(numeric_base)) {
      # Invalid base, display as log(e)
      log.base <- "e"
    }
  }
  
  # Create a formatted label for logFC that includes the base
  if (length_label == "logFC") {
    if (log.base == "e") {
      # PDF-compatible format for base e (with parentheses)
      formatted_label <- "log(e)FC"
    } else if (log.base == "2") {
      # PDF-compatible format for base 2 (without parentheses)
      formatted_label <- "log2FC"
    } else if (log.base == "10") {
      # PDF-compatible format for base 10 (without parentheses)
      formatted_label <- "log10FC"
    } else {
      # Custom numeric base with PDF-compatible format (without parentheses)
      if (is.numeric(log.base) || (is.character(log.base) && !is.na(as.numeric(log.base)))) {
        formatted_label <- paste0("log", log.base, "FC")
      } else {
        # Fallback for any invalid base
        formatted_label <- "logFC"
      }
    }
  } else {
    formatted_label <- length_label
  }
  
  if(is.null(y.label)) {
    y.label <- formatted_label
    if(isTRUE(flip)) {
      # Use more common arrow symbols that work in both PNG and PDF
      y.label <- paste0(ident.2, " <- ", y.label, " -> ", ident.1)
    }
  }
  return(list(title, y.label))
}
