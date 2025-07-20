#' @title Enhanced Dot Plot for Single-Cell Data Visualization
#' @description Creates an enhanced dot plot for visualizing gene expression across different cell types or clusters in single-cell data, with support for split visualization.
#' @param seu A Seurat object.
#' @param features A vector of gene names or a list of named vectors for grouped features.
#' @param group.by Column name in seu@meta.data for grouping cells. Default: NULL (uses current Idents).
#' @param cells Cell identifiers to be used. Defaults to all cells.
#' @param split.by Column name in seu@meta.data for splitting the groups. Default: NULL.
#' @param split.by.method Method for visualizing the split groups. Options are "border" or "color". Default: "border".
#'   - "border": Uses different border colors to represent different split groups, while the fill color represents the expression level.
#'
#'   - "color": Uses different dot colors to represent different split groups, while the transparency represents the expression level.
#' @param nudge_factor Factor to adjust the spacing between split groups. Default: 0.35.
#' @param color_scheme Color scheme for the plot. When split.by is NULL, or when split.by is specified and split.by.method is "border", this color scheme is used to represent the relative expression level (zscore). Default: 'A'.
#'    This parameter accepts multiple input formats to provide flexibility in defining color schemes:
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
#' @param split.by.colors Colors for split groups. When split.by.method is "border", this sets the border colors; when split.by.method is "color", this sets the dot colors.
#' Flexible color settings for the plot, accepting a variety of inputs:
#'
#'     - Eight color_pro styles: "default", "light", "pro_red", "pro_yellow", "pro_green", "pro_blue", "pro_purple", "bright" (default).
#'
#'     - Five color_iwh styles: "iwh_default", "iwh_intense", "iwh_pastel", "iwh_all", "iwh_all_hard".
#'
#'     - Brewer color scales as specified by `brewer.pal.info`.
#'
#'     - Any manually specified colors.
#'
#' @param center_color Center color for diverging color schemes. Default: NULL.
#' @param angle Angle of x-axis labels. Default: NULL (auto-determined).
#' @param hjust Horizontal justification of x-axis labels. Default: NULL (auto-determined).
#' @param vjust Vertical justification of x-axis labels. Default: NULL (auto-determined).
#' @param legend_position Position of the legend. Default: 'right'.
#' @param plot.margin Margins around the plot. Default: margin(t = 5.5, r = 5.5, b = 5.5, l = 5.5).
#' @param panel.spacing Spacing between facet panels. Default: unit(5, "pt").
#' @param strip.placement Placement of facet labels. Default: 'outside'.
#' @param border Whether to draw borders around points when split.by is NULL. Default: TRUE.
#' @param border.width Width of the border around points. Default: 0.6.
#' @param flip Whether to flip the coordinates of the plot. Default: FALSE.
#' @param free_space Whether to allow free space in facets. Default: TRUE.
#' @param show_grid Whether to show grid lines. Default: TRUE.
#' @param scale_percent Whether to scale the percentage values to 0-100 range (TRUE) or keep them in 0-1 range (FALSE). Default: TRUE.
#' @param legend_order Order of legends in the plot. Default: NULL (uses ggplot2 default order). When specified, must exactly match the available legend types for the current plot configuration. Use get_available_legends() helper function to see available options.
#' @param ... Additional arguments passed to theme().
#' @return A ggplot object representing the dot plot.
#' @details
#' This function creates a dot plot where the size of each dot represents the percentage of cells
#' expressing the gene, and the color represents the average expression level (corrected expression level, zscore). It supports
#' grouped features, coordinate flipping, and various customization options. The function also
#' allows for splitting the visualization by a specified variable, offering two methods of
#' representation:
#' 1. Using border colors (split.by.method = "border"): In this method, the fill color of the dots represents
#'    the expression level, while the border color represents the split group. This allows for easy
#'    comparison of expression levels across different split groups.
#' 2. Using dot colors (split.by.method = "color"): In this method, the color of the dots represents
#'    the split group, while the transparency represents the expression level. This can be useful when
#'    the focus is on comparing the distribution of split groups across different cell types or genes.
#' @examples
#' # Basic usage
#' genes <- VariableFeatures(pbmc)[1:10]
#' DotPlot2(pbmc, features = genes)
#'
#' # Grouped features
#' DotPlot2(pbmc, features = list(group1 = genes[1:3], group2 = genes[4:10]))
#'
#' # Split visualization by sample
#' DotPlot2(pbmc, features = genes, group.by = "cluster", split.by = "orig.ident", show_grid = FALSE)
#'
#' # Split visualization using colors instead of borders
#' DotPlot2(pbmc, features = genes, group.by = "cluster", split.by = "orig.ident", split.by.method = "color", show_grid = FALSE)
#'
#' # Use specific cells only
#' b_cells <- colnames(pbmc)[pbmc$cluster == "B cell"]
#' DotPlot2(pbmc, features = genes, cells = b_cells)
#'
#' # Custom settings
#' DotPlot2(pbmc, features = genes, color_scheme = "OrRd", show_grid = FALSE, border = FALSE, flip = TRUE)
#'
#' # Keep percentage values in 0-1 range
#' DotPlot2(pbmc, features = genes, scale_percent = FALSE)
#'
#' # Control legend order
#' DotPlot2(pbmc, features = genes, legend_order = c("fill", "size"))
#'
#' @rdname DotPlot2
#' @export

DotPlot2 <- function(
    seu,
    features,
    group.by = NULL,
    cells = NULL,
    split.by = NULL,
    split.by.method = "border",
    nudge_factor = 0.35,
    color_scheme = "A",
    split.by.colors = "bright",
    center_color = NULL,
    angle = NULL,
    hjust = NULL,
    vjust = NULL,
    legend_position = "right",
    plot.margin = margin(t = 5.5, r = 5.5, b = 5.5, l = 5.5),
    panel.spacing = unit(5, "pt"),
    strip.placement = "outside",
    border = TRUE,
    border.width = 0.6,
    flip = FALSE,
    free_space = TRUE,
    show_grid = TRUE,
    scale_percent = TRUE,
    legend_order = NULL,
    ...
) {
  library(ggplot2)
  library(reshape2)
  library(dplyr)

  # Validate features first
  features <- validate_features(features, seu)

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

  # Handle cells parameter for metadata subsetting
  if (!is.null(cells)) {
    # Validate cells parameter (similar to feature_percent implementation)
    if(is.logical(cells)) {
      if(length(cells) != ncol(seu)) {
        stop("Logical value of 'cells' should be the same length as cells in Seurat object")
      } else {
        cells <- colnames(seu)[cells]
      }
    } else {
      if(all(!cells %in% colnames(seu))) {
        stop("'cells' not found in Seurat object")
      }else if(any(!cells %in% colnames(seu))) {
        cells.out <- setdiff(cells, colnames(seu))
        stop(length(cells.out), " cell(s) not found in Seurat object: '",
             cells.out[1], "'...")
      }
    }
  } else {
    cells <- colnames(seu)
  }

  # Determine group levels based on the selected cells
  if (is.null(group.by)) {
    groups_subset <- Idents(seu)[cells]
    group_levels <- levels(factor(groups_subset))
  } else {
    group_values_subset <- seu@meta.data[cells, group.by]
    group_levels <- levels(factor(group_values_subset))
  }

  # Create combined group if split.by is provided
  if (!is.null(split.by)) {
    # Only consider the specified cells when determining split levels
    split_values_subset <- seu@meta.data[cells, split.by]
    split_levels <- levels(factor(split_values_subset))
    
    # Handle NULL group.by by using current Idents for combined group
    if (is.null(group.by)) {
      group_values <- as.character(Idents(seu))
    } else {
      group_values <- seu@meta.data[[group.by]]
    }
    combined_group <- paste(group_values, seu@meta.data[[split.by]], sep = "___")
    seu[["__internal_combined_group__"]] <- combined_group
    calc_group.by <- "__internal_combined_group__"
  } else {
    calc_group.by <- group.by
  }

  pct <- feature_percent(seu, tp, group.by = calc_group.by, cells = cells)
  if (scale_percent) {
    pct <- pct * 100
  }
  pct.m <- melt(pct, value.name = "pct")
  if(ncol(pct) == 1) {
    warning("Only one identity present, the mean expression values will be used")
    z <- CalcStats(seu, tp, group.by = calc_group.by, cells = cells, method = "mean") %>% as.matrix %>% melt(value.name = "zscore")
    lab_value <- "Average Expression"
  } else {
    z <- CalcStats(seu, tp, group.by = calc_group.by, cells = cells) %>% as.matrix %>% melt(value.name = "zscore")
    lab_value <- "zscore"
  }

  ToPlot <- inner_join(pct.m, z, by = c("Var1","Var2"))

  # Split combined group back into original groups
  if (!is.null(split.by)) {
    ToPlot <- ToPlot %>%
      tidyr::separate(Var2, c("group", "split"), sep = "___")
    ToPlot$group <- factor(ToPlot$group, levels = group_levels)
    ToPlot$split <-  factor(ToPlot$split, levels = split_levels)

    # Calculate nudge values
    n_splits <- length(split_levels)
    nudge_values <- seq(-nudge_factor/2, nudge_factor/2, length.out = n_splits)
    ToPlot$nudge <- nudge_values[ToPlot$split]
  } else {
    ToPlot$group <- ToPlot$Var2
    ToPlot$split <- NA
    ToPlot$nudge <- 0
  }

  if (!is.null(feature_groups)) {
    ToPlot$FeatureGroup <- rep(feature_groups, times = ncol(pct))
    ToPlot$FeatureGroup <- factor(ToPlot$FeatureGroup, levels = unique(ToPlot$FeatureGroup))
  }

  if (flip) {
    ToPlot$group <- factor(ToPlot$group, levels = rev(unique(ToPlot$group)))
    ToPlot$Var1 <- factor(ToPlot$Var1, levels = unique(ToPlot$Var1))
  } else {
    ToPlot$group <- factor(ToPlot$group, levels = unique(ToPlot$group))
    ToPlot$Var1 <- factor(ToPlot$Var1, levels = rev(unique(ToPlot$Var1)))
  }

  value_range <- range(ToPlot$zscore)

  # Determine default angle based on label lengths
  if (is.null(angle)) {
    if (flip) {
      max_label_length <- max(nchar(levels(ToPlot$Var1)))
    } else {
      max_label_length <- max(nchar(levels(ToPlot$group)))
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

  # Create the plot based on split.by, split.by.method, and border
  if (!is.null(split.by)) {
    if (split.by.method == "border") {
      p <- ggplot(ToPlot, aes(x = group, y = Var1, size = pct, fill = zscore, color = split)) +
        geom_point(shape = 21, stroke = border.width, position = position_nudge(x = ToPlot$nudge))
      color_scale <- scale_fill_cont_auto(color_scheme, center_color = center_color, value_range = value_range)
      split_scale <- scale_color_disc_auto(split.by.colors, n_splits)
      color_lab <- lab_value
    } else if (split.by.method == "color") {
      p <- ggplot(ToPlot, aes(x = group, y = Var1, size = pct, color = split, alpha = zscore)) +
        geom_point(position = position_nudge(x = ToPlot$nudge))
      color_scale <- scale_color_disc_auto(split.by.colors, n_splits)
      alpha_scale <- scale_alpha(range = c(0.1, 1))
      color_lab <- split.by
    }
  } else {
    if (border) {
      p <- ggplot(ToPlot, aes(x = group, y = Var1, size = pct, fill = zscore)) +
        geom_point(shape = 21, color = "black", stroke = border.width)
      color_scale <- scale_fill_cont_auto(color_scheme, center_color = center_color, value_range = value_range)
      color_lab <- lab_value
    } else {
      p <- ggplot(ToPlot, aes(x = group, y = Var1, size = pct, color = zscore)) +
        geom_point()
      color_scale <- scale_color_cont_auto(color_scheme, center_color = center_color, value_range = value_range)
      color_lab <- lab_value
    }
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

  if (!is.null(split.by) && split.by.method == "border") {
    p <- p + labs(color = split.by) + split_scale
  }

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

  # Apply custom legend order if specified
  if (!is.null(legend_order)) {
    available_legends <- get_available_legends(split.by, split.by.method, border)
    
    # Strict validation: must exactly match available_legends
    if (!identical(sort(legend_order), sort(available_legends))) {
      stop("legend_order must exactly match available legend types: ", 
           paste(available_legends, collapse = ", "), 
           ". You provided: ", paste(legend_order, collapse = ", "))
    }
    
    # Apply user-specified order
    guide_list <- list()
    for (i in seq_along(legend_order)) {
      legend_type <- legend_order[i]
      if (legend_type == "size") {
        guide_list[[legend_type]] <- guide_legend(order = i)
      } else if (legend_type == "fill") {
        guide_list[[legend_type]] <- guide_colorbar(order = i)
      } else if (legend_type == "color") {
        if (is.null(split.by)) {
          guide_list[[legend_type]] <- guide_colorbar(order = i)
        } else {
          guide_list[[legend_type]] <- guide_legend(order = i)
        }
      } else if (legend_type == "alpha") {
        guide_list[[legend_type]] <- guide_legend(order = i)
      }
    }
    p <- p + do.call(guides, guide_list)
  }

  return(p)
}

validate_features <- function(features, seu) {
  if (is.list(features)) {
    # Process each group
    validated_features <- lapply(features, function(f) {
      existing <- intersect(f, rownames(seu))
      if (length(existing) == 0) return(NULL)

      # Warn about missing features
      missing <- setdiff(f, rownames(seu))
      if (length(missing) > 0) {
        warning(sprintf("The following requested variables were not found: %s",
                        paste(missing, collapse = ", ")))
      }
      return(existing)
    })

    # Remove empty groups
    empty_groups <- names(validated_features)[sapply(validated_features, is.null)]
    if (length(empty_groups) > 0) {
      warning(sprintf("The following groups had no valid features and were removed: %s",
                      paste(empty_groups, collapse = ", ")))
      validated_features <- validated_features[!sapply(validated_features, is.null)]
    }

    if (length(validated_features) == 0) {
      stop("No valid features found in any group")
    }

    # Check for duplicates across all groups
    all_features <- unlist(validated_features)
    duplicates <- all_features[duplicated(all_features)]
    if (length(duplicates) > 0) {
      warning(sprintf("Removing duplicate features (keeping first occurrence): %s",
                      paste(unique(duplicates), collapse = ", ")))
      # Keep first occurrence of each feature
      seen <- c()
      validated_features <- lapply(validated_features, function(f) {
        unique_f <- setdiff(f, seen)
        seen <<- c(seen, unique_f)
        return(unique_f)
      })
      # Remove any groups that became empty after duplicate removal
      validated_features <- validated_features[sapply(validated_features, length) > 0]
    }

    return(validated_features)

  } else {
    # Process single vector of features
    existing <- intersect(features, rownames(seu))

    if (length(existing) == 0) {
      stop("No valid features found")
    }

    # Warn about missing features
    missing <- setdiff(features, rownames(seu))
    if (length(missing) > 0) {
      warning(sprintf("The following requested variables were not found: %s",
                      paste(missing, collapse = ", ")))
    }

    # Handle duplicates
    duplicates <- features[duplicated(features)]
    if (length(duplicates) > 0) {
      warning(sprintf("Removing duplicate features (keeping first occurrence): %s",
                      paste(unique(duplicates), collapse = ", ")))
      existing <- unique(existing)
    }

    return(existing)
  }
}

#' @title Get Available Legend Types for DotPlot2
#' @description Helper function to determine which legend types are available for the current DotPlot2 configuration.
#' @param split.by Split variable (same as in DotPlot2)
#' @param split.by.method Split method: "border" or "color" (same as in DotPlot2)
#' @param border Whether borders are used (same as in DotPlot2)
#' @return Character vector of available legend types
#' @examples
#' # No split, with border
#' get_available_legends(split.by = NULL, split.by.method = "border", border = TRUE)
#' # Returns: c("size", "fill")
#'
#' # Split by variable, using border method
#' get_available_legends(split.by = "condition", split.by.method = "border", border = TRUE)
#' # Returns: c("size", "fill", "color")
#'
#' # Split by variable, using color method
#' get_available_legends(split.by = "condition", split.by.method = "color", border = TRUE)
#' # Returns: c("size", "color", "alpha")
#' @export
get_available_legends <- function(split.by, split.by.method, border) {
  if (is.null(split.by)) {
    if (border) {
      return(c("size", "fill"))
    } else {
      return(c("size", "color"))
    }
  } else {
    if (split.by.method == "border") {
      return(c("size", "fill", "color"))
    } else if (split.by.method == "color") {
      return(c("size", "color", "alpha"))
    }
  }
}
