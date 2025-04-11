#' @include generics.R
#'
NULL

#' @param seu A Seurat object. Only applicable for the Seurat method.
#' @param group.by A variable from `meta.data` for grouping or a character vector of equal length as the number of cells. Only applicable for the Seurat method.
#' @param split.by A variable from `meta.data` to bifurcate the violin plots. Only applicable for the Seurat method.
#' @param cells Cell identifiers for use. Defaults to all cells. Only applicable for the Seurat method.
#' @param slot Slot to retrieve feature data from. Only applicable for the Seurat method.
#' @param assay Name of the assay to employ. Defaults to the active assay. Only applicable for the Seurat method.
#' @param priority If set to "expr", extracts data from the expression matrix over `meta.data`. Only applicable for the Seurat method.
#' @param load.cols When TRUE, automatically loads pre-stored color information for variables from `seu@misc[["var_colors"]]`.
#' @param angle Angle for label rotation. If NULL, automatically determined based on label length. Default: NULL.
#' @param hjust Horizontal justification of labels. If NULL, automatically determined based on angle. Default: NULL.
#' @param vjust Vertical justification of labels. If NULL, automatically determined based on angle. Default: NULL.
#' @rdname VlnPlot2
#' @export

VlnPlot2.Seurat <- function(
  seu,
  features,
  group.by = NULL,
  split.by = NULL,
  cells = NULL,
  slot = "data",
  assay = NULL,
  priority = c("expr","none"),
  cols = "auto",
  load.cols = TRUE,
  ncol = NULL,
  lab_fill = "group",
  scales = "free_y",
  violin = T,
  box = T,
  width = 0.9,
  show.mean = FALSE,
  mean_colors = c("red", "blue"),
  pt = T,
  hide.outlier = F,
  pt.style = c("jitter","quasirandom"),
  pt.size = 0.2,
  pt.alpha = 1,
  strip.position = "top",
  stat.method = c("none", "wilcox.test", "t.test"),
  stats.method = NULL,
  p.adjust.method = "holm",
  label = c("p.signif","p","p.adj","p.format"),
  comparisons = NULL,
  hide.ns = TRUE,
  step.increase = 0.12,
  tip.length = 0.03,
  angle = NULL,
  hjust = NULL,
  vjust = NULL,
  ...
) {
  Std.matr <- Seu2Matr(
    seu = seu,
    features = features,
    group.by = group.by,
    split.by = split.by,
    cells = cells,
    slot = slot,
    assay = assay,
    priority = priority
  )

  cols <- VlnPlot2_SelColDisc(
    seu = seu,
    group.by = group.by,
    split.by = split.by,
    cols = cols,
    load.cols = load.cols
  )

  p <- VlnPlot2.default(
    matr = Std.matr$matr,
    f = Std.matr$f,
    f2 = Std.matr$f2,
    t = T,
    cols = cols,
    ncol = ncol,
    lab_fill = lab_fill,
    scales = scales,
    violin = violin,
    box = box,
    width = width,
    show.mean = show.mean,
    mean_colors = mean_colors,
    pt = pt,
    hide.outlier = hide.outlier,
    pt.style = pt.style,
    pt.size = pt.size,
    pt.alpha = pt.alpha,
    strip.position = strip.position,
    stat.method = stat.method,
    stats.method = stats.method,
    p.adjust.method = p.adjust.method,
    label = label,
    comparisons = comparisons,
    hide.ns = hide.ns,
    step.increase = step.increase,
    tip.length = tip.length,
    angle = angle,
    hjust = hjust,
    vjust = vjust,
    ...)
  return(p)
}

#' @param matr A matrix or data frame with rows as features and columns as cells.
#' @param f A factor or vector indicating the identity of each cell. Should match the column length of `matr`.
#' @param f2 A factor or vector akin to `f` for splitting the violin plots. Default: NULL.
#' @param features Features to depict, such as gene expression, metrics, PC scores, or any data obtainable via `FetchData()`. Default: NULL (all features in matrix).
#' @param t If the matrix has features in columns and cells in rows, transpose the matrix first. Default: FALSE.
#' @param cols Flexible color settings for the plot, accepting a variety of inputs:
#'
#'     - Seven color_pro styles: "default", "light", "pro_red", "pro_yellow", "pro_green", "pro_blue", "pro_purple".
#'
#'     - Five color_iwh styles: "iwh_default", "iwh_intense", "iwh_pastel", "iwh_all", "iwh_all_hard".
#'
#'     - Brewer color scales as specified by `brewer.pal.info`.
#'
#'     - Any manually specified colors.
#' @param ncol Specifies the number of columns for display if multiple plots are shown. Default: NULL.
#' @param lab_fill Label for the figure legend. Default: 'group'.
#' @param scales Scales parameter passed to \code{\link[ggplot2:facet_wrap]{ggplot2::facet_wrap()}}. Default: 'free_y'.
#' @param violin Indicates whether to generate a violin plot. Default: TRUE.
#' @param box Indicates whether to depict a box plot. Default: TRUE.
#' @param width Width of the box plot. Default: 0.9.
#' @param show.mean Logical value indicating whether to show mean and median lines. This is particularly useful for genes with low expression levels where the median and box plot quartiles might overlap at zero, making it difficult to interpret differences between groups. Default is FALSE.
#' @param mean_colors Vector of two colors for mean and median lines respectively. Default is c("red", "blue").
#' @param pt Indicates if points should be plotted. Default: TRUE.
#' @param hide.outlier Conceals outlier points from the box plot. Default: FALSE.
#' @param pt.style Position adjustment. Default choices: "jitter", "quasirandom".
#' @param pt.size Point size setting. Default: 0.2.
#' @param pt.alpha Adjusts the transparency of points. Default: 1.
#' @param strip.position Positions the strip ("top" (default), "bottom", "left", or "right"). Only used when `f2 = NULL`.
#' @param stat.method Determines if pairwise statistics are added to the plot. Either "wilcox.test" or "t.test". Default: "none".
#' @param stats.method Alias for \code{stat.method}. Provided for convenience but \code{stat.method} is preferred.
#'   When both are provided, \code{stats.method} takes precedence.
#' @param p.adjust.method Method for adjusting p-values, especially when conducting multiple pairwise tests or dealing with multiple grouping variables. Options include "holm", "hochberg", "hommel", "bonferroni", "BH", "BY", "fdr", and "none". Note: Adjustments are independently conducted for each variable in formulas containing multiple variables. Default: 'holm'.
#' @param label Specifies label type. Options include "p.signif" (showing significance levels), "p.format" (formatted p value), or "p", "p.adj". Default: "p.signif".
#' @param comparisons List of length-2 vectors, each containing either names of two x-axis values or two integers pointing to groups of interest for comparison. Default: all groups.
#' @param hide.ns If TRUE, the 'ns' symbol is concealed when displaying significance levels. Default: TRUE.
#' @param step.increase Numeric vector denoting the increase in fraction of total height for each additional comparison, minimizing overlap. Default: 0.12.
#' @param tip.length Numeric vector indicating the fraction of total height the bar descends to specify the exact column. For a line display instead of a bracket, set to 0. Default: 0.03.
#' @param angle Angle for label rotation. If NULL, automatically determined based on label length. Default: NULL.
#' @param hjust Horizontal justification of labels. If NULL, automatically determined based on angle. Default: NULL.
#' @param vjust Vertical justification of labels. If NULL, automatically determined based on angle. Default: NULL.
#' @param ... Further arguments passed to the \code{\link[ggpubr:stat_pvalue_manual]{ggpubr::stat_pvalue_manual()}}.
#' @rdname VlnPlot2
#' @export

VlnPlot2.default <- function(
  matr, f, f2 = NULL,
  features = NULL,
  t = F,
  cols = "pro_default",
  ncol = NULL,
  lab_fill = "group",
  scales = "free_y",
  violin = T,
  box = T,
  width = 0.9,
  show.mean = FALSE,
  mean_colors = c("red", "blue"),
  pt = T,
  hide.outlier = F,
  pt.style = c("jitter","quasirandom"),
  pt.size = 0.2,
  pt.alpha = 1,
  strip.position = "top",
  stat.method = c("none", "wilcox.test", "t.test"),
  stats.method = NULL,
  p.adjust.method = "holm",
  label = c("p.signif","p","p.adj","p.format"),
  comparisons = NULL,
  hide.ns = TRUE,
  step.increase = 0.12,
  tip.length = 0.03,
  angle = NULL,
  hjust = NULL,
  vjust = NULL,
  ...
) {

  # Handle stats.method alias
  if (!is.null(stats.method)) {
    # Issue gentle warning
    warning("Parameter 'stats.method' is an alias for 'stat.method'. Please consider using 'stat.method' in future code.", call. = FALSE)
    stat.method <- stats.method
  }

  scores <- VlnPlot2_Calc(
    matr = matr,
    f = f,
    f2 = f2,
    features = features,
    t = t
  )

  p <- VlnPlot2_Plot(
    scores = scores,
    cols = cols,
    ncol = ncol,
    lab_fill = lab_fill,
    scales = scales,
    violin = violin,
    box = box,
    width = width,
    show.mean = show.mean,
    mean_colors = mean_colors,
    pt = pt,
    hide.outlier = hide.outlier,
    pt.style = pt.style,
    pt.size = pt.size,
    pt.alpha = pt.alpha,
    strip.position = strip.position,
    angle = angle,
    hjust = hjust,
    vjust = vjust
  )

  p <- vlnplot2_Stat(
    p = p,
    stat.method = stat.method,
    p.adjust.method = p.adjust.method,
    label = label,
    comparisons = comparisons,
    hide.ns = hide.ns,
    step.increase = step.increase,
    tip.length = tip.length,
    ...
  )

  return(p)
}

#' @title StackedViolin
#' @description Alias of \code{\link[SeuratExtend:VlnPlot2]{VlnPlot2()}}
#' @seealso \code{\link[SeuratExtend:VlnPlot2]{VlnPlot2()}}
#' @rdname StackedViolin
#' @export

StackedViolin <- VlnPlot2.default


#' @rdname StackedViolin
#' @export

StackedViolin_v3 <- VlnPlot2.Seurat

# Internal ----------------------------------------------------------------

VlnPlot2_Calc <- function(
    matr,
    f,
    f2,
    features,
    t
) {
  library(reshape2)
  if(!t) matr <- t(matr)
  features <- features %||% colnames(matr)
  if(length(setdiff(features, colnames(matr))) > 0){
    message(paste0(setdiff(features, colnames(matr)), collapse = ", "), " not found")
    features <- intersect(features, colnames(matr))
  }
  f <- factor(f)
  f2 <- f2 %||% data.frame(row.names = rownames(matr))
  scores <- cbind(f, f2, as.data.frame(matr[,features,drop = F]))
  scores <- melt(scores, measure.vars = features, variable.name = "feature")
  return(scores)
}

VlnPlot2_Plot <- function(
    scores,
    cols,
    ncol,
    lab_fill,
    scales,
    violin,
    box,
    width,
    show.mean,
    mean_colors,
    pt,
    hide.outlier,
    pt.style,
    pt.size,
    pt.alpha,
    strip.position,
    angle = NULL,
    hjust = NULL,
    vjust = NULL
) {
  library(ggplot2)
  x <- ifelse(!"f2" %in% colnames(scores), "f", "f2")
  p <- ggplot(scores, aes(x = .data[[x]], y = value))
  n <- nlevels(factor(scores[[x]]))

  # Auto-determine angle, hjust, and vjust if not provided
  if (is.null(angle)) {
    max_label_length <- max(nchar(levels(scores[[x]])))
    if ("f2" %in% colnames(scores)) {
      # For split version, use -90 degrees for longer labels
      angle <- if (max_label_length <= 2) 0 else -90
    } else {
      # For non-split version, use 45 degrees for longer labels
      angle <- if (max_label_length <= 2) 0 else 45
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

  # Add violin plot if requested
  if(violin) {
    p <- p + geom_violin(mapping = aes(fill = .data[[x]]), scale = "width", width = width)
  }

  # Add box plot without violin
  if(box & !violin) {
    if(pt | hide.outlier) {
      p <- p + geom_boxplot(mapping = aes(fill = .data[[x]]), outlier.shape = NA, width = width)
    } else {
      p <- p + geom_boxplot(mapping = aes(fill = .data[[x]]), outlier.size = pt.size, width = width)
    }

    # Add mean/median lines for box without violin
    if(show.mean) {
      p <- p +
        stat_summary(
          fun = "mean",
          geom = "crossbar",
          aes(color = "Mean"),
          width = width) +
        stat_summary(
          fun = "median",
          geom = "crossbar",
          aes(color = "Median"),
          width = width) +
        scale_colour_manual(
          values = c(
            Mean = mean_colors[1],
            Median = mean_colors[2]),
          name = "")
    }
  }

  # Add points if requested
  if(pt) {
    pt.style <- pt.style[1]
    if(!pt.style %in% c("quasirandom", "jitter")) stop('"pt.style" shoule be "quasirandom" or "jitter"')
    if(pt.style == "jitter") p <- p + geom_jitter(width = width/2.2, size = pt.size, alpha= pt.alpha)
    if(pt.style == "quasirandom") {
      import("ggbeeswarm")
      p <- p + geom_quasirandom(size = pt.size, width = width/2, alpha= pt.alpha)
    }
  }

  # Add box plot with violin
  if(box & violin) {
    box_width <- if(show.mean) 0.3 else 0.12

    if(pt | hide.outlier) {
      p <- p + geom_boxplot(outlier.shape = NA, width = box_width, fill = "white")
    } else {
      p <- p + geom_boxplot(fill = "white", outlier.size = pt.size, width = box_width, outlier.alpha = pt.alpha)
    }

    # Add mean/median lines for box with violin
    if(show.mean) {
      p <- p +
        stat_summary(fun = "mean",
                     geom = "crossbar",
                     aes(color = "Mean"),
                     width = box_width) +
        stat_summary(fun = "median",
                     geom = "crossbar",
                     aes(color = "Median"),
                     width = box_width) +
        scale_colour_manual(values = c(Mean = mean_colors[1],
                                       Median = mean_colors[2]),
                            name = "")
    }
  }

  p <- p + scale_fill_disc_auto(color_scheme = cols, n = n)

  if(x == "f"){
    p <- p +
      facet_wrap(vars(feature), ncol = ncol, strip.position=strip.position, scales = scales)+
      ylab(NULL) +
      xlab(NULL) +
      theme_classic() +
      theme(strip.background = element_blank(),
            strip.placement = "outside",
            legend.position = if(show.mean) "right" else "none",
            axis.text.x = element_text(angle = angle, hjust = hjust, vjust = vjust),
            strip.text = element_text(face = "bold", size = 10)) +
      labs(fill = lab_fill) +
      scale_y_continuous(expand = expansion(mult = c(0,0.08)))
  }else{
    p <- p +
      facet_grid(vars(feature), vars(f), switch = c("both"), scales = scales)+
      ylab(NULL) +
      xlab(NULL) +
      theme_classic() +
      theme(strip.background = element_blank(),
            strip.placement = "outside",
            axis.text.x = element_blank(),
            strip.text.x = element_text(angle = angle, hjust = hjust, vjust = vjust, face = "bold", size = 10)) +
      labs(fill = lab_fill) +
      scale_y_continuous(expand = expansion(mult = c(0,0.08)))
  }
  return(p)
}

vlnplot2_Stat <- function(
    p,
    stat.method = c("none", "wilcox.test", "t.test"),
    p.adjust.method = "holm",
    label = c("p.signif","p","p.adj","p.format"),
    comparisons = NULL,
    hide.ns = TRUE,
    step.increase = 0.12,
    tip.length = 0.03,
    ...
) {
  stat.method <- stat.method[1]
  if(stat.method %in% c("wilcox.test", "t.test")) {
    library(ggpubr)
    library(dplyr)

    scores <- p$data

    if (!"f2" %in% colnames(scores)) {
      formula <- value ~ f
      group_by_arg <- "feature"
    } else {
      formula <- value ~ f2
      group_by_arg <- c("feature", "f")
    }
    stat.test <- compare_means(
      formula,
      data = scores,
      method = stat.method,
      group.by = group_by_arg,
      p.adjust.method = p.adjust.method
    )

    if(hide.ns == TRUE) {
      stat.test <- filter(stat.test, p.signif != "ns")
      if(nrow(stat.test) == 0) {
        message("No statistical significance.")
        return(p)
      }
    }
    stat.test$groups <- apply(stat.test, 1, function(x) c(x[["group1"]], x[["group2"]]), simplify = F)
    if(!is.null(comparisons)) {
      level_comparisons <- lapply(comparisons, function(pair) {
        if(is.numeric(pair)) {
          if (!"f2" %in% colnames(scores)) {
            levels(scores$f)[pair]
          } else {
            levels(scores$f2)[pair]
          }
        } else if(is.character(pair)) {
          pair
        } else {
          stop("Invalid comparison type")
        }
      })
      stat.test <- filter(stat.test, groups %in% level_comparisons)
    }
    stat.test <- vlnplot2_Stat_add_y(stat.test, scores = scores, step.increase = step.increase)
    p <- p +
      stat_pvalue_manual(
        stat.test,
        label = label[1],
        tip.length = tip.length,
        ...)
  }
  return(p)
}

vlnplot2_Stat_add_y <- function(stat.test, scores, step.increase) {
  summary_data <- scores %>%
    group_by(feature) %>%
    dplyr::summarize(
      min_value = min(value, na.rm = TRUE),
      max_value = max(value, na.rm = TRUE)
    )
  summary_data <- summary_data %>%
    dplyr::mutate(
      step = (max_value - min_value) * step.increase,
      start = (max_value - min_value) * 0.1 + max_value
    )
  grouping_cols <- if ("f2" %in% colnames(scores)) c("feature", "f") else "feature"
  stat.test <- stat.test %>%
    dplyr::left_join(summary_data, by = "feature") %>%
    dplyr::group_by(across(all_of(grouping_cols))) %>%
    dplyr::mutate(y.position = start + step * (row_number() - 1)) %>%
    dplyr::select(-min_value, -max_value, -step, -start) %>%
    ungroup()
  return(stat.test)
}

VlnPlot2_SelColDisc <- function(
    seu,
    group.by,
    split.by,
    cols,
    load.cols
) {
  if(is.null(cols)) return(NULL)
  if(cols[1] != "auto") return(cols)
  if(is.null(group.by)) group.by <- "ident"
  if(!is.null(split.by)) {
    var <- split.by
  } else {
    var <- group.by
  }
  if(length(var) == 1) {
    load_var <- seu@misc[["var_colors"]][[var]]
    if(!is.null(load_var)) {
      cols <- load_var
    } else {
      cols <- "pro_default"
    }
  }
  return(cols)
}
