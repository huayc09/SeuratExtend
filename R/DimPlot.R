#' @title Create an Enhanced Dimensional Reduction Plot
#' @description This function creates a dimension reduction plot that can handle both discrete and continuous variables seamlessly. It incorporates additional customization options for visual representation and automatically recognizes input variable types to optimize visualization.
#' @param seu Seurat object containing single-cell data for visualization.
#' @param features Variables to be visualized in the plot, accepting both discrete and continuous variables. Default: NULL.
#' @param group.by Alias for `features`. Default: NULL.
#' @param split.by A metadata column name by which to split the plot, creating separate plots for each unique value. This can be useful for visualizing differences across conditions or experiments.
#'   Default: NULL.
#' @param cells A vector specifying a subset of cells to include in the plot.
#'   Default: all cells are included.
#' @param slot Which data slot to use for pulling expression data. Accepts 'data', 'scale.data', or 'counts'. Default: 'data'.
#' @param assay Specify the assay from which to retrieve data.
#'   Default: NULL, which will use the default assay.
#' @param dims A two-length numeric vector specifying which dimensions to use for the x and y axes, typically from a PCA, tSNE, or UMAP reduction.
#'   Default: c(1, 2).
#' @param reduction Which dimensionality reduction to use. If not specified, will search in order of 'umap', 'tsne', then 'pca'.
#'   Default: NULL.
#' @param priority Specifies which to prioritize when metadata column names conflict with gene names: 'expr' for expression, 'none' for metadata.
#'   Default: c("expr", "none").
#' @param ncol Number of columns to display when combining multiple plots into a single patchworked ggplot object.
#'   Default: NULL.
#' @param nrow Number of rows to display when combining multiple plots into a single patchworked ggplot object.
#'   Default: NULL.
#' @param nrow.each Specifies the number of rows each split plot should have when using the split.by parameter.
#'   Default: NULL.
#' @param ncol.legend Integer specifying the number of columns in the plot legend of discrete variables. Default: NULL.
#' @param cols Flexible color settings for the plot, accepting a variety of inputs:
#'
#'   - A vector specifying a global color setting similar to Seurat's `DimPlot`/`FeaturePlot`.
#'
#'   - A list specifying colors for each variable type (discrete/continuous) or for each individual variable. For example, `list(discrete = "auto", continuous = "A")` applies automatic styling from `color_pro()` for discrete variables and `viridis` "A" style for continuous variables. More detailed setups can be included. For example, `list("cluster" = "light", "CD14" = "OrRd")`.
#'
#'   For continuous variables:
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
#'     - Append "-rev" to any RColorBrewer palette name to reverse the color order. Example: "Spectral-rev".
#'
#'     - Custom color gradient using a vector of colors.
#'
#'   For discrete variables:
#'
#'     - Eight color_pro styles: "default", "light", "pro_red", "pro_yellow", "pro_green", "pro_blue", "pro_purple", "bright".
#'
#'     - Five color_iwh styles: "iwh_default", "iwh_intense", "iwh_pastel", "iwh_all", "iwh_all_hard".
#'
#'     - Brewer color scales as specified by `brewer.pal.info`.
#'
#'     - Any manually specified colors.
#'
#'   Default: list().
#' @param load.cols When TRUE, automatically loads pre-stored color information for variables from `seu@misc[["var_colors"]]`.
#'   Default: TRUE.#' @param pt.size Point size for plotting, adjusts the size of each cell in the plot.
#'   Default: NULL.
#' @param shape.by Metadata column or expression data used to specify different shapes for cells in the plot, allowing for additional visual distinctions.
#'   Default: NULL.
#' @param alpha.by Transparency of points in the plot, which can be helpful in densely plotted areas.
#'   Default: NULL.
#' @param order Boolean determining whether to plot cells in order of expression. Can be useful if points representing higher expression levels are being buried.
#'   Default: c(discrete = FALSE, continuous = TRUE).
#' @param shuffle Randomly shuffles the order of points to prevent burial under other points.
#'   Default: c(discrete = TRUE, continuous = FALSE).
#' @param label Whether to label clusters or other features in the plot.
#'   Default: FALSE.
#' @param label.color Customize the color of labels; defaults to the same as cluster colors unless specified, such as "black".
#'   Default: NULL.
#' @param box Whether to draw a box around labels to enhance visibility.
#'   Default: FALSE.
#' @param index.title Specify a prefix for cluster indices when labels are replaced by numerical indices to simplify the plot.
#'   Default: NULL.
#' @param repel Whether to use a repelling algorithm to avoid overlapping text labels.
#'   Default: FALSE.
#' @param label.size Size of the text labels used for clusters or features.
#'   Default: 4.
#' @param theme Customizes ggplot themes and other plot elements. Since DimPlot2 returns a combined plot object (using plot_grid), standard ggplot syntax like \code{+ theme()} cannot be applied directly. This parameter provides flexible ways to customize plot appearance:
#'
#'   \strong{Single theme elements:}
#'   \itemize{
#'     \item \code{theme = NoAxes()} - Remove axes (Seurat function)
#'     \item \code{theme = labs(title = "My Title", color = "Expression")} - Add labels
#'     \item \code{theme = theme_minimal()} - Apply ggplot2 themes
#'   }
#'
#'   \strong{Simple combinations (limited cases):}
#'   \itemize{
#'     \item \code{theme = NoAxes() + NoLegend()} - Works for some Seurat theme functions
#'   }
#'
#'   \strong{Complex combinations (recommended for multiple elements):}
#'
#'   For multiple theme elements that cannot be directly combined with \code{+}, provide a list:
#'   \itemize{
#'     \item \code{theme = list(NoAxes(), labs(title = "My Title", color = "Expression"))}
#'     \item \code{theme = list(theme_minimal(), labs(color = "Gene Expression"), theme(legend.position = "bottom"))}
#'   }
#'
#'   Each element in the list will be applied sequentially to individual ggplot objects before combining.
#'   Default: NULL.
#' @param cells.highlight A vector of cell names to highlight; simpler input than Seurat's approach, focusing on ease of use.
#'   Default: NULL.
#' @param cols.highlight A color or vector of colors to use for highlighting specified cells; will repeat to match the number of groups in cells.highlight.
#'   Default: '#DE2D26'.
#' @param sizes.highlight Size of highlighted points, providing emphasis where needed.
#'   Default: 1.
#' @param na.value Color to use for NA points when using a custom color scale.
#'   Default: 'grey80'.
#' @param raster Whether to convert the plot points to a raster format, which can help with performance on large datasets.
#'   Default: NULL.
#' @param raster.dpi The resolution for rasterized plots, useful for maintaining detail in dense plots.
#'   Default: NULL.
#' @param combine Whether to combine multiple plots into a single ggplot object using patchwork.
#'   Default: TRUE.
#' @param align Specifies how plots should be aligned if combined, accepting 'h' for horizontal, 'v' for vertical, or 'hv' for both.
#'   Default: 'hv'.
#' @return A ggplot object if `combine` is TRUE; otherwise, a list of ggplot objects, allowing for flexible plot arrangements or combined visualizations.
#' @details `DimPlot2` extends the functionality of Seurat's visualization tools by combining the features of `DimPlot` and `FeaturePlot` into a single, more versatile function. It automatically recognizes whether the input features are discrete or continuous, adjusting the visualization accordingly. This makes `DimPlot2` ideal for exploring complex scRNA-seq data without the need to switch between different plotting functions based on variable types. The function also offers advanced customization options for colors, themes, and labeling, making it highly adaptable to various data visualization needs.
#' @examples
#' library(Seurat)
#' library(SeuratExtend)
#'
#' # Create a basic dimensional reduction plot with default settings
#' DimPlot2(pbmc)
#'
#' # Visualize different variables, including both discrete and continuous types
#' DimPlot2(pbmc, features = c("cluster", "orig.ident", "CD14", "CD3D"))
#'
#' # Split the visualization by a specific variable for comparative analysis
#' DimPlot2(pbmc, features = c("cluster", "CD14"), split.by = "orig.ident", ncol = 1)
#'
#' # Highlight specific cells, such as a particular cluster
#' b_cells <- colnames(pbmc)[pbmc$cluster == "B cell"]
#' DimPlot2(pbmc, cells.highlight = b_cells)
#'
#' # Apply advanced customization for colors and themes
#' DimPlot2(
#'   pbmc,
#'   features = c("cluster", "orig.ident", "CD14", "CD3D"),
#'   cols = list(
#'     "cluster" = "pro_blue",
#'     "CD14" = "D",
#'     "CD3D" = c("#EEEEEE", "black")
#'   ),
#'   theme = NoAxes())
#'
#' # Enhance the plot with labels and bounding boxes
#' DimPlot2(pbmc, label = TRUE, box = TRUE, label.color = "black", repel = TRUE, theme = NoLegend())
#'
#' # Use multiple theme elements together
#' DimPlot2(
#'   pbmc,
#'   features = c("CD14", "CD3D"),
#'   theme = list(NoAxes(), labs(title = "Gene Expression", color = "Expression Level")))
#'
#' # Complex theme customization
#' DimPlot2(
#'   pbmc,
#'   features = c("cluster", "CD14"),
#'   theme = list(
#'     theme_minimal(),
#'     labs(color = "Expression"),
#'     theme(legend.position = "bottom", plot.title = element_text(size = 16))
#'   ))
#'
#' # Use indices instead of long cluster names to simplify labels in the plot
#' DimPlot2(pbmc, index.title = "C", box = TRUE, label.color = "black")
#' @rdname DimPlot2
#' @export

DimPlot2 <- function(
    seu,
    features = NULL,
    group.by = NULL,
    split.by = NULL,
    cells = NULL,
    slot = "data",
    assay = NULL,
    dims = c(1, 2),
    reduction = NULL,
    priority = c("expr", "none"),
    ncol = NULL,
    nrow = NULL,
    nrow.each = NULL,
    ncol.legend = NULL,
    cols = list(),
    load.cols = TRUE,
    pt.size = NULL,
    shape.by = NULL,
    alpha.by = NULL,
    order = c(discrete = FALSE, continuous = TRUE),
    shuffle = c(discrete = TRUE, continuous = FALSE),
    label = FALSE,
    label.color = NULL,
    box = FALSE,
    index.title = NULL,
    repel = FALSE,
    label.size = 4,
    theme = NULL,
    cells.highlight = NULL,
    cols.highlight = "#DE2D26",
    sizes.highlight = 1,
    na.value = "grey80",
    raster = NULL,
    raster.dpi = NULL,
    combine = TRUE,
    align = "hv"
) {
  plot_data <- DimPlot2_GetData(
    seu = seu,
    features = features,
    group.by = group.by,
    split.by = split.by,
    cells = cells,
    slot = slot,
    assay = assay,
    dims = dims,
    reduction = reduction,
    shape.by = shape.by,
    alpha.by = alpha.by,
    cells.highlight = cells.highlight
  )

  data.dim <- plot_data$data.dim
  data.var <- plot_data$data.var
  vars <- colnames(data.var)
  n_features <- length(vars)
  raster <- raster %||% (nrow(x = data.dim) > 1e+05)
  if ((nrow(x = data.dim) > 1e+05) & !isFALSE(raster)) {
    message("Rasterizing points since number of points exceeds 100,000.",
            "\nTo disable this behavior set `raster=FALSE`")
  }

  if (!is.null(x = raster.dpi)) {
    if (!is.numeric(x = raster.dpi) || length(x = raster.dpi) != 2)
      stop("'raster.dpi' must be a two-length numeric vector")
  } else raster.dpi <- c(512, 512)

  pt.size <- pt.size %||% DimPlot2_AutoPointSize(
    data = data.dim,
    raster = raster,
    n_features = n_features
  )

  get_disc_cont_par <- function(par, type, default) {
    if (is.vector(par) && is.logical(par) && type %in% names(par)) {
      par.d <- par[type]
    } else if (is.logical(par) && length(par) == 1) {
      par.d <- par
    } else {
      par.d <- default
    }
    return(par.d)
  }

  library(ggplot2)
  p <- list()
  for (i in vars) {
    data.single <- data.dim
    data.single$var <- data.var[[i]]
    if(is_continuous(data.single$var)) {
      scale_color <- DimPlot2_SelColCont(
        seu = seu,
        var = i,
        cols = cols,
        load.cols = load.cols
      )
      order.c <- get_disc_cont_par(order, "continuous", TRUE)
      shuffle.c <- get_disc_cont_par(shuffle, "continuous", FALSE)

      p[[i]] <- DimPlot2_PlotSingle(
        data.plot = data.single,
        title = i,
        cols = scale_color,
        pt.size = pt.size,
        order = order.c,
        shuffle = shuffle.c,
        label = FALSE,
        label.color = NULL,
        repel = FALSE,
        box = FALSE,
        label.size = 4,
        cols.highlight = "#DE2D26",
        sizes.highlight = 1,
        na.value = na.value,
        nrow.each = nrow.each,
        theme = theme,
        raster = raster,
        raster.dpi = raster.dpi)
    } else {
      if(i != "Selected_cells") data.single$var <- factor(data.single$var)
      n <- nlevels(data.single$var)
      if(!is.null(index.title)) {
        data.single$var_orig <- factor(data.single$var)
        index <- paste0(index.title, seq(nlevels(data.single$var_orig)))
        index_cluster <- paste(index, levels(data.single$var_orig))
        data.single$var <- factor(index[data.single$var_orig], levels = index)
        data.single$index_cluster <- factor(index_cluster[data.single$var_orig], levels = index_cluster)
        label <- TRUE
        labels <- index_cluster
      } else {labels <- waiver()}
      if(n > 100) stop("Discrete variable '",i,"' has > 100 values. Unable to plot.")
      if(i != "Selected_cells") {
        scale_color <- DimPlot2_SelColDisc(
          seu = seu,
          n = n,
          var = i,
          cols = cols,
          load.cols = load.cols,
          label = label,
          labels = labels,
          box = box
        )
      } else scale_color <- NULL

      order.d <- get_disc_cont_par(order, "discrete", FALSE)
      shuffle.d <- get_disc_cont_par(shuffle, "discrete", TRUE)

      p[[i]] <- DimPlot2_PlotSingle(
        data.plot = data.single,
        title = i,
        cols = scale_color,
        pt.size = pt.size,
        order = order.d,
        shuffle = shuffle.d,
        label = label,
        label.color = label.color,
        repel = repel,
        box = box,
        label.size = label.size,
        cols.highlight = cols.highlight,
        sizes.highlight = sizes.highlight,
        na.value = na.value,
        nrow.each = nrow.each,
        theme = theme,
        raster = raster,
        raster.dpi = raster.dpi,
        ncol.legend = ncol.legend)
    }
  }
  if(!combine) {
    return(p)
  } else {
    import("cowplot")
    p <- plot_grid(plotlist = p, ncol = ncol, nrow = nrow, align = align)
    return(p)
  }
}

is_continuous <- function(vec) {
  return(is.numeric(vec) & !is.factor(vec))
}

DimPlot2_GetData <- function(
    seu,
    features = NULL,
    group.by = NULL,
    split.by = NULL,
    cells = NULL,
    slot = "data",
    assay = NULL,
    dims = c(1, 2),
    reduction = NULL,
    shape.by = NULL,
    alpha.by = NULL,
    cells.highlight = NULL
) {
  if(!require(SeuratObject)) library(Seurat)

  if(!is.null(assay)) DefaultAssay(seu) <- assay

  # cells
  cells <- cells %||% colnames(seu)
  if(is.logical(cells)) {
    if(length(cells) != ncol(seu)) {
      stop("Logical value of 'cells' should be the same length as cells in ",
           "Seurat object")
    } else {
      cells.l <- cells
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
    cells.l <- colnames(seu) %in% cells
  }

  # reduction
  if (length(x = dims) != 2 || !is.numeric(x = dims)) {
    stop("'dims' must be a two-length integer vector")
  }
  reduction <- reduction %||% DefaultDimReduc(object = seu)
  data.dim <- Embeddings(object = seu[[reduction]])[cells, dims]
  data.dim <- as.data.frame(x = data.dim)
  dims <- colnames(data.dim)

  # vars
  if(!is.null(cells.highlight)) {
    seu@meta.data[["Selected_cells"]] <- factor(
      ifelse(
        colnames(seu) %in% cells.highlight,
        "Selected", "Unselected"),
      levels = c("Unselected","Selected")
    )
  }
  if(is.null(features) & is.null(group.by) & is.null(cells.highlight)) {
    data.var <- data.frame(Idents = factor(Idents(seu)[cells]))
  } else {
    vars <- c(features, group.by)
    if(!is.null(cells.highlight)) vars = c(vars, "Selected_cells")
    if (utils::packageVersion("SeuratObject") >= "5.0.0") {
      data.var <- FetchData(object = seu, vars = vars, cells = cells, layer = slot, clean = "none")
    } else {
      data.var <- FetchData(object = seu, vars = vars, cells = cells, slot = slot)
    }
  }

  # split.by, shape.by, alpha.by
  plot_vars <- list(
    split.by = split.by,
    shape.by = shape.by,
    alpha.by = alpha.by)
  for (i in names(plot_vars)) {
    if(!is.null(plot_vars[[i]])) {
      data.dim[[i]] <- DimPlot2_PlotVars(
        seu = seu,
        var = plot_vars[[i]],
        var_name = i,
        cells = cells,
        cells.l = cells.l)
    }
  }

  return(list(
    data.dim = data.dim,
    data.var = data.var
  ))
}

DimPlot2_PlotVars <- function(
    seu,
    var,
    var_name,
    cells,
    cells.l
) {
  if(is.null(var)) {
    f2 <- NULL
  }else if(length(var) == 1) {
    if(var %in% colnames(seu@meta.data)){
      f2 <- factor(seu[[var]][cells,])
      names(f2) <- cells
    }else{
      stop("Cannot find '", var, "' in meta.data")
    }
  }else if(length(var) == ncol(seu)) {
    f2 <- factor(var[cells.l])
    names(f2) <- cells
  }else if(length(var) == length(cells)) {
    f2 <- factor(var)
    names(f2) <- cells
  }else{
    stop("'",var_name,"' should be variable name in 'meta.data' or ",
         "string with the same length of cells")
  }
  return(f2)
}

DimPlot2_PlotSingle <- function (
    data.plot,
    title = NULL,
    cols = NULL,
    pt.size = NULL,
    order = NULL,
    shuffle = NULL,
    label = FALSE,
    label.color = NULL,
    repel = FALSE,
    box = FALSE,
    label.size = 4,
    cols.highlight = "#DE2D26",
    sizes.highlight = 1,
    na.value = "grey80",
    nrow.each = NULL,
    theme = NULL,
    raster = NULL,
    raster.dpi = NULL,
    ncol.legend = NULL)
{
  library(ggplot2)
  dims <- colnames(data.plot)[1:2]

  if (isTRUE(x = shuffle)) {
    data.plot <- data.plot[sample(x = 1:nrow(x = data.plot)), ]
  }
  istrue.cell.highlight <- (title == "Selected_cells" & identical(levels(data.plot$var), c("Unselected","Selected")))
  if (isTRUE(x = order | istrue.cell.highlight)) {
    data.plot <- data.plot[order(data.plot$var), ]
  }

  if(!"shape.by" %in% colnames(data.plot)) shape.by <- NULL
  if(!"alpha.by" %in% colnames(data.plot)) alpha.by <- NULL

  plot <- ggplot(data = data.plot, mapping = aes(x = .data[[dims[1]]], y = .data[[dims[2]]], color = var, shape = shape.by, alpha = alpha.by))
  if(istrue.cell.highlight) pt.size <- sizes.highlight
  plot <-
    if (isTRUE(x = raster)) {
      import("scattermore")
      plot + geom_scattermore(pointsize = pt.size, pixels = raster.dpi)
    } else {
      plot + geom_point(size = pt.size)
    }
  if(!is_continuous(data.plot$var)) {
    plot <- plot +
      guides(color = guide_legend(override.aes = list(size = 3), ncol = ncol.legend))
  }
  plot <- plot + labs(color = NULL, title = title)
  if (label & title != "Selected_cells") {
    plot <- DimPlot2_LabelClusters(
      plot = plot,
      id = "var",
      repel = repel,
      label.color = label.color,
      box = box,
      size = label.size)
  }
  if("split.by" %in% colnames(data.plot)) {
    plot <- plot + facet_wrap(vars(split.by), nrow = nrow.each)
  }
  import("cowplot")
  plot <- plot + theme_cowplot() + CenterTitle()
  
  # Apply theme elements
  plot <- apply_theme_elements(plot, theme)
  
  if(istrue.cell.highlight) {
    plot <- plot + scale_color_manual(values = c("#C3C3C3", cols.highlight))
  } else {
    plot <- plot + cols
  }
  return(plot)
}

DimPlot2_AutoPointSize <- function(data, raster = NULL, n_features = 1) {
  # When raster is TRUE, return 1
  if (isTRUE(x = raster)) {
    return(1)
  }

  # Calculate base point size
  base_size <- min(1583/nrow(x = data), 1)

  # Adjust size based on number of features
  # Using square root of n_features to scale down point size
  adjusted_size <- base_size / sqrt(n_features)

  return(adjusted_size)
}

DimPlot2_LabelClusters <- function(
    plot,
    id,
    clusters = NULL,
    labels = NULL,
    label.color = NULL,
    split.by = NULL,
    repel = TRUE,
    box = FALSE,
    geom = 'GeomPoint',
    position = "median",
    ...
) {
  if(repel) import("ggrepel")
  xynames <- unlist(x = GetXYAesthetics(plot = plot, geom = geom), use.names = TRUE)
  if (!id %in% colnames(x = plot$data)) {
    stop("Cannot find variable ", id, " in plotting data")
  }
  if (!is.null(x = split.by) && !split.by %in% colnames(x = plot$data)) {
    warning("Cannot find splitting variable ", id, " in plotting data")
    split.by <- NULL
  }
  data <- plot$data[, c(xynames, id, split.by)]
  possible.clusters <- as.character(x = na.omit(object = unique(x = data[, id])))
  groups <- clusters %||% as.character(x = na.omit(object = unique(x = data[, id])))
  if (any(!groups %in% possible.clusters)) {
    stop("The following clusters were not found: ", paste(groups[!groups %in% possible.clusters], collapse = ","))
  }
  pb <- ggplot_build(plot = plot)
  if (geom == 'GeomSpatial') {
    xrange.save <- layer_scales(plot = plot)$x$range$range
    yrange.save <- layer_scales(plot = plot)$y$range$range
    data[, xynames["y"]] = max(data[, xynames["y"]]) - data[, xynames["y"]] + min(data[, xynames["y"]])
    if (!pb$plot$plot_env$crop) {
      y.transform <- c(0, nrow(x = pb$plot$plot_env$image)) - pb$layout$panel_params[[1]]$y.range
      data[, xynames["y"]] <- data[, xynames["y"]] + sum(y.transform)
    }
  }
  data <- cbind(data, color = pb$data[[1]][[1]])
  labels.loc <- lapply(
    X = groups,
    FUN = function(group) {
      data.use <- data[data[, id] == group, , drop = FALSE]
      data.medians <- if (!is.null(x = split.by)) {
        do.call(
          what = 'rbind',
          args = lapply(
            X = unique(x = data.use[, split.by]),
            FUN = function(split) {
              medians <- apply(
                X = data.use[data.use[, split.by] == split, xynames, drop = FALSE],
                MARGIN = 2,
                FUN = median,
                na.rm = TRUE
              )
              medians <- as.data.frame(x = t(x = medians))
              medians[, split.by] <- split
              return(medians)
            }
          )
        )
      } else {
        as.data.frame(x = t(x = apply(
          X = data.use[, xynames, drop = FALSE],
          MARGIN = 2,
          FUN = median,
          na.rm = TRUE
        )))
      }
      data.medians[, id] <- group
      data.medians$color <- data.use$color[1]
      return(data.medians)
    }
  )
  if (position == "nearest") {
    labels.loc <- lapply(X = labels.loc, FUN = function(x) {
      group.data <- data[as.character(x = data[, id]) == as.character(x[3]), ]
      nearest.point <- nn2(data = group.data[, 1:2], query = as.matrix(x = x[c(1,2)]), k = 1)$nn.idx
      x[1:2] <- group.data[nearest.point, 1:2]
      return(x)
    })
  }
  labels.loc <- do.call(what = 'rbind', args = labels.loc)
  labels.loc[, id] <- factor(x = labels.loc[, id], levels = levels(data[, id]))
  labels <- labels %||% groups
  if (length(x = unique(x = labels.loc[, id])) != length(x = labels)) {
    stop("Length of labels (", length(x = labels),  ") must be equal to the number of clusters being labeled (", length(x = labels.loc), ").")
  }
  names(x = labels) <- groups
  for (group in groups) {
    labels.loc[labels.loc[, id] == group, id] <- labels[group]
  }
  if (box) {
    geom.use <- ifelse(test = repel, yes = geom_label_repel, no = geom_label)
    if(is.null(label.color)) {
      plot <- plot + geom.use(
        data = labels.loc,
        mapping = aes_string(x = xynames['x'], y = xynames['y'], label = id),
        show.legend = FALSE,
        ...
      )
    } else {
      plot <- plot + geom.use(
        data = labels.loc,
        mapping = aes_string(x = xynames['x'], y = xynames['y'], label = id),
        show.legend = FALSE,
        color = label.color,
        ...
      )
    }
  } else {
    geom.use <- ifelse(test = repel, yes = geom_text_repel, no = geom_text)
    plot <- plot + geom.use(
      data = labels.loc,
      mapping = aes_string(x = xynames['x'], y = xynames['y'], label = id),
      show.legend = FALSE,
      color = "black",
      ...
    )
  }
  # restore old axis ranges
  if (geom == 'GeomSpatial') {
    plot <- suppressMessages(expr = plot + coord_fixed(xlim = xrange.save, ylim = yrange.save))
  }
  return(plot)
}

GetXYAesthetics <- function (plot, geom = "GeomPoint", plot.first = TRUE) {
  geoms <- sapply(X = plot$layers, FUN = function(layer) {
    return(class(x = layer$geom)[1])
  })
  if (geom == "GeomPoint" && "GeomScattermore" %in% geoms) {
    geom <- "GeomScattermore"
  }
  geoms <- which(x = geoms == geom)
  if (length(x = geoms) == 0) {
    stop("Cannot find a geom of class ", geom)
  }
  geoms <- min(geoms)
  if (plot.first) {
    x <- as_label(x = plot$mapping$x %||% plot$layers[[geoms]]$mapping$x)
    y <- as_label(x = plot$mapping$y %||% plot$layers[[geoms]]$mapping$y)
  } else {
    x <- as_label(x = plot$layers[[geoms]]$mapping$x %||%
                    plot$mapping$x)
    y <- as_label(x = plot$layers[[geoms]]$mapping$y %||%
                    plot$mapping$y)
  }
  return(list(x = x, y = y))
}

DimPlot2_SelColCont <- function(
  seu,
  var,
  cols,
  load.cols
) {
  list_l <- is.list(cols)
  cols_var_l <- var %in% names(cols)
  load_var <- seu@misc[["var_colors"]][[var]]
  cols_cont_l <- "continuous" %in% names(cols)
  load_cont <- seu@misc[["var_colors"]][["continuous"]]
  if(list_l & cols_var_l) {
    cols <- cols[[var]]
  } else if(load.cols & !is.null(load_var)) {
    cols <- load_var
  } else if(list_l & cols_cont_l) {
    cols <- cols[["continuous"]]
  } else if(load.cols & !is.null(load_cont)) {
    cols <- load_cont
  } else if (list_l) {
    cols <- "Blues"
  }
  scale_color <- scale_color_cont_auto(cols, center_color = FALSE, value_range = NULL)
  return(scale_color)
}

DimPlot2_SelColDisc <- function(
    seu,
    n,
    var,
    cols,
    load.cols,
    label,
    labels = waiver(),
    box
) {
  list_l <- is.list(cols)
  cols_var_l <- var %in% names(cols)
  load_var <- seu@misc[["var_colors"]][[var]]
  cols_disc_l <- "discrete" %in% names(cols)
  load_disc <- seu@misc[["var_colors"]][["discrete"]]
  label_l <- label & !box
  if(list_l & cols_var_l) {
    cols <- cols[[var]]
  } else if(load.cols & !is.null(load_var)) {
    cols <- load_var
  } else if(list_l & cols_disc_l) {
    cols <- cols[["discrete"]]
  } else if(load.cols & !is.null(load_disc)) {
    cols <- load_disc
  } else if(list_l) {
    cols <- "pro_light"
  }
  if(is.null(cols)) return(NULL)
  if(cols[1] == "light") cols <- ifelse(label_l, "pro_light","pro_light")
  scale_color <- scale_color_disc_auto(cols, n = n, labels = labels)
  return(scale_color)
}

CenterTitle <- function () {
  return(theme(plot.title = element_text(hjust = 0.5), validate = TRUE))
}

#' @title Add Simplified Axis Indicators to ggplot
#' @description Adds simplified axis indicators (arrows and labels) to a ggplot object, typically used for dimension reduction plots like UMAP or t-SNE.
#' @param anchor_x X-coordinate for the anchor point of the arrows, Default: unit(6, "mm")
#' @param anchor_y Y-coordinate for the anchor point of the arrows, Default: unit(6, "mm")
#' @param line_length Length of the arrow lines, Default: unit(12, "mm")
#' @param arrow_length Length of the arrow heads, Default: unit(2.5, "mm")
#' @param text_offset_x X-offset for the axis labels, Default: unit(2.5, "mm")
#' @param text_offset_y Y-offset for the axis labels, Default: unit(2.5, "mm")
#' @param text_size Font size for the axis labels, Default: 10
#' @param line_width Width of the arrow lines, Default: 1
#' @param x_label Custom label for the x-axis, Default: NULL (auto-detected)
#' @param y_label Custom label for the y-axis, Default: NULL (auto-detected)
#' @return A ggplot theme object that can be added to an existing ggplot
#' @details This function creates a theme that adds simplified axis indicators to the bottom-left corner of a ggplot.
#' It consists of two perpendicular arrows representing the x and y axes, along with their respective labels.
#' This approach is particularly useful for dimension reduction plots like UMAP or t-SNE where traditional axes are often removed to reduce visual clutter.
#' The function attempts to auto-detect axis labels, but custom labels can be provided.
#' @examples
#' library(Seurat)
#' library(SeuratExtend)
#'
#' features <- c("cluster", "orig.ident", "CD3D", "CD14")
#'
#' # Add arrows to the overall plot
#' DimPlot2(pbmc, features = features, theme = NoAxes()) +
#'   theme_umap_arrows()
#'
#' # Add arrows to each subplot
#' DimPlot2(pbmc, features = features, theme = theme_umap_arrows())
#' @rdname theme_umap_arrows
#' @export

theme_umap_arrows <- function(
    anchor_x = unit(6, "mm"),
    anchor_y = unit(6, "mm"),
    line_length = unit(12, "mm"),
    arrow_length = unit(2.5, "mm"),
    text_offset_x = unit(2.5, "mm"),
    text_offset_y = unit(2.5, "mm"),
    text_size = 10,
    line_width = 1,
    x_label = NULL,
    y_label = NULL
) {

  library(grid)
  ensure_unit <- function(x, default_unit = "mm") {
    if (!is(x, "unit")) {
      x <- unit(x, default_unit)
    }
    x
  }

  anchor_x <- ensure_unit(anchor_x)
  anchor_y <- ensure_unit(anchor_y)
  line_length <- ensure_unit(line_length)
  arrow_length <- ensure_unit(arrow_length)
  text_offset_x <- ensure_unit(text_offset_x)
  text_offset_y <- ensure_unit(text_offset_y)

  no_axes <- theme(
    axis.line = element_blank(),
    axis.text = element_blank(),
    axis.ticks = element_blank(),
    axis.title = element_blank()
  )

  x_arrow <- arrow(length = arrow_length, ends = "last", type = "closed")
  y_arrow <- arrow(length = arrow_length, ends = "last", type = "closed")

  structure(
    list(
      no_axes = no_axes,
      anchor_x = anchor_x,
      anchor_y = anchor_y,
      line_length = line_length,
      x_arrow = x_arrow,
      y_arrow = y_arrow,
      text_offset_x = text_offset_x,
      text_offset_y = text_offset_y,
      text_size = text_size,
      line_width = line_width,
      x_label = x_label,
      y_label = y_label
    ),
    class = "theme_umap_arrows"
  )
}

#' @exportS3Method ggplot2::ggplot_add
ggplot_add.theme_umap_arrows <- function(object, plot, object_name) {
  x_label <- object$x_label %||% extract_label(plot, "x")
  y_label <- object$y_label %||% extract_label(plot, "y")

  plot +
    object$no_axes +
    annotation_custom(
      grob = segmentsGrob(
        x0 = object$anchor_x,
        x1 = object$anchor_x + object$line_length,
        y0 = object$anchor_y,
        y1 = object$anchor_y,
        arrow = object$x_arrow,
        gp = gpar(col = "black", fill = "black", lwd = object$line_width)
      )
    ) +
    annotation_custom(
      grob = segmentsGrob(
        x0 = object$anchor_x,
        x1 = object$anchor_x,
        y0 = object$anchor_y,
        y1 = object$anchor_y + object$line_length,
        arrow = object$y_arrow,
        gp = gpar(col = "black", fill = "black", lwd = object$line_width)
      )
    ) +
    annotation_custom(
      grob = textGrob(
        label = x_label,
        x = object$anchor_x,
        y = object$anchor_y - object$text_offset_y,
        just = c(0, 1),
        gp = gpar(fontsize = object$text_size)
      )
    ) +
    annotation_custom(
      grob = textGrob(
        label = y_label,
        x = object$anchor_x - object$text_offset_x,
        y = object$anchor_y,
        just = c(0, 0),
        rot = 90,
        gp = gpar(fontsize = object$text_size)
      )
    )
}

extract_label <- function(plot, axis) {
  # First, try to extract from plot$labels
  label <- plot$labels[[axis]]

  # If the above fails, try to extract from the built plot
  if (is.null(label)) {
    # Force ggplot to compute the plot
    built_plot <- ggplot2::ggplot_build(plot)

    # Try to access grobs from the first layer
    grobs <- tryCatch({
      built_plot$plot$layers[[1]]$computed_geom_params$grob$grobs
    }, error = function(e) NULL)

    if (!is.null(grobs)) {
      for (grob in grobs) {
        if (grepl(paste0("axis.title.", axis), grob$name)) {
          # Try to extract label, with error handling
          label <- tryCatch({
            if (is.null(grob$children)) {
              grob$label
            } else if (length(grob$children) > 0) {
              grob$children[[1]]$label
            } else {
              NULL
            }
          }, error = function(e) NULL)

          if (!is.null(label)) break
        }
      }
    }
  }

  # If still not found, return the default value
  return(label %||% ifelse(axis == "x", "UMAP_1", "UMAP_2"))
}

# Helper function to apply theme elements (single element or list)
apply_theme_elements <- function(plot, theme) {
  if (is.null(theme)) {
    return(plot)
  }
  
  # If it's a list with only "list" class, treat as multiple elements
  if (is.list(theme) && identical(class(theme), "list")) {
    # Apply each theme element sequentially
    for (theme_element in theme) {
      if (!is.null(theme_element)) {
        plot <- plot + theme_element
      }
    }
  } else {
    # Single element (including theme(), labs(), NoAxes(), etc.)
    plot <- plot + theme
  }
  
  return(plot)
}
