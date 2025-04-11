#' @title Simultaneous Visualization of Three Features in a Single Plot
#' @description This function visualizes three distinct features on a single dimension reduction plot using a color blending system. It allows for the quantitative display of gene expressions or other continuous variables by mixing colors according to the RYB or RGB color models, providing a unique perspective on feature interactions and expression levels within individual cells.
#' @param seu A Seurat object that contains the data for plotting. This object should have precomputed dimensionality reduction coordinates.
#' @param feature.1 The name of the first feature (gene or other variable) to be plotted. Default: NA.
#' @param feature.2 The name of the second feature. Default: NA.
#' @param feature.3 The name of the third feature. Default: NA.
#' @param color The color model used to blend the expression data of the three features. Options include "ryb" (red-yellow-blue) and "rgb" (red-green-blue), affecting how expression intensities are represented through color. Default: c("ryb", "rgb").
#' @param color.range The range of expression intensity that is represented by the color spectrum in the plot, helping to enhance visibility of lower expressions and prevent oversaturation at high expression levels. Default: c(0.1, 0.9).
#' @param reduction The type of dimension reduction used to display the data, such as 'umap' or 'tsne'. This choice determines the underlying plot layout. Default: 'umap'.
#' @param order A logical value indicating whether to plot cells with higher expressions on top of those with lower expressions, which can help prevent significant data points from being obscured in dense areas of the plot. Default: TRUE.
#' @param pt.size Point size for plotting individual cells in the grid. Smaller values are typically used for large datasets or dense plots, whereas larger values enhance visibility for plots with fewer cells or less overlap. Default: 0.1.
#' @param dark.theme A logical value indicating whether to apply Seurat's DarkTheme to the plot. Default: FALSE.
#' @return A ggplot object that represents a dimension reduction plot incorporating three features with color blending, showing how each feature contributes to the overall expression patterns observed.
#' @details `FeaturePlot3` is designed for detailed exploratory analysis where understanding the interplay between multiple variables is crucial. This function is particularly useful for researchers looking to explore gene expressions in complex datasets, such as those involving interactions between different cell types or conditions.
#' @examples
#' library(Seurat)
#' library(SeuratExtend)
#'
#' # Basic usage
#' FeaturePlot3(
#'   pbmc,
#'   feature.1 = "CD3D",
#'   feature.2 = "CD14",
#'   feature.3 = "CD79A",
#'   color = "ryb"
#' )
#'
#' # With dark theme
#' FeaturePlot3(
#'   pbmc,
#'   feature.1 = "CD3D",
#'   feature.2 = "CD14",
#'   feature.3 = "CD79A",
#'   color = "rgb",
#'   dark.theme = TRUE,
#'   pt.size = 1,
#'   color.range = c(0.1, 1)
#' )
#'
#' @rdname FeaturePlot3
#' @export

FeaturePlot3 <- function(
  seu,
  feature.1 = NA,
  feature.2 = NA,
  feature.3 = NA,
  color = c("ryb","rgb"),
  color.range = c(0.1,0.9),
  reduction = "umap",
  order = T,
  pt.size = 0.1,
  dark.theme = FALSE
){
  library(Seurat)
  library(reshape2)
  library(ggpubr)
  library(dplyr)
  features <- c(feature.1, feature.2, feature.3)
  color <- color[1]
  if(all(is.na(features))) stop("No feature to plot.")
  tp <- matrix(rep(color.range[1], times = (ncol(seu) * 3)), ncol = 3) %>%
    as.data.frame()
  l <- color.range[1]
  h <- color.range[2]
  colors <- switch (color,
    "ryb" = c(
      ryb2rgb(c(h,l,l)),
      ryb2rgb(c(l,h,l)),
      ryb2rgb(c(l,l,h)),
      ryb2rgb(c(l,l,l))
    ),
    "rgb" = c(
      rgb(h,l,l),
      rgb(l,h,l),
      rgb(l,l,h),
      rgb(l,l,l)
    )
  )
  lgd <- function(title, col) {
    value = tp[,title]
    df <- data.frame(
      value = c(min(value),max(value))
    )
    if(max(value) == min(value)) col[2] <- col[1]
    p.tmp <-
      ggplot(df, aes(x = 1, y = 1, color = value)) +
      geom_point() +
      scale_color_gradient(low = col[1], high = col[2]) +
      labs(color = title) +
      theme(legend.justification = c(0,1))

    # Apply DarkTheme to the legend if requested
    if(dark.theme) p.tmp <- p.tmp + DarkTheme2()

    p.tmp <- get_legend(p.tmp)
    return(p.tmp)
  }
  p.leg <- list()
  for (i in 1:3) {
    if(!is.na(features[i])) {
      tp[,i] <- FetchData(seu, features[i])
      colnames(tp)[i] <- features[i]
      p.leg[[features[i]]] <- lgd(features[i], col = c(colors[4], colors[i]))
    }
  }
  tp.c <- as.data.frame(apply(tp, 2, function(x){
    if(max(x) == min(x)){
      y <- rep(color.range[1], length(x))
    } else {
      y <- sca(x, color.range)
    }
  }))
  tp.c$color <- apply(tp.c, 1, function(x){
    col <- switch(
      color,
      "ryb" = ryb2rgb(x),
      "rgb" = rgb(x[1],x[2],x[3]))
    return(col)
  })
  tp.c <- cbind(Embeddings(seu, reduction = reduction), tp.c)
  if(order) tp.c <- tp.c[order(rowMeans(tp.c[3:5])), ]
  p <-
    ggplot(tp.c, aes_string(
    x = colnames(tp.c)[1],
    y = colnames(tp.c)[2])) +
    geom_point(color = tp.c$color, size = pt.size) +
    theme_classic()

  # Apply DarkTheme if requested
  if(dark.theme) p <- p + DarkTheme2()

  p <- ggarrange(
    p, ggarrange(plotlist = p.leg, ncol = 1),
    widths = c(8,1)
  )

  if(dark.theme) p <- p + theme(plot.background = black.background)
  return(p)
}

#' @title Simultaneous Visualization of Three Features on a Grid of Plots
#' @description This function allows for the simultaneous visualization of three features across multiple plots, utilizing a grid layout. It supports two color blending systems (RYB or RGB) to represent the intensity of each feature within a single plot, providing a comprehensive overview of expression patterns across a dataset.
#' @param seu A Seurat object containing single-cell RNA sequencing data.
#' @param features A vector of feature names (genes or other continuous variables) to be displayed. This vector should be divisible by three, with each triplet of features being plotted in a separate subplot within the grid. If a triplet includes `NA`, that position will not display a feature, allowing for flexibility in visualization.
#' @param color Specifies the color blending system used to display the features. The available options are "ryb" for red-yellow-blue and "rgb" for red-green-blue. This parameter controls how the three features are visually represented based on their expression levels.
#' @param color.range Adjusts the luminance range used for feature visualization to ensure that low expressions are visible and not obscured by the background color. Default: c(0.1, 0.95), where 0.1 prevents the lowest expressions from being pure white and 0.95 keeps the highest expressions from saturating to pure color.
#' @param reduction The dimensionality reduction technique to use for the plot layout. Typically 'umap' or 'tsne' are used, with 'umap' being the default.
#' @param order Controls whether cells with higher feature expressions are plotted above those with lower expressions. This is useful for ensuring that cells with significant expression levels are visible and not obscured by those with lower levels. Default: TRUE.
#' @param pt.size Point size for plotting individual cells in the grid. Smaller values are typically used for large datasets or dense plots, whereas larger values enhance visibility for plots with fewer cells or less overlap. Default: 0.1.
#' @param legend Determines whether to display a legend describing the features and color scales. Default: FALSE.
#' @param dark.theme A logical value indicating whether to apply Seurat's DarkTheme to the plot. Default: FALSE.
#' @return A ggplot object displaying a grid of dimension reduction plots, each illustrating the expression patterns of three features using the specified color blending system.
#' @details The `FeaturePlot3.grid` function is particularly useful for exploratory data analysis where visualization of multiple gene interactions or expression patterns across different cell populations is required. It effectively combines data from multiple features into a single coherent visual representation.
#' @examples
#' library(Seurat)
#' library(SeuratExtend)
#'
#' # Basic usage
#' FeaturePlot3.grid(
#'   pbmc,
#'   features = c("CD3D", "CD14", "CD79A", "FCGR3A", NA, "LYZ"),
#'   color = "ryb",
#'   pt.size = 0.5
#' )
#'
#' # With dark theme
#' FeaturePlot3.grid(
#'   pbmc,
#'   features = c("CD3D", "CD14", "CD79A", "FCGR3A", NA, "LYZ"),
#'   color = "rgb",
#'   pt.size = 1,
#'   dark.theme = TRUE,
#'   legend = TRUE,
#'   color.range = c(0.1, 1)
#' )
#'
#' @rdname FeaturePlot3-grid
#' @export

FeaturePlot3.grid <- function(
  seu,
  features,
  color = c("ryb","rgb"),
  color.range = c(0.1,0.95),
  reduction = "umap",
  order = T,
  pt.size = 0.1,
  legend = F,
  dark.theme = FALSE
){
  library(Seurat)
  library(reshape2)
  library(ggpubr)
  library(dplyr)
  library(magrittr)
  library(rlist)
  import("ggtext")

  # features
  add3 <- function(x) {
    y <- switch(length(x) %% 3 + 1, 0, 2, 1)
    x <- c(x, rep(NA, times = y))
    return(x)
  }
  trim3 <- function(x) {
    if(length(x) > 3) {
      warning("Vector length > 3. Only use the first 3 elements")
      x <- x[1:3]
    }
    if(length(x) < 3) x <- add3(x)
    return(x)
  }
  if(is.data.frame(features)){
    ft <- features[1:3,] %>% set_rownames(NULL)
  }else if(is.list(features)){
    ft <- lapply(features, trim3)
    ft <- as.data.frame(ft)
  }else if(is.vector(features)){
    ft <- as.data.frame(matrix(add3(features), nrow = 3))
  }else stop("'features' should be class 'vector', 'data.frame' or 'list'.")

  # colors
  color <- color[1]
  l <- color.range[1]
  h <- color.range[2]
  colors <- switch (
    color,
    "ryb" = c(
      ryb2rgb(c(h,l,l)),
      ryb2rgb(c(l,h,l)),
      ryb2rgb(c(l,l,h)),
      ryb2rgb(c(l,l,l))
    ),
    "rgb" = c(
      rgb(h,l,l),
      rgb(l,h,l),
      rgb(l,l,h),
      rgb(l,l,l)
    )
  )

  # data for plot
  tp <- list()
  tp.c <- list()
  for (i in 1:ncol(ft)) {
    tp[[i]] <-
      matrix(rep(color.range[1], times = (ncol(seu) * 3)), ncol = 3) %>%
      as.data.frame()
    for (j in 1:3) {
      f <- ft[j,i]
      if(!is.na(f)) {
        tp[[i]][,j] <- FetchData(seu, f)
        colnames(tp[[i]])[j] <- f
      }
    }
    tp.c[[i]] <- as.data.frame(
      apply(tp[[i]], 2, function(x){
        if(max(x) == min(x)){
          y <- rep(color.range[1], length(x))
        } else {
          y <- sca(x, color.range)
        }
      })
    )
    tp.c[[i]]$color <- apply(
      tp.c[[i]], 1,
      function(x){
        col <- switch(
          color,
          "ryb" = ryb2rgb(x),
          "rgb" = rgb(x[1],x[2],x[3]))
        return(col)
      }
    )
    tp.c[[i]] <- cbind(Embeddings(seu, reduction = reduction), tp.c[[i]])
    if(order) tp.c[[i]] <- tp.c[[i]][order(rowMeans(tp.c[[i]][,3:5])), ]
    tp.c[[i]]$index <- i
    tp.c[[i]] <- tp.c[[i]][,-c(3:5)]
  }
  tp.c <- list.rbind(tp.c)

  # strip title
  ft.col.str <- function(f, col){
    if(is.na(f)){
      s <- NULL
    }else{
      s <- paste0(
        "<span style='color:",
        col, ";'>",
        f, "</span>")
    }
    return(s)
  }
  lab.str <-
    apply(ft, 2, function(x){
      tmp <- c(
        ft.col.str(x[1], colors[1]),
        ft.col.str(x[2], colors[2]),
        ft.col.str(x[3], colors[3])
      )
      tmp <- paste(tmp, collapse = " ")
      return(tmp)
    })
  names(lab.str) <- 1:ncol(ft)

  # main plot
  p <-
    ggplot(tp.c, aes_string(
    x = colnames(tp.c)[1],
    y = colnames(tp.c)[2])) +
    facet_wrap(~index, labeller = labeller(index = lab.str)) +
    geom_point(color = tp.c$color, size = pt.size) +
    theme_bw() +
    theme(panel.grid = element_blank(),
          axis.text = element_blank(),
          axis.ticks = element_blank(),
          strip.text = element_markdown(face = "bold", size = 11),
          strip.background = element_rect(fill = NA))

  # Apply dark theme if requested
  if(dark.theme) {
    p <- p +
      theme(
        panel.background = element_rect(fill = "black", color = "white", linewidth = 0.6),
        plot.background = black.background.no.border,
        panel.border = element_rect(color = "white"),
        strip.background = element_rect(fill = "black", color = "white", linewidth = 0.6),
        axis.title = element_text(color = "white")
      )
  }

  # figure legend
  if(legend){
    p.leg <- list()
    for (i in 1:3) {
      if(!all(is.na(ft[i,]))){
        p.tmp <-
          ggplot(data.frame(value = 0:1), aes(x = 1, y = 1, color = value)) +
          geom_point() +
          scale_color_gradient(
            low = colors[4], high = colors[i],
            breaks = c(0,1), labels = c("min", "max")) +
          labs(color = NULL) +
          theme(legend.justification = c(0,1))

        # Apply dark theme to the legend if requested
        if(dark.theme) {
          p.tmp <- p.tmp +
            theme(
              legend.background = black.background,
              legend.key = black.background.no.border,
              legend.text = element_text(color = "white")
            )
        }

        p.leg[[as.character(i)]] <- get_legend(p.tmp)
      }
    }
    p <- ggarrange(
      p, ggarrange(plotlist = p.leg, ncol = 1),
      widths = c(8,1)
    )

    # Apply dark theme to the overall combined plot
    if(dark.theme) {
      p <- p + theme(plot.background = black.background.no.border)
    }
  }
  return(p)
}

black.background <- element_rect(fill = "black")
black.background.no.border <- element_rect(fill = "black", linewidth = 0)
font.margin <- 4
white.text <- element_text(
  colour = "white", margin = margin(
    t = font.margin,
    r = font.margin,
    b = font.margin,
    l = font.margin)
)
DarkTheme2 <- function (strip.text = "white", ...) {
  white.line <- element_line(colour = "white", size = 1)
  no.line <- element_line(size = 0)
  dark.theme <- theme(
    plot.background = black.background.no.border,
    panel.background = black.background, legend.background = black.background,
    legend.box.background = black.background.no.border,
    legend.key = black.background.no.border,
    plot.title = white.text, plot.subtitle = white.text,
    axis.title = white.text, axis.text = white.text, legend.title = white.text,
    legend.text = white.text, axis.line.x = white.line,
    axis.line.y = white.line, panel.grid = no.line, panel.grid.minor = no.line,
    validate = TRUE, ...)
  if(strip.text == "white") dark.theme <- dark.theme + theme(strip.text = white.text)
  return(dark.theme)
}

