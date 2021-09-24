#' @title FUNCTION_TITLE
#' @description FUNCTION_DESCRIPTION
#' @param seu PARAM_DESCRIPTION
#' @param feature.1 PARAM_DESCRIPTION, Default: NA
#' @param feature.2 PARAM_DESCRIPTION, Default: NA
#' @param feature.3 PARAM_DESCRIPTION, Default: NA
#' @param color PARAM_DESCRIPTION, Default: c("ryb", "rgb")
#' @param color.range PARAM_DESCRIPTION, Default: c(0.1, 0.9)
#' @param reduction PARAM_DESCRIPTION, Default: 'umap'
#' @param order PARAM_DESCRIPTION, Default: T
#' @return OUTPUT_DESCRIPTION
#' @details DETAILS
#' @examples
#' \dontrun{
#' if(interactive()){
#'  #EXAMPLE1
#'  }
#' }
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
  order = T
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
    geom_point(color = tp.c$color) +
    theme_classic()
  p <- ggarrange(
    p, ggarrange(plotlist = p.leg, ncol = 1),
    widths = c(8,1)
  )
  return(p)
}

#' @title FUNCTION_TITLE
#' @description FUNCTION_DESCRIPTION
#' @param seu PARAM_DESCRIPTION
#' @param features PARAM_DESCRIPTION
#' @param color PARAM_DESCRIPTION, Default: c("ryb", "rgb")
#' @param color.range PARAM_DESCRIPTION, Default: c(0.1, 0.95)
#' @param reduction PARAM_DESCRIPTION, Default: 'umap'
#' @param order PARAM_DESCRIPTION, Default: T
#' @param pt.size PARAM_DESCRIPTION, Default: 0.1
#' @param legend PARAM_DESCRIPTION, Default: F
#' @return OUTPUT_DESCRIPTION
#' @details DETAILS
#' @examples
#' \dontrun{
#' if(interactive()){
#'  #EXAMPLE1
#'  }
#' }
#' @rdname FeaturePlot3.grid
#' @export
FeaturePlot3.grid <- function(
  seu,
  features,
  color = c("ryb","rgb"),
  color.range = c(0.1,0.95),
  reduction = "umap",
  order = T,
  pt.size = 0.1,
  legend = F
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
        p.leg[[as.character(i)]] <- get_legend(p.tmp)
      }
    }
    p <- ggarrange(
      p, ggarrange(plotlist = p.leg, ncol = 1),
      widths = c(8,1)
    )
  }
  return(p)
}


