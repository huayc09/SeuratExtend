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

