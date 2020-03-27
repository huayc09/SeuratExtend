WaterfallPlot <-
  function(matr, f, ident.1 = NULL, ident.2 = NULL,
           length = "tscore", color = "tscore", threshold = 1, flip = T){
    library(rlang)
    library(dplyr)
    library(mosaic)
    library(rlist)
    library(scales)

    f <- factor(f)
    ident.1 <- ident.1 %||% levels(f)[1]
    cell.1 <- (f == ident.1)
    cell.2 <- if(is.null(ident.2)) f != ident.1 else f == ident.2

    scores <- list()
    if("tscore" %in% c(length, color)){
      scores[["tscore"]] <-
        matr %>%
        apply(1, function(x) t.test(x[cell.1], x[cell.2])[["statistic"]]) %>%
        unlist()
    }
    scores <- scores %>%
      list.cbind() %>%
      .[abs(.[,length]) > threshold, ,drop = FALSE] %>%
      .[order(.[,length], decreasing = !flip), ,drop = FALSE] %>%
      as.data.frame() %>%
      mutate(rank = factor(rownames(.), levels = rownames(.)),
             length = .[[length]],
             color = .[[color]])

    p <- ggplot(scores, aes(x = rank, y = length, fill = color)) +
      geom_bar(stat = "identity") +
      theme_classic() +
      labs(fill = length, x = element_blank(), y = ident.1) +
      scale_fill_gradient2(low = muted("blue"), high = muted("red"))
    if(flip) p <- p + coord_flip() + scale_x_discrete(position = "top")
    return(p)
  }


# matr <- GetAssayData(seu)[1:10,]
# f <- seu$seurat_clusters
# ident.1 <- "0"
# ident.2 <- NULL
# length = "tscore"
# color = "tscore"
# threshold = 0.5
# flip = F
