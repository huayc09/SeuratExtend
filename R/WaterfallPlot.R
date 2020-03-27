#' @title FUNCTION_TITLE
#' @description FUNCTION_DESCRIPTION
#' @param matr PARAM_DESCRIPTION
#' @param f PARAM_DESCRIPTION
#' @param ident.1 PARAM_DESCRIPTION, Default: NULL
#' @param ident.2 PARAM_DESCRIPTION, Default: NULL
#' @param exp.transform PARAM_DESCRIPTION, Default: F
#' @param length PARAM_DESCRIPTION, Default: 'tscore'
#' @param color PARAM_DESCRIPTION, Default: 'tscore'
#' @param len.threshold PARAM_DESCRIPTION, Default: 0
#' @param col.threshold PARAM_DESCRIPTION, Default: 0
#' @param flip PARAM_DESCRIPTION, Default: T
#' @param y.label PARAM_DESCRIPTION, Default: length
#' @param angle PARAM_DESCRIPTION, Default: NULL
#' @param hjust PARAM_DESCRIPTION, Default: NULL
#' @param vjust PARAM_DESCRIPTION, Default: NULL
#' @return OUTPUT_DESCRIPTION
#' @details DETAILS
#' @examples
#' \dontrun{
#' if(interactive()){
#'  #EXAMPLE1
#'  }
#' }
#' @rdname WaterfallPlot
#' @export

WaterfallPlot <-
  function(matr, f, ident.1 = NULL, ident.2 = NULL, exp.transform = F,
           length = "tscore", color = "tscore", len.threshold = 0, col.threshold = 0,
           flip = T, y.label = length, angle = NULL, hjust = NULL, vjust = NULL){
    library(rlang)
    library(dplyr)
    library(mosaic)
    library(rlist)
    library(scales)
    library(tidyr)

    f <- factor(f)
    ident.1 <- ident.1 %||% levels(f)[1]
    cell.1 <- (f == ident.1)
    cell.2 <- if(is.null(ident.2)) f != ident.1 else f == ident.2
    if(exp.transform) matr <- expm1(matr)

    scores <- list()
    if("tscore" %in% c(length, color)){
      scores[["tscore"]] <-
        matr %>%
        apply(1, function(x) t.test(x[cell.1], x[cell.2])[["statistic"]]) %>%
        unlist()
    }
    if("p" %in% c(length, color)){
      scores[["p"]] <-
        matr %>%
        apply(1, function(x){
          t.test(x[cell.1], x[cell.2])[["p.value"]] * ifelse(mean(x[cell.1]) > mean(x[cell.2]), 1, -1)
          }) %>%
        unlist()
    }
    if("logFC" %in% c(length, color)){
      scores[["logFC"]] <-
        matr %>%
        apply(1, function(x) log(mean(x[cell.1]+1)/mean(x[cell.2]+1))) %>%
        unlist()
    }
    scores <- scores %>%
      list.cbind() %>%
      drop_na() %>%
      .[abs(.[,length]) > len.threshold, ,drop = FALSE] %>%
      .[abs(.[,color]) > col.threshold, ,drop = FALSE] %>%
      .[order(.[,length], decreasing = !flip), ,drop = FALSE] %>%
      as.data.frame() %>%
      mutate(rank = factor(rownames(.), levels = rownames(.)),
             length = .[[length]],
             color = .[[color]])
    if(flip){
      angle <- angle %||% 0
      hjust <- hjust %||% 0
      vjust <- vjust %||% 0.5
    }else{
      angle <- angle %||% -90
      hjust <- hjust %||% 0
      vjust <- vjust %||% 0.5
    }
    p <- ggplot(scores, aes(x = rank, y = length, fill = color)) +
      geom_bar(stat = "identity") +
      theme_classic() +
      labs(fill = color, x = element_blank(), y = y.label) +
      scale_fill_gradient2(low = muted("blue"), high = muted("red")) +
      theme(axis.text.x=element_text(angle = angle, hjust = hjust, vjust = vjust))
    if(flip) p <- p + coord_flip() + scale_x_discrete(position = "top")
    return(p)
  }

#' @title FUNCTION_TITLE
#' @description FUNCTION_DESCRIPTION
#' @param Seu PARAM_DESCRIPTION
#' @param features PARAM_DESCRIPTION
#' @param group.by PARAM_DESCRIPTION
#' @param assay PARAM_DESCRIPTION, Default: 'RNA'
#' @param ident.1 PARAM_DESCRIPTION, Default: NULL
#' @param ident.2 PARAM_DESCRIPTION, Default: NULL
#' @param exp.transform PARAM_DESCRIPTION, Default: T
#' @param length PARAM_DESCRIPTION, Default: 'logFC'
#' @param color PARAM_DESCRIPTION, Default: 'p'
#' @param len.threshold PARAM_DESCRIPTION, Default: 0
#' @param col.threshold PARAM_DESCRIPTION, Default: 0
#' @param flip PARAM_DESCRIPTION, Default: F
#' @param y.label PARAM_DESCRIPTION, Default: length
#' @param angle PARAM_DESCRIPTION, Default: NULL
#' @param hjust PARAM_DESCRIPTION, Default: NULL
#' @param vjust PARAM_DESCRIPTION, Default: NULL
#' @return OUTPUT_DESCRIPTION
#' @details DETAILS
#' @examples
#' \dontrun{
#' if(interactive()){
#'  #EXAMPLE1
#'  }
#' }
#' @rdname WaterfallPlot_v3
#' @export

WaterfallPlot_v3 <-
  function(Seu, features, group.by, assay = "RNA", ident.1 = NULL, ident.2 = NULL, exp.transform = T,
           length = "logFC", color = "p", len.threshold = 0, col.threshold = 0,
           flip = F, y.label = length, angle = NULL, hjust = NULL, vjust = NULL){
    library(Seurat)
    if(any(!features %in% rownames(Seu))) {
      WrongGenes <- setdiff(features, rownames(Seu))
      warning(paste0(paste(WrongGenes, collapse = ", "), " are not detected in the matrix"))
      features <- intersect(features, rownames(Seu))
    }
    matr <- GetAssayData(Seu, assay = assay)[features,]
    f <- Seu@meta.data[, group.by]
    return(WaterfallPlot(matr, f, ident.1, ident.2, exp.transform,
                         length, color, len.threshold, col.threshold,
                         flip, y.label, angle, hjust, vjust))
  }
# matr <- GetAssayData(seu)[1:10,]
# f <- seu$seurat_clusters
# ident.1 <- "0"
# ident.2 <- NULL
# length = "logFC"
# color = "p"
# threshold = 0.5
# flip = F
# exp.transform = T
# WaterfallPlot(matr, seu$seurat_clusters,
#               length = "logFC", color = "p", exp.transform = T, flip = F)
# WaterfallPlot_v3(seu, group.by = "seurat_clusters", features = rownames(seu)[1:10])

