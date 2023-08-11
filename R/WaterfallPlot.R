#' @include generics.R
#'
NULL

#' @param seu A Seurat object. Only applicable when using the Seurat method.
#' @param features Features to be plotted, which can include gene expression, metrics, PC scores, or any other data that can be retrieved using the `FetchData()` function. Only applicable for the Seurat method.
#' @param group.by A variable from `meta.data` for grouping, or a character vector of the same length as the number of cells. Only applicable for the Seurat method.
#' @param cell Cell identifiers to be used in the plot. Defaults to all cells. Only applicable for the Seurat method.
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
    top.n = NULL,
    flip = FALSE,
    y.label = length,
    angle = NULL,
    hjust = NULL,
    vjust = NULL,
    title = NULL
) {

  Std.matr <- Seu2Matr(
    seu = seu,
    features = features,
    group.by = group.by,
    cells = cell,
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
    top.n = top.n,
    flip = flip,
    y.label = y.label,
    angle = angle,
    hjust = hjust,
    vjust = vjust,
    title = title
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
#' @param top.n Retains only the top `n` bars in both positive and negative directions. If `length(top.n)` is 1, the function retains `top.n` bars for both positive and negative directions. If `length(top.n)` is 2, it retains `top.n[1]` positive bars and `top.n[2]` negative bars. Defaults to NULL.
#' @param flip Determines whether the plot should be flipped. Defaults to TRUE for matrix input and FALSE for Seurat object input.
#' @param y.label Label for the y-axis. Defaults to "length".
#' @param angle Angle of the x-axis labels. This argument is passed to `element_text()`.
#' @param hjust Horizontal justification for the x-axis labels. This argument is passed to `element_text()`.
#' @param vjust Vertical justification for the x-axis labels. This argument is passed to `element_text()`.
#' @param title Title of the plot. Defaults to NULL.
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
    top.n = NULL,
    flip = TRUE,
    y.label = length,
    angle = NULL,
    hjust = NULL,
    vjust = NULL,
    title = NULL
){

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
    top.n = top.n
  )

  if(is.null(title)) {
    title <- WaterfallPlot_Title(
      f = f,
      ident.1 = ident.1,
      ident.2 = ident.2)
  }

  p <- WaterfallPlot_Plot(
    scores = scores,
    color = color,
    flip = flip,
    y.label = y.label,
    angle = angle,
    hjust = hjust,
    vjust = vjust,
    title = title
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
    top.n
) {
  library(dplyr)
  library(rlist)
  library(tidyr)

  f <- factor(f)
  ident.1 <- ident.1 %||% levels(f)[1]
  cell.1 <- (f == ident.1)
  cell.2 <- if(is.null(ident.2)) f != ident.1 else f == ident.2
  if(exp.transform) matr <- expm1(matr)

  scores <- list()
  if("tscore" %in% c(length, color)){
    scores[["tscore"]] <-
      apply(matr, 1, function(x) t.test(x[cell.1], x[cell.2])[["statistic"]], simplify = TRUE)
  }
  if("p" %in% c(length, color)){
    scores[["p"]] <-
      apply(matr, 1, function(x){
        -log10(t.test(x[cell.1], x[cell.2])[["p.value"]]) * ifelse(mean(x[cell.1]) > mean(x[cell.2]), 1, -1)
      }, simplify = TRUE)
  }
  if("logFC" %in% c(length, color)){
    scores[["logFC"]] <-
      apply(matr, 1, function(x) log(mean(x[cell.1]+1)/mean(x[cell.2]+1)), simplify = TRUE)
  }
  # scores <- scores %>%
  #   list.cbind() %>%
  #   as.data.frame() %>%
  #   drop_na() %>%
  #   .[abs(.[,length]) > len.threshold, ,drop = FALSE] %>%
  #   .[abs(.[,color]) > col.threshold, ,drop = FALSE]
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
  # scores <- scores %>%
  #   as.data.frame() %>%
  #   mutate(rank = factor(rownames(.), levels = rownames(.)),
  #          length = .[[length]],
  #          color = .[[color]])
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
    flip,
    y.label,
    angle,
    hjust,
    vjust,
    title
) {
  library(ggplot2)
  library(scales)
  if(flip){
    scores <- scores[nrow(scores):1, ]
    angle <- angle %||% 0
    hjust <- hjust %||% 0
    vjust <- vjust %||% 0.5
  }else{
    angle <- angle %||% -90
    hjust <- hjust %||% 0
    vjust <- vjust %||% 0.5
  }
  scores$rank <- factor(scores$rank, levels = unique(scores$rank))
  if(color == "p") lab_fill <- "-log10(p)" else lab_fill <- color
  p <- ggplot(scores, aes(x = rank, y = length, fill = color)) +
    geom_bar(stat = "identity") +
    theme_classic() +
    labs(fill = lab_fill, x = element_blank(), y = y.label, title = title) +
    scale_fill_gradient2(low = muted("blue"), high = muted("red")) +
    theme(axis.text.x=element_text(angle = angle, hjust = hjust, vjust = vjust),
          plot.title = element_text(hjust = 0.5, face = "bold"))
  if(flip) p <- p + scale_x_discrete(position = "top") + coord_flip()
  return(p)
}

WaterfallPlot_Title <- function(f, ident.1, ident.2) {
  f <- factor(f)
  ident.1 <- ident.1 %||% levels(f)[1]
  ident.2 <- ident.2 %||% paste0("non-", ident.1)
  title <- paste0(ident.1, " vs. ", ident.2)
  return(title)
}
