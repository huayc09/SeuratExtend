#' @include generics.R
#'
NULL

#' @param seu (Seurat version) Seurat object
#' @param group.by (Seurat version) A variable name in meta.data to
#' group the violin plots by
#' @param split.by (Seurat version) A variable name in meta.data to
#' split the violin plots by
#' @param cell (Seurat version) Cell names to use, Default: all cells
#' @rdname VlnPlot2
#' @export

VlnPlot2.Seurat <- function(
  seu,
  features,
  group.by = NULL,
  split.by = NULL,
  cell = NULL,
  ...
) {
  require(rlang)
  cell <- cell %||% colnames(seu)
  matr <- t(FetchData(seu, vars = features, cells = cell))
  if(is.null(group.by)) {
    f <- factor(Idents(seu)[cell])
  }else{
    f <- factor(seu[[group.by]][cell,])
  }
  f2 <- seu[[split.by]][cell,]
  p <- VlnPlot2.default(matr, f, f2, features, ...)
  return(p)
}

#' @param matr Matrix or data frame.Row - features; columns - cells
#' @param f Factor or vector. Identity of each cell. Should be the
#' same length of cells
#' @param f2 Factor or vector. Similar to \code{f}. A variable to split
#' the violin plots by. Default: NULL
#' @param features Features to plot (gene expression, metrics, PC scores,
#' anything that can be retreived by FetchData), Default: NULL (All features
#' in matrix)
#' @param ncol Number of columns if multiple plots are displayed, Default: 1
#' @param lab_fill Title of figure legend, Default: 'group'
#' @param scales scales parameter passed to \code{\link[ggplot2:facet_wrap]{ggplot2::facet_wrap()}}
#' , Default: 'free_y'
#' @param violin Whether to plot violin plot, Default: T
#' @param box Whether to plot box plot, Default: T
#' @param width Width of box plot, Default: 0.9
#' @param pt Whether to plot points, Default: T
#' @param hide.outlier Whether to hide outlier points of boxplot, Default: F
#' @param pt.style Position adjustment, Default: c("jitter", "quasirandom")
#' @param pt.size Point size, Default: 1
#' @param pt.alpha Point transparency, Default: 0.35
#' @param strip.position Were to put the strip ("top", "bottom", "left" (default)
#' or "right"). Only use when \code{f2 = NULL}
#' @rdname VlnPlot2
#' @export

VlnPlot2.default <- function(
  matr, f, f2 = NULL,
  features = NULL,
  ncol = 1,
  lab_fill = "group",
  scales = "free_y",
  violin = T,
  box = T,
  width = 0.9,
  pt = T,
  hide.outlier = F,
  pt.style = c("jitter","quasirandom"),
  pt.size = 1,
  pt.alpha = 0.35,
  strip.position = "left"
) {
  library(ggplot2)
  library(rlang)
  library(dplyr)
  library(reshape2)
  features <- features %||% rownames(matr)
  if(!is_empty(setdiff(features, rownames(matr)))){
    message(paste0(setdiff(features, rownames(matr)), collapse = ", "), " not found")
    features <- intersect(features, rownames(matr))
  }
  f2 <- f2 %||% data.frame(row.names = colnames(matr))
  ToPlot <-
    cbind(f, f2, as.data.frame(t(matr[features,,drop = F]))) %>%
    melt(measure.vars = features)
  if(!violin & !box & !pt) stop("No plot type defined (violin/boxplot/point)")
  x <- ifelse(is_empty(f2), "f", "f2")
  p <- ggplot(ToPlot, aes_string(x = x, y = "value", fill = x))
  if(violin) {
    p <- p +
      geom_violin(scale = "width", width = width)
  }
  if(box & !violin) {
    if(pt | hide.outlier) {
      p <- p +
        geom_boxplot(outlier.shape = NA, width = width)
    } else {
      p <- p +
        geom_boxplot(outlier.size = pt.size, width = width)
    }
  }
  if(pt) {
    pt.style <- pt.style[1]
    if(!pt.style %in% c("quasirandom", "jitter")) stop('"pt.style" shoule be "quasirandom" or "jitter"')
    if(pt.style == "jitter") p <- p + geom_jitter(width = width/2, size = pt.size, alpha= pt.alpha)
    if(pt.style == "quasirandom") {
      import("ggbeeswarm")
      p <- p + geom_quasirandom(size = pt.size, width = width/2, alpha= pt.alpha)
    }
  }
  if(box & violin) {
    if(pt | hide.outlier) {
      p <- p +
        geom_boxplot(outlier.shape = NA, width = 0.1, fill = "white")
    } else {
      p <- p +
        geom_boxplot(fill = "white", outlier.size = pt.size, width = 0.1, outlier.alpha = pt.alpha)
    }
  }
  if(is_empty(f2)){
    p <- p +
      facet_wrap( ~variable, ncol = ncol, strip.position=strip.position, scales = scales)+
      ylab(NULL) +
      xlab(NULL) +
      theme_classic() +
      theme(strip.background = element_blank(),
            strip.placement = "outside",
            legend.position = "none",
            axis.text.x = element_text(angle = 45,hjust = 1),
            strip.text = element_text(face = "bold", size = 10)) +
      labs(fill = lab_fill) +
      scale_y_continuous(expand = expansion(mult = c(0,0.05)))
  }else{
    p <- p +
      facet_grid(vars(variable), vars(f), switch = c("both"), scales = scales)+
      ylab(NULL) +
      xlab(NULL) +
      theme_classic() +
      theme(strip.background = element_blank(),
            strip.placement = "outside",
            axis.text.x = element_blank(),
            strip.text.y = element_text(face = "bold", size = 10)) +
      labs(fill = lab_fill) +
      scale_y_continuous(expand = expansion(mult = c(0,0.05)))
  }
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



