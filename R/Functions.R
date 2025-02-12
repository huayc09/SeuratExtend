#' @title Proportion of features
#' @description Check the Proportion of positive cells (default: expression above 0)
#' in certain clusters
#' @param seu Seurat object
#' @param feature Features to plot (gene expression, metrics, PC scores, anything that can be
#' retreived by FetchData)
#' @param ident cluster name, Default: all clusters
#' @param group.by A variable name in meta.data to group by, if you don't want to use
#' default Idents of Seurat object
#' @param DefaultAssay Assay to use, Default: NULL (uses object's default assay)
#' @param above Above which value will be considered as positive cell, Default: 0
#' @param total If to calculate proportion in total cells of selected clusters, Default: F
#' @param if.expressed If only return logical value, Default: F
#' @param min.pct Minimum fraction of min.pct cells in the cluster will be considered "expressed",
#' only use this parameter when "if.expressed" is set to TRUE. Default: 0.1
#' @return A matrix
#' @details None.
#' @examples
#' genes <- VariableFeatures(pbmc)[1:10]
#'
#' feature_percent(pbmc, feature = genes)
#'
#' # Count one cell as positive only the expression is above 2
#' feature_percent(pbmc, feature = genes, above = 2)
#'
#' # Only check subset of clusters
#' feature_percent(pbmc, feature = genes, ident = c("B cell", "CD8 T cell"))
#'
#' # Group by a different variable
#' feature_percent(pbmc, feature = genes, group.by = "orig.ident")
#'
#' # Also check the proportion of expressed cells in total clusters
#' feature_percent(pbmc, feature = genes, total = T)
#'
#' # only show logical value; 'TRUE' if more than 20% cells are positive
#' feature_percent(pbmc, feature = genes, if.expressed = T, min.pct = 0.2)
#'
#' @rdname feature_percent
#' @export

feature_percent <- function(
  seu,
  feature,
  ident = NULL,
  group.by = NULL,
  DefaultAssay = NULL,
  above = 0,
  total = F,
  if.expressed = F,
  min.pct = 0.1
) {
  library(Seurat)
  library(rlist)

  # Handle default assay and warning
  if (is.null(DefaultAssay)) {
    DefaultAssay <- DefaultAssay(seu)
  }
  if (DefaultAssay == "TF") {
    warning("Current assay is set to 'TF'. If you want to examine gene expression, consider changing to 'RNA' assay.")
  }

  DefaultAssay(seu) <- DefaultAssay

  if(is.null(group.by)){
    f <- Idents(seu)
  } else {
    f <- factor(seu@meta.data[,group.by])
  }
  ident <- ident %||% levels(f)
  cells <- colnames(seu)[f %in% ident]
  m <- FetchData(seu, vars = feature, cells = cells)
  f2 <- f[f %in% ident]
  f2 <- factor(f2)
  m2 <- split(m, f2)
  if(total) m2[["total"]] <- m
  m2 <- lapply(m2, function(x){
    apply(x, 2, function(y) sum(y > above)) / nrow(x)
  })
  m2 <- list.cbind(m2)
  if(if.expressed) m2 <- m2 > min.pct
  return(m2)
}

#' Display Multiple Color Palettes in a Grid Layout
#'
#' @description
#' Creates a grid visualization of multiple color palettes using ggplot2. Each palette is displayed
#' as a row of colored tiles with white borders, making it easy to compare different color schemes
#' or visualize multiple variations of a palette.
#'
#' @param palette_list A named list of color palettes. Each element should be a character vector
#'        of colors (hex codes or color names).
#' @param ncol Integer specifying the number of columns in the grid. Determines how many colors
#'        are displayed per row before wrapping.
#'
#' @return A ggplot object displaying the color palettes in a grid format.
#'
#' @details
#' The function is particularly useful for visualizing:
#' - Different variations of color_pro() outputs
#' - Multiple color schemes side by side
#' - Comparisons between different palette generation methods
#' Each row is labeled with the name from the input list, and colors are displayed
#' as tiles with white borders for clear separation.
#'
#' @examples
#' \dontrun{
#' library(SeuratExtend)
#'
#' # Example 1: Compare different numbers of colors
#' palette_list <- list(
#'   "n=2" = color_pro(2),
#'   "n=5" = color_pro(5),
#'   "n=10" = color_pro(10),
#'   "n=20" = color_pro(20)
#' )
#' show_col2(palette_list, ncol = 20)
#'
#' # Example 2: Compare different color schemes
#' palette_list <- list(
#'   "default" = color_pro(10, 1),
#'   "light" = color_pro(10, 2),
#'   "red" = color_pro(10, 3),
#'   "yellow" = color_pro(10, 4),
#'   "green" = color_pro(10, 5),
#'   "blue" = color_pro(10, 6),
#'   "purple" = color_pro(10, 7)
#' )
#' show_col2(palette_list, ncol = 10)
#'
#' # Example 3: Compare different sorting methods
#' palette_list <- list(
#'   "sort_hue" = color_pro(10, 1, sort = "hue"),
#'   "sort_diff" = color_pro(10, 1, sort = "diff")
#' )
#' show_col2(palette_list, ncol = 10)
#'
#' # Example 4: Compare different random sets
#' palette_list <- list(
#'   "Set1" = color_pro(10, 1, 1, 1),
#'   "Set2" = color_pro(10, 1, 1, 2),
#'   "Set3" = color_pro(10, 1, 1, 3),
#'   "Set4" = color_pro(10, 1, 1, 4),
#'   "Set5" = color_pro(10, 1, 1, 5)
#' )
#' show_col2(palette_list, ncol = 10)
#'
#' # Example 5: Compare color_iwh styles
#' palette_list <- list(
#'   "default" = color_iwh(10, 1),
#'   "intense" = color_iwh(10, 2),
#'   "pastel" = color_iwh(10, 3),
#'   "all" = color_iwh(10, 4),
#'   "all_hard" = color_iwh(30, 5)
#' )
#' show_col2(palette_list, ncol = 10)
#' }
#'
#' @import ggplot2
#' @importFrom reshape2 melt
#' @export

show_col2 <- function(palette_list, ncol) {
  library(ggplot2)
  pal_mtx_list <- lapply(palette_list, function(x) {
    if (length(x) < ncol) x <- x[1:ncol]
    return(matrix(x, ncol = ncol))
  })
  facet_row <-
    lapply(names(pal_mtx_list), function(x) {
      rep(x, nrow(pal_mtx_list[[x]]))
    }) %>% unlist
  pal_mtx <- list.rbind(pal_mtx_list)
  ToPlot <- melt(pal_mtx)
  ToPlot$facet_row <- factor(rep(facet_row, times = ncol), levels = unique(facet_row))
  p <-
    ggplot(ToPlot, aes(Var2, Var1)) +
    geom_tile(aes(fill = value), colour = "white", linewidth = 1) +
    theme_classic()+
    scale_fill_identity() + NoAxes() +
    scale_y_reverse() +
    theme(strip.background = element_blank(),
          strip.text.y.left = element_text(angle = 0, size = 12)) +
    facet_grid(rows = vars(facet_row), scales = "free", space = "free", switch = "y")
  return(p)
}


