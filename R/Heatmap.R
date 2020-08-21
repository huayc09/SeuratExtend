#' @title FUNCTION_TITLE
#' @description FUNCTION_DESCRIPTION
#' @param score PARAM_DESCRIPTION
#' @param color_scheme PARAM_DESCRIPTION, Default: c(low = muted("blue"), mid = "white", high = muted("red"))
#' @param border_color PARAM_DESCRIPTION, Default: 'white'
#' @param lab_fill PARAM_DESCRIPTION, Default: 'score'
#' @param angle PARAM_DESCRIPTION, Default: 45
#' @param hjust PARAM_DESCRIPTION, Default: 1
#' @param vjust PARAM_DESCRIPTION, Default: 1
#' @param legend_position PARAM_DESCRIPTION, Default: 'right'
#' @param feature_text_subset PARAM_DESCRIPTION, Default: NULL
#' @param nudge_x PARAM_DESCRIPTION, Default: ncol(score)/10
#' @param hide_axis_line PARAM_DESCRIPTION, Default: TRUE
#' @param expand_limits_x PARAM_DESCRIPTION, Default: NULL
#' @param facet_col PARAM_DESCRIPTION, Default: NULL
#' @param facet_row PARAM_DESCRIPTION, Default: NULL
#' @param panel.spacing PARAM_DESCRIPTION, Default: unit(5, "pt")
#' @param strip.placement PARAM_DESCRIPTION, Default: 'outside'
#' @return OUTPUT_DESCRIPTION
#' @details DETAILS
#' @examples
#' \dontrun{
#' if(interactive()){
#'  #EXAMPLE1
#'  }
#' }
#' @rdname Heatmap
#' @export

Heatmap <-
  function(score, color_scheme = c(low = muted("blue"), mid = "white", high = muted("red")), border_color = "white",
           lab_fill = "score", angle = 45, hjust = 1, vjust = 1, legend_position = "right",
           feature_text_subset = NULL, nudge_x = ncol(score)/10, hide_axis_line = TRUE, expand_limits_x = NULL,
           facet_col = NULL, facet_row = NULL, panel.spacing = unit(5, "pt"), strip.placement = "outside"){
  library(ggplot2)
  library(reshape2)
  library(viridis)
  library(scales)

  ToPlot <-
    data.frame(score, id = factor(rownames(score), levels=unique(rev(rownames(score))))) %>%
    melt()
  ToPlot$variable <- colnames(score)[ToPlot$variable] %>% factor(levels = unique(.))
  if(!is.null(facet_col)){
    if(length(facet_col)==ncol(score)){
      ToPlot$facet_col <- rep(facet_col, each = nrow(score))
    }else{
      stop('"facet_col" must be the same length as the number of input matrix columns')
    }
  }
  if(!is.null(facet_row)){
    if(length(facet_row)==nrow(score)){
      ToPlot$facet_row <- rep(facet_row, times = ncol(score))
    }else{
      stop('"facet_row" must be the same length as the number of input matrix rows')
    }
  }
  p <- ggplot(ToPlot, aes(variable, id)) +
    geom_tile(aes(fill = value), colour = border_color) +
    theme_classic()+
    labs(x = "", y = "", fill = lab_fill)+
    scale_y_discrete(position = "right", expand = c(0, 0)) +
    scale_x_discrete(expand = c(0, 0))+
    theme(axis.text.x=element_text(angle = angle, hjust = hjust, vjust = vjust),
          legend.position = legend_position)

  if(any(color_scheme %in% c("A","B","C","D","E"))) {
    p <- p + scale_fill_viridis(option = color_scheme)
  }else if(!is.na(color_scheme["mid"])) {
    p <- p + scale_fill_gradient2(low = color_scheme["low"],
                                  mid = color_scheme["mid"],
                                  high = color_scheme["high"])
  }else if(all(!is.na(color_scheme[c("low","high")]))) {
    p <- p + scale_fill_gradient(low = color_scheme["low"],
                                 high = color_scheme["high"])
  }else{
    p <- p + scale_fill_gradientn(colors = color_scheme)
  }

  if(hide_axis_line) {
    p <- p + theme(axis.line = element_blank(),
                   axis.ticks = element_blank())
  }

  if(!is.null(expand_limits_x)) {
    p <- p + expand_limits(x = expand_limits_x)
  }

  if(!is.null(feature_text_subset)) {
    library(ggrepel)
    p <- p +
      theme(axis.line.y = element_blank(),
            axis.text.y = element_blank(),
            axis.ticks.y = element_blank()) +
      geom_text_repel(data = data.frame(id = unique(ToPlot$id %>% .[. %in% feature_text_subset])),
                      aes(x = nlevels(ToPlot$variable)+0.5, y = id, label = id),
                      direction = "y", nudge_x = nudge_x, hjust = 0)
  }

  if(!is.null(facet_col) | !is.null(facet_row)){
    p <- p +
      facet_grid(
        rows = if(is.null(facet_row)) NULL else vars(facet_row),
        cols = if(is.null(facet_col)) NULL else vars(facet_col),
        scales = "free", space = "free",
        switch = "y") +
      theme(panel.grid = element_blank(),
            strip.background = element_rect(size = 0),
            panel.spacing = panel.spacing,
            strip.placement = "outside")
  }

  return(p)
}

# score <- CalcScoreGeneral_v3(pbmc, features = VariableFeatures(pbmc)[1:10], group.by = "seurat_clusters", method = "zscore")
# colnames(score) <- rep(1:3, times = 3)
# colnames(score) %>% class
# Heatmap(score, facet_col = rep(letters[1:3], each = 3), facet_row = c(rep("a",4), rep("b",6)), border_color = NA) +
#   theme(strip.placement = "outside")

# p2 <-
#   ggplot(ToPlot, aes(x = variable, y = 0.001, fill = facet_col)) +
#   geom_bar(stat = "identity",
#            width = 1) +
#   theme_void()+
#   theme(panel.spacing.x = unit(5, "pt")) +
#   scale_x_discrete(expand = c(0,0)) +
#   facet_grid(cols = vars(facet_col), scales = "free_x")
# legend <- plot_grid(get_legend(p), get_legend(p2), ncol = 1)
# p <- p + theme(legend.position = "none")
# p2 <- p2 + theme(legend.position = "none")
# plot_grid(p2, p, align = "v", ncol = 1, axis = "lr", rel_heights = c(0.5, 15))
