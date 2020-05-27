#' @title FUNCTION_TITLE
#' @description FUNCTION_DESCRIPTION
#' @param score PARAM_DESCRIPTION
#' @param color_scheme PARAM_DESCRIPTION, Default: c(low = muted("blue"), mid = "white", high = muted("red"))
#' @param border_color Tile border color, set NA to remove the border, Default: 'white'
#' @param lab_fill PARAM_DESCRIPTION, Default: 'score'
#' @param angle PARAM_DESCRIPTION, Default: 45
#' @param hjust PARAM_DESCRIPTION, Default: 1
#' @param vjust PARAM_DESCRIPTION, Default: 1
#' @param legend_position PARAM_DESCRIPTION, Default: 'right'
#' @param feature_text_subset PARAM_DESCRIPTION, Default: NULL
#' @param nudge_x PARAM_DESCRIPTION, Default: ncol(score)/10
#' @param hide_axis_line PARAM_DESCRIPTION, Default: TRUE
#' @param expand_limits_x PARAM_DESCRIPTION, Default: NULL
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
           feature_text_subset = NULL, nudge_x = ncol(score)/10, hide_axis_line = TRUE, expand_limits_x = NULL){
  library(ggplot2)
  library(reshape2)
  library(viridis)
  library(scales)

  ToPlot <-
    data.frame(score, id = factor(rownames(score), levels=unique(rev(rownames(score))))) %>%
    melt()
  p <- ggplot(ToPlot, aes(variable, id)) +
    geom_tile(aes(fill = value), colour = border_color) +
    theme_classic()+
    labs(x = "", y = "", fill = lab_fill)+
    scale_y_discrete(position = "right")+
    scale_x_discrete(labels = colnames(score))+
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

  return(p)
}

# matr <- as.matrix(GetAssayData(seu))[1:20,]
# f <- seu@meta.data$cluster
# score <- matr %>% CalcScoreGeneral(f, "zscore")
# lab_fill = "score"
# angle = 45
# hjust = 1
# vjust = 1
# Heatmap(score)
# Heatmap(score, hide_axis_line = F)
# Heatmap(score, feature_text_subset = rownames(matr)[c(1,3,5,8)], expand_limits_x = c(0,8),
#         legend_position = "bottom")
# Heatmap(score, feature_text_subset = rownames(matr)[c(1,3,5,8)], expand_limits_x = c(0,8),
#         legend_position = "bottom")
# Heatmap(score, color_scheme = c("E"), lab_fill = "zscore")

