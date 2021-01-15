#' @title FUNCTION_TITLE
#' @description FUNCTION_DESCRIPTION
#' @param score PARAM_DESCRIPTION
#' @param color_scheme PARAM_DESCRIPTION, Default: c(low = muted("blue"), mid = "white", high = muted("red"))
#' @param border_color PARAM_DESCRIPTION, Default: NULL
#' @param lab_fill PARAM_DESCRIPTION, Default: 'score'
#' @param angle PARAM_DESCRIPTION, Default: 45
#' @param hjust PARAM_DESCRIPTION, Default: 1
#' @param vjust PARAM_DESCRIPTION, Default: 1
#' @param legend_position PARAM_DESCRIPTION, Default: 'right'
#' @param feature_text_subset PARAM_DESCRIPTION, Default: NULL
#' @param segment.width PARAM_DESCRIPTION, Default: c(1, 2.5, 1)
#' @param segment.size PARAM_DESCRIPTION, Default: 0.2
#' @param text.spacing PARAM_DESCRIPTION, Default: 2
#' @param text.size PARAM_DESCRIPTION, Default: 2.5
#' @param hide_axis_line PARAM_DESCRIPTION, Default: TRUE
#' @param expand_limits_x PARAM_DESCRIPTION, Default: NULL
#' @param facet_col PARAM_DESCRIPTION, Default: NULL
#' @param facet_row PARAM_DESCRIPTION, Default: NULL
#' @param panel.spacing PARAM_DESCRIPTION, Default: unit(5, "pt")
#' @param strip.placement PARAM_DESCRIPTION, Default: 'outside'
#' @param ncol PARAM_DESCRIPTION, Default: NULL
#' @param nrow PARAM_DESCRIPTION, Default: NULL
#' @param ... PARAM_DESCRIPTION
#' @return OUTPUT_DESCRIPTION
#' @details DETAILS
#' @examples
#' \dontrun{
#' if(interactive()){
#'  #EXAMPLE1
#'  }
#' }
#' @seealso
#'  \code{\link[grid]{grid.draw}}
#' @rdname Heatmap
#' @export
#' @importFrom grid grid.draw

Heatmap <-
  function(score, color_scheme = c(low = muted("blue"), mid = "white", high = muted("red")), border_color = NULL,
           lab_fill = "score", angle = 45, hjust = 1, vjust = 1, legend_position = "right",
           feature_text_subset = NULL, segment.width = c(1,2.5,1), segment.size = 0.2, text.spacing = 2, text.size = 2.5,
           hide_axis_line = TRUE, expand_limits_x = NULL,
           facet_col = NULL, facet_row = NULL, panel.spacing = unit(5, "pt"), strip.placement = "outside",
           ncol = NULL, nrow = NULL, ...){
    library(ggplot2)
    library(reshape2)
    library(viridis)
    library(scales)
    library(rlang)
    library(dplyr)

    ToPlot <-
      data.frame(score, id = factor(rownames(score), levels=unique(rev(rownames(score))))) %>%
      melt()
    ToPlot$variable <- colnames(score)[ToPlot$variable] %>% factor(levels = unique(.))
    if(!is.null(facet_col)){
      if(length(facet_col)==ncol(score)){
        if(!is.null(feature_text_subset)) facet_col <- factor(facet_col) %>% `levels<-`(c(levels(.),""))
        ToPlot$facet_col <- rep(facet_col, each = nrow(score))
      }else{
        stop('"facet_col" must be the same length as the number of input matrix columns')
      }
    }
    if(!is.null(facet_row)){
      if(!is.null(feature_text_subset)) stop('"facet_row" and "feature_text_subset" cannot be set at the same time')
      if(length(facet_row)==nrow(score)){
        ToPlot$facet_row <- rep(facet_row, times = ncol(score))
      }else{
        stop('"facet_row" must be the same length as the number of input matrix rows')
      }
    }
    border_color <- border_color %||% ifelse(ncol(score) > 80 | nrow(score) > 80, NA, "white")

    p <- ggplot(ToPlot, aes(variable, id)) +
      geom_tile(aes(fill = value), colour = border_color) +
      theme_classic()+
      labs(x = "", y = "", fill = lab_fill)+
      scale_y_discrete(position = "right", expand = c(0, 0)) +
      scale_x_discrete(expand = c(0, 0))+
      theme(axis.text.x=element_text(angle = angle, hjust = hjust, vjust = vjust),
            legend.position = legend_position) +
      theme(...)

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
      text.spacing = nrow(score) * text.spacing / 100
      text_table <-
        data.frame(text = levels(ToPlot$id) %>% .[. %in% feature_text_subset],
                   y.orig = which(levels(ToPlot$id) %in% feature_text_subset))
      text_table[nrow(text_table),"y1"] <- text_table[nrow(text_table),"y.orig"]
      if(text_table[nrow(text_table),"y1"] > nrow(score) - text.spacing / 3) {
        text_table[nrow(text_table),"y1"] <- nrow(score) - text.spacing / 3
      }
      for (i in (nrow(text_table)-1):1) {
        if(text_table[i,"y.orig"] > text_table[(i+1),"y1"] - text.spacing) {
          text_table[i,"y1"] <- text_table[(i+1),"y1"] - text.spacing
        } else {
          text_table[i,"y1"] <- text_table[i,"y.orig"]
        }
      }
      if(text_table[1, "y1"] < 1 + text.spacing / 3) {
        i <- 1
        text_table[1, "y1"] <- 1 + text.spacing / 3
        while (text_table[(i+1), "y1"] < text_table[i, "y1"] + text.spacing) {
          text_table[(i+1), "y1"] <- text_table[i, "y1"] + text.spacing
          i <- i + 1
        }
      }
      if(!is.null(facet_col)) {
        last.level <- factor("", levels = levels(facet_col))
        text_table$facet_col <- last.level
        x0 <- 0
      } else x0 <- ncol(score) + 0.5 + ncol(score)*0.01

      text_table$x1 <- x0
      text_table$x2 <- text_table$x1 + ncol(score) * segment.width[1] / 100
      text_table$x3 <- text_table$x2 + ncol(score) * segment.width[2] / 100
      text_table$x4 <- text_table$x3 + ncol(score) * segment.width[3] / 100


      p <- p +
        theme(axis.line.y = element_blank(),
              axis.text.y = element_blank(),
              axis.ticks.y = element_blank()) +
        geom_segment(data = text_table, aes(x = x1, y = y.orig, xend = x2, yend = y.orig), lineend = "round", size = segment.size) +
        geom_segment(data = text_table, aes(x = x2, y = y.orig, xend = x3, yend = y1), lineend = "round", size = segment.size) +
        geom_segment(data = text_table, aes(x = x3, y = y1, xend = x4, yend = y1), lineend = "round", size = segment.size) +
        geom_text(data = text_table, aes(label = text, x = x4 + ncol(score)*0.01, y = y1), hjust = 0, size = text.size)
    }

    if(!is.null(facet_col) | !is.null(facet_row)){
      if(!is.null(ncol) | !is.null(nrow)) {
        p <- p +
          facet_rep_wrap(
            facets = vars(facet_col), ncol = ncol, nrow = nrow,
            scales = "fixed")
      } else {
        p <- p +
          facet_grid(
            rows = if(is.null(facet_row)) NULL else vars(facet_row),
            cols = if(is.null(facet_col)) NULL else vars(facet_col),
            scales = "free", space = "free",
            switch = "y")
      }
      p <- p +
        theme(panel.grid = element_blank(),
              strip.background = element_rect(size = 0),
              panel.spacing = panel.spacing,
              strip.placement = strip.placement)
      if(!is.null(feature_text_subset)){
        g <- ggplot_gtable(ggplot_build(p))
        stripr <- which(grepl('strip-r', g$layout$name) | grepl('strip-t', g$layout$name)) %>% tail(1)
        j <- which(grepl('rect', g$grobs[[stripr]]$grobs[[1]]$childrenOrder))
        g$grobs[[stripr]]$grobs[[1]]$children[[j]]$gp$col <- NA
        p <- grid::grid.draw(g)
      }
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

# Heatmap(matrix(sample(-5000:4999), ncol = 10) %>% as.data.frame(),
#         feature_text_subset = c(1:10, 30:40, 990:999), facet_col = c("Z",1,1,2,2,2,3,3,3,"Z"), expand_limits_x = c(1))

