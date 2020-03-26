#' @title FUNCTION_TITLE
#' @description FUNCTION_DESCRIPTION
#' @param score PARAM_DESCRIPTION
#' @param color_scheme PARAM_DESCRIPTION, Default: '2'
#' @param lab_fill PARAM_DESCRIPTION, Default: 'score'
#' @param angle PARAM_DESCRIPTION, Default: 45
#' @param hjust PARAM_DESCRIPTION, Default: 1
#' @param vjust PARAM_DESCRIPTION, Default: 1
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

Heatmap<-function(score, color_scheme="2", lab_fill = "score",
                  angle = 45, hjust = 1, vjust = 1){
  library(ggplot2)
  library(reshape2)
  library(viridis)
  library(scales)

  ToPlot <-
    data.frame(score, id = factor(rownames(score), levels=unique(rev(rownames(score))))) %>%
    melt()
  p <- ggplot(ToPlot, aes(variable, id)) +
    geom_tile(aes(fill = value), colour = "white") +
    theme_classic()+
    labs(x = "", y = "", fill = lab_fill)+
    scale_y_discrete(position = "right")+
    scale_x_discrete(labels = colnames(score))+
    theme(axis.text.x=element_text(angle = angle, hjust = hjust, vjust = vjust))
  if(color_scheme == "2"){
    p <- p + scale_fill_gradient2(low = muted("blue"), high = muted("red"))
  }else if(color_scheme %in% c("A","B","C","D","E")){
    p <- p + scale_fill_gradientn(colors = viridis_pal(option = color_scheme)(20))
  }
  return(p)
}
# matr <- as.matrix(GetAssayData(seu))[1:20,]
# f <- seu@meta.data$cluster
# matr %>% CalcScoreGeneral(f, "zscore") %>% Heatmap()
# matr %>% CalcScoreGeneral(f, "zscore") %>% Heatmap(angle = -90, hjust = 0, vjust = 0.5)
