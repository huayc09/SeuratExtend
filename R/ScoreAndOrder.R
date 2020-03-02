#' @title FUNCTION_TITLE
#' @description FUNCTION_DESCRIPTION
#' @param matr PARAM_DESCRIPTION
#' @param f PARAM_DESCRIPTION
#' @param method PARAM_DESCRIPTION, Default: 'zscore'
#' @param order PARAM_DESCRIPTION, Default: 'p'
#' @param n PARAM_DESCRIPTION, Default: Inf
#' @param p.threshold PARAM_DESCRIPTION, Default: 0.05
#' @return OUTPUT_DESCRIPTION
#' @details DETAILS
#' @examples
#' \dontrun{
#' if(interactive()){
#'  #EXAMPLE1
#'  }
#' }
#' @rdname ScoreAndOrder
#' @export
ScoreAndOrder<-function(matr, f, method = "zscore", order = "p", n = Inf, p.threshold = 0.05){
  library(purrr)
  f <- factor(f)
  score <- CalcScoreGeneral(matr, f, method)
  if(order=="value") {
    score <-
      score %>%
      split(factor(levels(f)[apply(., 1, which.max)], levels = levels(f))) %>%
      map2(names(.), function(x, i) x[order(x[[i]], decreasing = T),] %>% head(n)) %>%
      setNames(NULL) %>%
      list.rbind()
  } else if(order=="p") {
    score <-
      score %>%
      split(factor(levels(f)[apply(., 1, which.max)], levels = levels(f))) %>%
      map2(names(.), function(x, i) {
        as.data.frame(matr)[rownames(x),] %>%
          apply(1, function(y) t.test(y[f==i], y[f!=i])[["p.value"]]) %>%
          .[order(.)] %>%
          .[.<p.threshold] %>%
          names() %>%
          x[.,] %>%
          head(n)
      }) %>%
      setNames(NULL) %>%
      list.rbind()
  }
  return(score)
}
