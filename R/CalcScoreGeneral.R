#' @title FUNCTION_TITLE
#' @description FUNCTION_DESCRIPTION
#' @param matr PARAM_DESCRIPTION
#' @param f PARAM_DESCRIPTION
#' @param method PARAM_DESCRIPTION, Default: 'zscore'
#' @return OUTPUT_DESCRIPTION
#' @details DETAILS
#' @examples
#' \dontrun{
#' if(interactive()){
#'  #EXAMPLE1
#'  }
#' }
#' @rdname CalcScoreGeneral
#' @export
CalcScoreGeneral<-function(matr, f, method = "zscore"){
  library(dplyr)
  library(rlist)
  library(mosaic)
  library(purrr)
  f <- factor(f)
  scores <-
    matr %>%
    t() %>%
    as.data.frame()

  if(method=="mean"){
    scores<-
      scores %>%
      split(f) %>%
      lapply(function(x) apply(x,2,mean)) %>%
      list.cbind() %>%
      as.data.frame()
  }else if(method=="median"){
    scores<-
      scores %>%
      split(f) %>%
      lapply(function(x) apply(x,2,median)) %>%
      list.cbind() %>%
      as.data.frame()
  }else if(method=="zscore"){
    scores<-
      CalcScoreGeneral(matr, f, "mean") %>%
      apply(1, zscore) %>%
      t() %>%
      as.data.frame()
  }else if(method=="tscore"){
    scores<-
      levels(f) %>%
      lapply(function(x) split(scores,f==x)) %>%
      lapply(function(x) map2(x$`TRUE`, x$`FALSE`, function(x,y) t.test(x,y)[["statistic"]])) %>%
      lapply(list.rbind) %>%
      list.cbind() %>%
      as.data.frame() %>%
      setNames(., levels(f))
  }else if(method=="p"){
    scores<-
      levels(f) %>%
      lapply(function(x) split(scores,f==x)) %>%
      lapply(function(x) map2(x$`TRUE`, x$`FALSE`, function(x,y) -log10(t.test(x,y)[["p.value"]]))) %>%
      lapply(list.rbind) %>%
      list.cbind() %>%
      as.data.frame() %>%
      setNames(., levels(f))
  }
  return(scores)
}

# CalcScoreGeneral(matr, f, "zscore") %>%
#   Heatmap()
