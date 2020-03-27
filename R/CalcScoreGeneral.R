#' @title FUNCTION_TITLE
#' @description FUNCTION_DESCRIPTION
#' @param matr PARAM_DESCRIPTION
#' @param f PARAM_DESCRIPTION
#' @param method PARAM_DESCRIPTION, Default: 'zscore'
#' @param exp.transform PARAM_DESCRIPTION, Default: F
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

CalcScoreGeneral<-function(matr, f, method = "zscore", exp.transform = F){
  library(dplyr)
  library(rlist)
  library(mosaic)
  library(purrr)
  library(tidyr)
  f <- factor(f)
  if(exp.transform) matr <- expm1(matr)
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
  }else if(method=="logFC"){
    scores <-
      levels(f) %>%
      lapply(function(x) split(scores,f==x)) %>%
      lapply(function(x) map2(x$`TRUE`, x$`FALSE`, function(x,y) log(mean(x+1)/mean(y+1)))) %>%
      lapply(list.rbind) %>%
      list.cbind() %>%
      as.data.frame() %>%
      setNames(., levels(f))
  }
  scores <- drop_na(scores)
  return(scores)
}

#' @title FUNCTION_TITLE
#' @description FUNCTION_DESCRIPTION
#' @param Seu PARAM_DESCRIPTION
#' @param features PARAM_DESCRIPTION
#' @param group.by PARAM_DESCRIPTION
#' @param method PARAM_DESCRIPTION
#' @param assay PARAM_DESCRIPTION, Default: 'RNA'
#' @param exp.transform PARAM_DESCRIPTION, Default: F
#' @return OUTPUT_DESCRIPTION
#' @details DETAILS
#' @examples
#' \dontrun{
#' if(interactive()){
#'  #EXAMPLE1
#'  }
#' }
#' @rdname CalcScoreGeneral_v3
#' @export

CalcScoreGeneral_v3<-function(Seu, features, group.by, method, assay = "RNA", exp.transform = F){
  library(Seurat)
  if(any(!features %in% rownames(Seu))) {
    WrongGenes <- setdiff(features, rownames(Seu))
    warning(paste0(paste(WrongGenes, collapse = ", "), " are not detected in the matrix"))
    features <- intersect(features, rownames(Seu))
  }
  matr <- as.matrix(GetAssayData(Seu, assay = assay)[features,])
  f <- factor(Seu@meta.data[, group.by])
  return(CalcScoreGeneral(matr, f, method, exp.transform))
}

# matr <- GetAssayData(seu)[c("Selp","Vwf","Chst4","Trp53i11","Cdh5"),]
# f <- seu@meta.data$cluster
# CalcScoreGeneral(matr, f, "LogFC") %>%
#   Heatmap()
