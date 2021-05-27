#' @include generics.R
#'
NULL

#' @param seu (Seurat version) Seurat object
#' @param features (Seurat version) Features to plot (gene expression, metrics, PC scores,
#' anything that can be retreived by FetchData), Default: NULL (All features
#' in matrix)
#' @param group.by (Seurat version) A variable name in meta.data to
#' group the violin plots by
#' @param cells (Seurat version) Cell names to use, Default: all cells
#' @param slot Slot to pull feature data for
#' @rdname CalcStats
#' @export

CalcStats.Seurat <-
  function(
    seu,
    features,
    group.by = NULL,
    cells = NULL,
    slot = "data",
    method = "zscore",
    exp.transform = F
  ) {
    library(Seurat)
    cells <- cells %||% colnames(seu)
    matr <- FetchData(object = seu, vars = features, cells = cells, slot = slot)
    if(is.null(group.by)) {
      f <- factor(Idents(seu)[cells])
    }else{
      f <- factor(seu[[group.by]][cells,])
    }
    scores <- CalcStats.default(
      matr,
      f = f,
      method = method,
      exp.transform = exp.transform,
      t = T)
    return(scores)
}

#' @param matr Matrix or data frame.Row - features; columns - cells
#' @param f Factor or vector. Identity of each cell. Should be the
#' same length of cells
#' @param method Should be either "mean", "median", "zscore", "tscore",
#' "p", or "logFC". Default: 'zscore'
#' @param exp.transform Whether to transform the data with
#' \code{\link[base:expm1]{expm1}}, Default: F
#' @param t If your matrix has features in columns and cells in rows,
#' you should transpose the matrix first. Default: F
#' @rdname CalcStats
#' @export

CalcStats.default <- function(
  matr,
  f,
  method = "zscore",
  exp.transform = F,
  t = F
) {
  library(dplyr)
  library(rlist)
  library(mosaic)
  library(purrr)
  library(tidyr)

  f <- factor(f)
  if(exp.transform) matr <- expm1(matr)
  if(!t) scores <- t(matr) else scores <- matr
  scores <- as.data.frame(scores)
  method <- tolower(method)
  if(nrow(scores) != length(f)) {
    stop("'f' should be the same length as number of matrix columns. \n",
         "Maybe need transpose? (t = TRUE)")
  }

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
    scores <-
      CalcStats.default(scores, f, "mean", t = T) %>%
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
  }else if(method=="logfc"){
    scores <-
      levels(f) %>%
      lapply(function(x) split(scores,f==x)) %>%
      lapply(function(x) map2(x$`TRUE`, x$`FALSE`, function(x,y) log(mean(x+1)/mean(y+1)))) %>%
      lapply(list.rbind) %>%
      list.cbind() %>%
      as.data.frame() %>%
      setNames(., levels(f))
  }else{
    stop("'method' should be one of 'mean', 'median', 'tscore', 'zscore', 'p' or 'logFC")
  }
  scores <- drop_na(scores)
  return(scores)
}

#' @title CalcScoreGeneral
#' @description Alias of \code{\link[SeuratExtend:CalcStats]{CalcStats()}}
#' @seealso \code{\link[SeuratExtend:CalcStats]{CalcStats()}}
#' @rdname CalcScoreGeneral
#' @export

CalcScoreGeneral <- CalcStats.default

#' @rdname CalcScoreGeneral
#' @export

CalcScoreGeneral_v3 <- CalcStats.Seurat
