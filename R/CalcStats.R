#' @include generics.R
#'
NULL

#' @param seu (Seurat version) Seurat object
#' @param features (Seurat version) Features to plot (gene expression, metrics, PC scores,
#' anything that can be retrieved by FetchData), Default: NULL (All features
#' in matrix)
#' @param group.by (Seurat version) A variable name in meta.data to
#' group by, or string with the same length of cells
#' @param cells (Seurat version) Cell names to use, Default: all cells
#' @param slot (Seurat version) Slot to pull feature data for
#' @param assay (Seurat version) Name of assay to use, defaults to the active assay
#' @param priority (Seurat version) If set to "expr", force fetch data from expression matrix
#' instead of meta.data
#' @rdname CalcStats
#' @export

CalcStats.Seurat <-
  function(
    seu,
    features,
    group.by = NULL,
    cells = NULL,
    slot = "data",
    assay = NULL,
    method = "zscore",
    exp.transform = F,
    order = NULL,
    n = Inf,
    p.threshold = 0.05,
    priority = c("expr","none")
  ) {
    Std.matr <- Seu2Matr(
      seu = seu,
      features = features,
      group.by = group.by,
      cells = cells,
      slot = slot,
      assay = assay,
      priority = priority
    )

    scores <- CalcStats.default(
      matr = Std.matr$matr,
      f = Std.matr$f,
      method = method,
      exp.transform = exp.transform,
      t = T,
      order = order,
      n = n,
      p.threshold = p.threshold)
    return(scores)
}

#' @param matr Matrix or data frame.Row - features; columns - cells
#' @param f Factor or vector. Identity of each cell. Should be the
#' same length of cells
#' @param method Should be either "mean", "median", "zscore", "tscore",
#' "p", or "logFC". Default: 'zscore'
#' @param exp.transform Whether to transform the data with
#' \code{\link[base:log]{expm1}}, Default: F
#' @param t If your matrix has features in columns and cells in rows,
#' you should transpose the matrix first. Default: F
#' @param order Re-order rows by "value" or "p" (in t.test)
#' @param n Top n rows of each cluster. Ignored when
#' \code{order} is \code{NULL}. Default: \code{Inf}
#' @param p.threshold When \code{order = "p"}, rows with pvalue greater
#' than \code{p.threshold} are removed.
#' @rdname CalcStats
#' @export

CalcStats.default <- function(
  matr,
  f,
  method = "zscore",
  exp.transform = F,
  t = F,
  order = NULL,
  n = Inf,
  p.threshold = 0.05
) {
  library(dplyr)
  library(mosaic)

  f <- factor(f)
  if(exp.transform) matr <- expm1(matr)
  if(!t) matr <- t(matr)
  matr <- as.data.frame(matr)
  method <- tolower(method)
  if(nrow(matr) != length(f)) {
    stop("'f' should be the same length as number of matrix columns. \n",
         "Maybe need transpose? (t = TRUE)")
  }

  if(method == "logfc" & !exp.transform) {
    warning("When calculating Log2FC of log-normalized gene expression data, ",
            "you should consider set 'exp.transform' to TRUE")
  }
  scores <- switch(
    method,
    mean = CalcMean(matr, f, "mean"),
    median = CalcMean(matr, f, "median"),
    zscore = CalcZscore(matr, f),
    tscore = CalcTscore(matr, f, "tscore"),
    p = CalcTscore(matr, f, "p"),
    logfc = CalcTscore(matr, f, "logfc")
  )
  scores <- tidyr::drop_na(scores)

  if(!is.null(order)) {
    row.max <- apply(scores, 1, which.max)
    scores.order <- split(scores, row.max)

    if(order == "value" | order == method) {
      for (i in names(scores.order)) {
        tmp <- scores.order[[i]]
        tmp <- tmp[order(tmp[[as.numeric(i)]], decreasing = T),]
        tmp <- head(tmp, n)
        scores.order[[i]] <- tmp
      }
      scores.order <- setNames(scores.order, NULL)
      scores <- rlist::list.rbind(scores.order)
    } else if (order == "p") {
      for (i in names(scores.order)) {
        tmp <- scores.order[[i]]
        scores2 <- matr[,rownames(tmp),drop = F]
        clust.id <- levels(f)[as.numeric(i)]

        scores2 <- CalcTscore(scores2, f, "p", idents = clust.id)
        scores2 <- scores2[scores2[[1]] > -log10(p.threshold), ,drop = F]
        feature.p <- rownames(scores2)[order(scores2[[1]], decreasing = T)]
        tmp <- tmp[head(feature.p, n),]
        scores.order[[i]] <- tmp
      }
      scores.order <- setNames(scores.order, NULL)
      scores <- rlist::list.rbind(scores.order)
    } else {
      stop("'order' must be 'value' or 'p'")
    }
  }
  return(scores)
}

#' @title CalcScoreGeneral, ScoreAndOrder
#' @description Alias of \code{\link[SeuratExtend:CalcStats]{CalcStats()}}
#' @seealso \code{\link[SeuratExtend:CalcStats]{CalcStats()}}
#' @rdname CalcScoreGeneral
#' @export

CalcScoreGeneral <- CalcStats.default

#' @rdname CalcScoreGeneral
#' @export

CalcScoreGeneral_v3 <- CalcStats.Seurat

#' @rdname CalcScoreGeneral
#' @export

ScoreAndOrder <- CalcStats.default

# Internal ----------------------------------------------------------------

CalcMean <- function(matr, f, method = c("mean","median")) {
  matr <- split(matr, f)
  method <- method[1]
  fun <- switch(
    method,
    mean = base::mean,
    median = base::median
    )
  matr <- lapply(matr, function(x) apply(x,2,fun))
  matr <- rlist::list.cbind(matr)
  matr <- as.data.frame(matr)
  return(matr)
}

CalcZscore <- function(matr, f) {
  matr <- CalcMean(matr, f)
  matr <- apply(matr, 1, mosaic::zscore, simplify = F)
  matr <- rlist::list.rbind(matr)
  matr <- as.data.frame(matr)
  return(matr)
}

CalcTscore <-
  function(
    matr, f,
    method = c("tscore","p","logfc"),
    idents = NULL
  ) {
    method <- method[1]
    idents <- idents %||% levels(f)
    t.list <- list()
    for (i in idents) {
      ident.info <- CheckIdent(f, ident.1 = i)
      t.list[[i]] <- apply(matr, 2, function(x) {
        x1 <- x[ident.info$cell.1]
        x2 <- x[ident.info$cell.2]
        switch(
          method,
          tscore = t.test(x1, x2)[["statistic"]],
          p = -log10(t.test(x1, x2)[["p.value"]]),
          logfc = log2( (mean(x1) + 1) / (mean(x2) + 1) )
        )
      })
    }
    matr <- rlist::list.cbind(t.list)
    matr <- as.data.frame(matr)
    return(matr)
  }
