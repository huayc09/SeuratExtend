#' @include generics.R
#'
NULL

#' @param seu A Seurat object. Only applicable when using the Seurat method.
#' @param features Features for computation, including gene expression, metrics, PC scores, or any other data retrievable via `FetchData()`. Defaults to NULL, implying all features in the matrix. Only applicable for the Seurat method.
#' @param group.by A variable from `meta.data` for grouping or a character vector of equal length as the number of cells. Only applicable for the Seurat method.
#' @param cells Cell identifiers to be used. Defaults to all cells. Only applicable for the Seurat method.
#' @param slot Slot from which to fetch feature data. Only applicable for the Seurat method.
#' @param assay Name of the assay to use. Defaults to the active assay. Only applicable for the Seurat method.
#' @param priority If set to "expr", fetches data from the expression matrix over `meta.data`. Only applicable for the Seurat method.
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

#' @param matr A matrix or data frame with rows as features and columns as cells.
#' @param f A factor or vector indicating cell identity. It should match the column length of `matr`.
#' @param method Computation method, one of: "mean", "median", "zscore", "tscore", "p", or "logFC". Default: "zscore".
#' @param exp.transform Indicates if data should be transformed with `expm1`. Default: FALSE.
#' @param t Transpose matrix if features are columns and cells are rows. Default: FALSE.
#' @param order Sort rows by either "value" or "p" (t.test).
#' @param n Top n rows for each cluster. Ignored if `order` is NULL. Default: Inf.
#' @param p.threshold Rows with p-value greater than this threshold are omitted when `order = "p"`.
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
