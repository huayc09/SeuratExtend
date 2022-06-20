#' @title FUNCTION_TITLE
#' @description FUNCTION_DESCRIPTION
#' @param pairs PARAM_DESCRIPTION
#' @param sep PARAM_DESCRIPTION, Default: '_'
#' @return OUTPUT_DESCRIPTION
#' @details DETAILS
#' @examples
#' \dontrun{
#' if(interactive()){
#'  #EXAMPLE1
#'  }
#' }
#' @rdname pair_rev
#' @export

pair_rev <- function(pairs, sep = "_") {
  library(dplyr)
  pairs_rev <-
    strsplit(pairs, split = sep, fixed = T) %>%
    sapply(function(x) paste(rev(x), collapse = sep))
  return(pairs_rev)
}

#' @title FUNCTION_TITLE
#' @description FUNCTION_DESCRIPTION
#' @param seu PARAM_DESCRIPTION
#' @param feature PARAM_DESCRIPTION
#' @param ident PARAM_DESCRIPTION
#' @param group.by PARAM_DESCRIPTION, Default: NULL
#' @param DefaultAssay PARAM_DESCRIPTION, Default: 'RNA'
#' @return OUTPUT_DESCRIPTION
#' @details DETAILS
#' @examples
#' \dontrun{
#' if(interactive()){
#'  #EXAMPLE1
#'  }
#' }
#' @rdname gene_percent
#' @export

gene_percent <- function(seu, feature, ident, group.by = NULL, DefaultAssay = "RNA") {
  library(Seurat)
  DefaultAssay(seu) <- DefaultAssay
  if(is.null(group.by)) {
    cells <- names(Idents(seu))[Idents(seu) == ident]
  }else{
    cells <- colnames(seu)[seu[[group.by]] == ident]
  }
  m <- FetchData(seu, vars = feature, cells = cells)
  return(apply(m, 2, function(x) sum(x>0)) / nrow(m))
}

#' @title FUNCTION_TITLE
#' @description FUNCTION_DESCRIPTION
#' @param ident PARAM_DESCRIPTION
#' @param seu PARAM_DESCRIPTION
#' @param feature PARAM_DESCRIPTION, Default: rownames(seu)
#' @param group.by PARAM_DESCRIPTION, Default: NULL
#' @param pct PARAM_DESCRIPTION, Default: 0.1
#' @param DefaultAssay PARAM_DESCRIPTION, Default: 'RNA'
#' @return OUTPUT_DESCRIPTION
#' @details DETAILS
#' @examples
#' \dontrun{
#' if(interactive()){
#'  #EXAMPLE1
#'  }
#' }
#' @rdname gene_expressed
#' @export

gene_expressed <- function(ident, seu, feature = rownames(seu), group.by = NULL, pct = 0.1,
                           DefaultAssay = "RNA") {
  gene_percent <- gene_percent(seu = seu, feature = feature, ident = ident, group.by = group.by,
                               DefaultAssay = DefaultAssay)
  return(names(gene_percent)[gene_percent >= pct])
}

#' @title FUNCTION_TITLE
#' @description FUNCTION_DESCRIPTION
#' @param fun PARAM_DESCRIPTION
#' @param next.line PARAM_DESCRIPTION, Default: T
#' @param ... PARAM_DESCRIPTION
#' @return OUTPUT_DESCRIPTION
#' @details DETAILS
#' @examples
#' \dontrun{
#' if(interactive()){
#'  #EXAMPLE1
#'  }
#' }
#' @rdname write.py.fun
#' @export

write.py.fun <- function(fun, next.line = T, ...) {
  par <- c(...)
  par <- paste(names(par), par, sep = " = ", collapse = ", ")
  text <- paste0(fun, "(", par, ")")
  if(next.line) text <- paste0(text, "\n")
  return(text)
}

#' @title FUNCTION_TITLE
#' @description FUNCTION_DESCRIPTION
#' @param df PARAM_DESCRIPTION
#' @param n PARAM_DESCRIPTION, Default: 3
#' @return OUTPUT_DESCRIPTION
#' @details DETAILS
#' @examples
#' \dontrun{
#' if(interactive()){
#'  #EXAMPLE1
#'  }
#' }
#' @rdname viewdf
#' @export

viewdf <- function(df, n = 3){
  df_orig <- df
  cl <- class(df)
  dm <- dim(df)
  if(ncol(df) > n * 2){
    ns <- 1:ncol(df)
    df <- cbind(df[,head(ns, n),drop = F],
                matrix(".", ncol = 3, nrow = nrow(df)) %>% `colnames<-`(rep(".",3)),
                df[,tail(ns, n),drop = F])
  }
  if(nrow(df) > n * 2){
    ns <- 1:nrow(df)
    df <- rbind(df[head(ns, n),,drop = F],
                matrix(".", nrow = 3, ncol = ncol(df)) %>% `colnames<-`(colnames(df)) %>% `rownames<-`(c(".",". ",".  ")),
                df[tail(ns, n),,drop = F])
  }
  print(list(df = df,
             class = cl,
             dim = dm),
        quote = F)
}

#' @title FUNCTION_TITLE
#' @description FUNCTION_DESCRIPTION
#' @param ... PARAM_DESCRIPTION
#' @return OUTPUT_DESCRIPTION
#' @details DETAILS
#' @examples
#' \dontrun{
#' if(interactive()){
#'  #EXAMPLE1
#'  }
#' }
#' @rdname assign.list
#' @export

assign.list <- function(...) {
  x <- list(...)
  for (i in names(x)) {
    assign(i, x[[i]],envir = globalenv())
  }
}

#' @title FUNCTION_TITLE
#' @description FUNCTION_DESCRIPTION
#' @param seu PARAM_DESCRIPTION
#' @param dr PARAM_DESCRIPTION
#' @param file PARAM_DESCRIPTION, Default: 'tmp/dr.csv'
#' @return OUTPUT_DESCRIPTION
#' @details DETAILS
#' @examples
#' \dontrun{
#' if(interactive()){
#'  #EXAMPLE1
#'  }
#' }
#' @rdname export.dr
#' @export

export.dr <- function(seu, dr, file = "tmp/dr.csv") {
  if(grepl("^tmp/", file) & !file.exists("tmp")) {
    dir.create("tmp")
  }
  write.csv(Embeddings(seu, reduction = dr), file = file, quote = F)
}

#' @title Proportion of features
#' @description Check the Proportion of positive cells (default: expression above 0)
#' in certain clusters
#' @param seu Seurat object
#' @param feature Features to plot (gene expression, metrics, PC scores, anything that can be
#' retreived by FetchData)
#' @param ident cluster name, Default: all clusters
#' @param group.by A variable name in meta.data to group by, if you don't want to use
#' default Idents of Seurat object
#' @param DefaultAssay Assay to use, Default: 'RNA'
#' @param above Above which value will be considered as positive cell, Default: 0
#' @param total If to calculate proportion in total cells of selected clusters, Default: F
#' @param if.expressed If only return logical value, Default: F
#' @param min.pct Minimum fraction of min.pct cells in the cluster will be considered "expressed",
#' only use this parameter when "if.expressed" is set to TRUE. Default: 0.1
#' @return A matrix
#' @details None.
#' @examples
#' genes <- VariableFeatures(pbmc)[1:10]
#'
#' feature_percent(pbmc, feature = genes)
#'
#' # Count one cell as positive only the expression is above 2
#' feature_percent(pbmc, feature = genes, above = 2)
#'
#' # Only check subset of clusters
#' feature_percent(pbmc, feature = genes, ident = c("B","CD8 T"))
#'
#' # Group by a different variable
#' feature_percent(pbmc, feature = genes, group.by = "orig.ident")
#'
#' # Also check the proportion of expressed cells in total clusters
#' feature_percent(pbmc, feature = genes, total = T)
#'
#' # only show logical value; 'TRUE' if more than 20% cells are positive
#' feature_percent(pbmc, feature = genes, if.expressed = T, min.pct = 0.2)
#'
#' @rdname feature_percent
#' @export

feature_percent <- function(
  seu,
  feature,
  ident = NULL,
  group.by = NULL,
  DefaultAssay = "RNA",
  above = 0,
  total = F,
  if.expressed = F,
  min.pct = 0.1
) {
  library(Seurat)
  library(rlang)
  library(rlist)
  DefaultAssay(seu) <- DefaultAssay
  if(is.null(group.by)){
    f <- Idents(seu)
  } else {
    f <- factor(seu@meta.data[,group.by])
  }
  ident <- ident %||% levels(f)
  cells <- colnames(seu)[f %in% ident]
  m <- FetchData(seu, vars = feature, cells = cells)
  f2 <- f[f %in% ident]
  f2 <- factor(f2)
  m2 <- split(m, f2)
  if(total) m2[["total"]] <- m
  m2 <- lapply(m2, function(x){
    apply(x, 2, function(y) sum(y > above)) / nrow(x)
  })
  m2 <- list.cbind(m2)
  if(if.expressed) m2 <- m2 > min.pct
  return(m2)
}



