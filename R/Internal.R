# convert Seurat object to standard expression matrix

Seu2Matr <-
  function(
    seu,
    features,
    group.by = NULL,
    split.by = NULL,
    cells = NULL,
    slot = "data",
    assay = NULL,
    priority = c("expr","none"),
    verbose = TRUE
  ) {
    if(!require(SeuratObject)) library(Seurat)

    if(!is.null(assay)) DefaultAssay(seu) <- assay

    cells <- cells %||% colnames(seu)
    if(is.logical(cells)) {
      if(length(cells) != ncol(seu)) {
        stop("Logical value of 'cells' should be the same length as cells in ",
             "Seurat object")
      } else {
        cells.l <- cells
        cells <- colnames(seu)[cells]
      }
    } else {
      if(all(!cells %in% colnames(seu))) {
        stop("'cells' not found in Seurat object")
      }else if(any(!cells %in% colnames(seu))) {
        cells.out <- setdiff(cells, colnames(seu))
        stop(length(cells.out), " cell(s) not found in Seurat object: '",
             cells.out[1], "'...")
      }
      cells.l <- colnames(seu) %in% cells
    }

    if(priority[1] == "expr") {
      check.rep <- intersect(rownames(seu), colnames(seu@meta.data))
      if(any(features %in% check.rep)) {
        feature.rep <- intersect(features, check.rep)
        if(verbose) {
          warning("'Features' found in both expression matrix and 'meta.data' (",
                  feature.rep[1],
                  "...). Ignore the 'meta.data'. If you need to fetch data from ",
                  "'meta.data' please set 'priority' to 'none'")
        }
        seu@meta.data[,feature.rep] <- NULL
      }
    }

    if (utils::packageVersion("SeuratObject") >= "5.0.0") {
      matr <- FetchData(object = seu, vars = features, cells = cells, layer = slot, clean = "none")
    } else {
      matr <- FetchData(object = seu, vars = features, cells = cells, slot = slot)
    }

    if(is.null(group.by)) {
      f <- factor(Idents(seu)[cells])
    }else if(length(group.by) == 1) {
      if(group.by %in% colnames(seu@meta.data)){
        f <- factor(seu[[group.by]][cells,])
      }else{
        stop("Cannot find '", group.by, "' in meta.data")
      }
    }else if(length(group.by) == ncol(seu)) {
      f <- factor(group.by[cells.l])
    }else if(length(group.by) == length(cells)) {
      f <- factor(group.by)
    }else{
      stop("'group.by' should be variable name in 'meta.data' or ",
           "string with the same length of cells")
    }
    names(f) <- cells

    if(is.null(split.by)) {
      f2 <- NULL
    }else if(length(split.by) == 1) {
      if(split.by %in% colnames(seu@meta.data)){
        f2 <- factor(seu[[split.by]][cells,])
        names(f2) <- cells
      }else{
        stop("Cannot find '", split.by, "' in meta.data")
      }
    }else if(length(split.by) == ncol(seu)) {
      f2 <- factor(split.by[cells.l])
      names(f2) <- cells
    }else if(length(split.by) == length(cells)) {
      f2 <- factor(split.by)
      names(f2) <- cells
    }else{
      stop("'split.by' should be variable name in 'meta.data' or ",
           "string with the same length of cells")
    }
    return(list(matr = matr, f = f, f2 = f2))
  }

# Check if ident.1 and ident.2 are in the factor levels and return
# logical vectors
CheckIdent <-
  function(
    f,
    ident.1 = NULL,
    ident.2 = NULL
  ) {
    lv <- levels(f)
    nlv <- 1:nlevels(f)

    ident.1 <- ident.1 %||% lv[1]
    if(ident.1 %in% lv) {
      cell.1 <- (f == ident.1)
      name.1 <- ident.1
    } else if(ident.1 %in% nlv) {
      name.1 <- lv[ident.1]
      cell.1 <- (f == name.1)
    } else stop("'ident.1' not found. Possible values: \n",
                paste(lv, collapse = ", "))

    if(is.null(ident.2)) {
      name.2 <- paste0("non-", name.1)
      cell.2 <- (f != name.1)
    } else if(ident.2 %in% lv) {
      cell.2 <- (f == ident.2)
      name.2 <- ident.2
    } else if(ident.2 %in% nlv) {
      name.2 <- lv[ident.2]
      cell.2 <- (f == name.2)
    } else stop("'ident.2' not found. Possible values: \n",
                paste(lv, collapse = ", "))

    if(identical(name.1, name.2))
      warning("'ident.1' and 'ident.2' are the same??")

    return(list(
      cell.1 = cell.1,
      name.1 = name.1,
      cell.2 = cell.2,
      name.2 = name.2
    ))
  }


# imported from rlang

`%||%` <- function (x, y)
{
  if (is.null(x))
    y
  else x
}

is_empty <- function (x) length(x) == 0

# Check species

check_spe <- function(spe){
  if(is.null(spe)) stop("species undefined: options(spe = c(\"mouse\", \"human\"))")
}

# Check package version
check_pkg_version <- function(package, version.min) {
  version.min <- package_version(version.min)
  version.pkg <- packageVersion(package)
  if(version.pkg < version.min) {
    stop("Current version of package '", package,"' is ", version.pkg,
         ", but version ", version.min, " is required")
  }
}

# message only once
message_only_once <- function(title, message) {
  option.title <- paste0("seuratextend_warning_",title)
  if(getOption(option.title, TRUE)) {
    message(message, "\nThis message is shown once per session")
    options(setNames(list(FALSE), option.title))
  }
}
