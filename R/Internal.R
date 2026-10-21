# convert Seurat object to standard expression matrix

Seu2Matr <-
  function(
    seu,
    features = NULL,
    group.by = NULL,
    split.by = NULL,
    cells = NULL,
    slot = "data",
    assay = NULL,
    priority = c("expr","none"),
    verbose = TRUE,
    external = NULL
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

    ext <- NULL
    ord <- NULL
    if (!is.null(external)) {
      if (anyDuplicated(cells)) {
        stop("'cells' contains duplicated cell names; each cell must be provided once")
      }
      ext <- .external_validate(external, colnames(seu))
      E <- rownames(ext)
      conflicts <- intersect(E, union(rownames(seu), colnames(seu@meta.data)))
      if (length(conflicts) > 0) {
        warning("Variables found in both 'external' and the Seurat object (",
                paste(conflicts, collapse = ", "),
                "). Values from 'external' are used for these variables.")
      }
      ord <- .external_order(features, E)
      features <- ord$order
    } else if (is.null(features)) {
      stop("'features' must be provided")
    }
    fetch_features <- if (is.null(external)) features else ord$native

    if(priority[1] == "expr") {
      check.rep <- intersect(rownames(seu), colnames(seu@meta.data))
      if(any(fetch_features %in% check.rep)) {
        feature.rep <- intersect(fetch_features, check.rep)
        if(verbose) {
          warning("'Features' found in both expression matrix and 'meta.data' (",
                  feature.rep[1],
                  "...). Ignore the 'meta.data'. If you need to fetch data from ",
                  "'meta.data' please set 'priority' to 'none'")
        }
        seu@meta.data[,feature.rep] <- NULL
      }
    }

    if (is.null(external)) {
      if (utils::packageVersion("SeuratObject") >= "5.0.0") {
        matr <- FetchData(object = seu, vars = features, cells = cells, layer = slot, clean = "none")
      } else {
        matr <- FetchData(object = seu, vars = features, cells = cells, slot = slot)
      }
    } else {
      native_df <- NULL
      if (length(ord$native) > 0) {
        if (utils::packageVersion("SeuratObject") >= "5.0.0") {
          native_df <- FetchData(object = seu, vars = ord$native, cells = cells,
                                 layer = slot, clean = "none")
        } else {
          native_df <- FetchData(object = seu, vars = ord$native, cells = cells,
                                 slot = slot)
        }
      }
      if (!is.null(native_df)) {
        native_df <- .external_align_native(native_df, cells)
      }
      ext_df <- as.data.frame(t(ext[E, cells, drop = FALSE]))
      matr <- data.frame(row.names = cells)
      for (nm in ord$order) {
        if (nm %in% E) {
          matr[[nm]] <- ext_df[[nm]]
        } else {
          matr[[nm]] <- native_df[[nm]]
        }
      }
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

# Validate an external feature matrix (rows = features, columns = cells) and
# align it to the cell order of the Seurat object.
#
# Checks, in order: accepted classes -> numeric content -> non-empty -> row
# names present and unique -> cell-axis alignment (4-tier rule). Returns a
# dense numeric matrix whose columns follow `cell_names`. Never intersects
# cells silently: a partial name match is an error.
.external_validate <- function(external, cell_names) {
  if (is.data.frame(external)) {
    nonnum <- names(external)[!vapply(external, is.numeric, logical(1))]
    if (length(nonnum) > 0) {
      stop("'external' data.frame contains non-numeric column(s): ",
           paste(nonnum, collapse = ", "))
    }
    mat <- as.matrix(external)
  } else if (is.matrix(external) || inherits(external, "Matrix")) {
    mat <- as.matrix(external)
  } else {
    stop("'external' must be a matrix, data.frame or sparse Matrix; got class '",
         paste(class(external), collapse = "/"), "'")
  }
  if (!is.numeric(mat)) {
    stop("'external' must be numeric; got a ", class(mat)[1], " matrix")
  }
  if (nrow(mat) == 0 || ncol(mat) == 0) {
    stop("'external' matrix must have at least one row and one column")
  }
  rn <- rownames(mat)
  if (is.null(rn)) {
    stop("'external' must have row names (the row names are the variable names)")
  }
  if (anyDuplicated(rn)) {
    dups <- unique(rn[duplicated(rn)])
    stop("'external' has duplicated row names: ", paste(dups, collapse = ", "))
  }

  cn <- colnames(mat)
  if (!is.null(cn)) {
    if (anyDuplicated(cn)) {
      dups <- unique(cn[duplicated(cn)])
      stop("'external' has duplicated column names: ",
           paste(dups, collapse = ", "))
    }
    if (identical(cn, cell_names)) {
      return(mat)
    }
    if (length(cn) == length(cell_names) && setequal(cn, cell_names)) {
      message("'external' columns were reordered by cell name to match the ",
              "cell order of the Seurat object")
      return(mat[, cell_names, drop = FALSE])
    }
    n_missing <- sum(!cell_names %in% cn)
    n_extra <- sum(!cn %in% cell_names)
    stop("'external' cell names do not match the Seurat object: ",
         n_missing, " Seurat cell(s) missing from 'external', ",
         n_extra, " extra cell(s) in 'external'")
  }
  if (ncol(mat) == length(cell_names)) {
    message("'external' matrix has no cell names; positional alignment to the ",
            "Seurat cell order is assumed. Please make sure the columns are ",
            "consistent.")
    colnames(mat) <- cell_names
    return(mat)
  }
  stop("'external' has no column names and ", ncol(mat),
       " column(s), which does not match the ", length(cell_names),
       " cell(s) in the Seurat object")
}

# Resolve the features x external interaction (F2/F3).
#   native: features not present in the matrix, in features order
#   order:  features order first, then unmentioned matrix rows in matrix order
.external_order <- function(features, E) {
  if (is.null(features)) {
    features <- E
  }
  features <- unique(as.character(features))
  list(
    native = features[!features %in% E],
    order = c(features, E[!E %in% features])
  )
}

# Seatbelt for the external path: FetchData is expected to return rows in the
# order of `cells`; if a future engine changes that, reorder by name or stop
# rather than silently mixing cell values.
.external_align_native <- function(native_df, cells) {
  if (identical(rownames(native_df), cells)) {
    return(native_df)
  }
  rn <- rownames(native_df)
  if (!is.null(rn) && length(rn) == length(cells) && setequal(rn, cells)) {
    return(native_df[cells, , drop = FALSE])
  }
  stop("Fetched data rows cannot be aligned to the requested cells by name")
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

check_spe <- function(spe) {
  if (is.null(spe)) {
    stop("Species not defined. Please set with options(spe = \"your_species\") or provide directly to the function. Default supported species are 'human' and 'mouse'. For custom species, please refer to the documentation on using custom data.")
  }
  return(TRUE)
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

# Helper function to handle GetAssayData compatibility between Seurat v4 and v5
.get_assay_data_compat <- function(object, slot = "data", assay = NULL) {
  # Check if we're using SeuratObject v5 or later
  if (utils::packageVersion("SeuratObject") >= "5.0.0") {
    return(GetAssayData(object, layer = slot, assay = assay))
  } else {
    return(GetAssayData(object, slot = slot, assay = assay))
  }
}

# Internal function to limit values to a specified range
# Similar to Seurat's MinMax function
# Used by DotPlot2, Heatmap, and potentially other visualization functions
# @param data Numeric vector or matrix to be clipped
# @param min Minimum cutoff value. Values below this will be set to min. Default: NULL (no lower limit)
# @param max Maximum cutoff value. Values above this will be set to max. Default: NULL (no upper limit)
# @return Data with values clipped to the specified range
MinMax_internal <- function(data, min = NULL, max = NULL) {
  if (!is.null(min)) {
    data[data < min] <- min
  }
  if (!is.null(max)) {
    data[data > max] <- max
  }
  return(data)
}