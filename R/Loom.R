#' @rdname adata.Load
#' @export

Seu2Loom <- function(
    seu,
    filename,
    add.normdata = FALSE,
    add.metadata = TRUE,
    layers = NULL,
    overwrite = FALSE
) {
  library(hdf5r)
  library(Seurat)
  library(tools)

  # Check file extension and modify filename if needed
  if (!grepl(pattern = "^loom$", x = file_ext(x = filename))) {
    filename <- paste0(filename, ".loom")
  }

  # Check if file exists and handle according to 'overwrite' parameter
  if (file.exists(filename)) {
    if (isTRUE(x = overwrite)) {
      warning("Overwriting previous file ", filename,
              call. = FALSE, immediate. = TRUE)
      file.remove(filename)
    }
    else {
      stop("Loom file at ", filename, " already exists",
           call. = FALSE)
    }
  }

  # Check dimensions of layers
  if (!is.null(layers)) {
    seu_dim <- dim(seu)
    for (layer_name in names(layers)) {
      layer_matrix <- layers[[layer_name]]
      if (!all(seu_dim == rev(dim(layer_matrix)))) {
        stop(paste("The dimensions of the", layer_name, "layer do not match the transposed dimensions of the Seurat object. The Seurat object has", seu_dim[2], "cells and", seu_dim[1], "features. Please adjust the", layer_name, "layer matrix to have", seu_dim[2], "rows (cells) and", seu_dim[1], "columns (features) and try again."))
      }
    }
  }

  # Create a string datatype object with UTF-8 encoding
  utf8_string_type <- H5T_STRING$new(size = Inf)$set_cset(cset = h5const$H5T_CSET_UTF8)

  # Open a new HDF5 file
  file <- H5File$new(filename, mode = "w")

  # Create and set attributes
  file$create_group("attrs")
  file[["attrs"]]$create_dataset("LOOM_SPEC_VERSION", "3.0.0", dtype = utf8_string_type)

  # Create necessary groups
  row_attrs <- file$create_group("row_attrs")
  col_attrs <- file$create_group("col_attrs")
  file$create_group("col_graphs")
  file$create_group("row_graphs")
  layers_group <- file$create_group("layers")

  # Extract count matrix and transpose to have cells as rows and genes as columns
  count_matrix <- t(as.matrix(.get_assay_data_compat(seu, slot = "counts")))

  # If add.normdata is TRUE, extract normalized data and set it as the main matrix
  if (isTRUE(add.normdata)) {
    norm_matrix <- t(as.matrix(.get_assay_data_compat(seu, slot = "data")))
    file$create_dataset("matrix", norm_matrix)
    layers_group$create_dataset("counts", count_matrix)
  } else {
    # Create the main matrix dataset with count matrix
    file$create_dataset("matrix", count_matrix)
  }

  # Add cell IDs to col_attrs group
  cells <- colnames(seu)
  col_attrs$create_dataset("CellID", cells, dtype = utf8_string_type)

  # Add gene names to row_attrs group
  genes <- rownames(seu)
  row_attrs$create_dataset("Gene", genes, dtype = utf8_string_type)

  # Add cell metadata to col_attrs group
  if(isTRUE(add.metadata)) {
    meta_data <- seu@meta.data
    for (colname in colnames(meta_data)) {
      # Convert factor columns to character
      if (is.factor(meta_data[[colname]])) {
        meta_data[[colname]] <- as.character(meta_data[[colname]])
      }
      # Add dataset to col_attrs group
      if (is.character(meta_data[[colname]])) {
        col_attrs$create_dataset(colname, meta_data[[colname]], dtype = utf8_string_type)
      } else {
        col_attrs$create_dataset(colname, meta_data[[colname]])
      }
    }
  }

  # Add layers to the "layers" group
  if (!is.null(layers)) {
    for (layer_name in names(layers)) {
      layers_group$create_dataset(layer_name, layers[[layer_name]])
    }
  }

  # Close the file
  file$close_all()
}
