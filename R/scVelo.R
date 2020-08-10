#' @title FUNCTION_TITLE
#' @description FUNCTION_DESCRIPTION
#' @param seu PARAM_DESCRIPTION
#' @param export.name PARAM_DESCRIPTION
#' @param prefix PARAM_DESCRIPTION
#' @param postfix PARAM_DESCRIPTION
#' @param velocyto.loom.path PARAM_DESCRIPTION
#' @param velocyto.loom.filenames PARAM_DESCRIPTION
#' @return OUTPUT_DESCRIPTION
#' @details DETAILS
#' @examples
#' \dontrun{
#' if(interactive()){
#'  #EXAMPLE1
#'  }
#' }
#' @rdname scVelo.SeuratToLoom
#' @export

scVelo.SeuratToLoom <-
  function(seu, export.name, prefix, postfix, velocyto.loom.path, velocyto.loom.filenames){
    library(Seurat)
    library(rlist)
    library(dplyr)
    library(loomR)

    check.metadata.na <- sapply(seu@meta.data, function(x) any(is.na(x)))
    if(any(check.metadata.na)) {
      warning(paste0("NA(s) exist in meta data: ",
                     paste0(check.metadata.na %>% names(.)[.], collapse = ", "),
                     "\nColumn(s) removed"))
      seu@meta.data[,check.metadata.na] <- NULL
    }

    check.loom.files <- velocyto.loom.filenames %in% list.files(velocyto.loom.path)
    if(!any(check.loom.files)) {
      stop(paste0("Loom file(s) not found in path ",
                  velocyto.loom.path, ": ",
                  paste0(velocyto.loom.filenames[!check.loom.files], collapse = ", ")))
    }

    pfile <- as.loom(x = seu, filename = export.name, overwrite = T)
    genes <- pfile[["row_attrs"]][["Gene"]][]
    cells <- pfile[["col_attrs"]][["CellID"]][]

    lfiles <- list()
    cellname <- list()
    # ambiguous <- list()
    spliced <- list()
    unspliced <- list()

    for (i in 1:length(velocyto.loom.filenames)) {
      lfiles[[i]] <- connect(filename = paste0(velocyto.loom.path, velocyto.loom.filenames[i]),
                             mode = "r+", skip.validate = T)
      cellname[[i]] <-
        strsplit(lfiles[[i]][["col_attrs"]][["CellID"]][], split = ":", fixed = T) %>%
        sapply(function(x) x[2])
      cellname[[i]] <- paste0(prefix[i],sub("x", "", cellname[[i]]),postfix[i])
      if(!any(cellname[[i]] %in% cells)) stop(paste0(velocyto.loom.filenames[i],": wrong prefix/postfix?"))
      g <- lfiles[[i]][["row_attrs"]][["Gene"]][] %>% make.unique()
      # ambiguous[[i]] <-
      #   lfiles[[i]][["layers"]][["ambiguous"]][,] %>%
      #   `colnames<-`(g) %>%
      #   `rownames<-`(cellname[[i]]) %>%
      #   .[rownames(.) %in% cells, genes]
      spliced[[i]] <-
        lfiles[[i]][["layers"]][["spliced"]][,] %>%
        `colnames<-`(g) %>%
        `rownames<-`(cellname[[i]]) %>%
        .[rownames(.) %in% cells, genes]
      unspliced[[i]] <-
        lfiles[[i]][["layers"]][["unspliced"]][,] %>%
        `colnames<-`(g) %>%
        `rownames<-`(cellname[[i]]) %>%
        .[rownames(.) %in% cells, genes]
    }

    # ambiguous <- ambiguous %>%
    #   list.rbind() %>%
    #   .[cells,]
    spliced <- spliced %>%
      list.rbind() %>%
      .[cells,]
    unspliced <- unspliced %>%
      list.rbind() %>%
      .[cells,]

    pfile$add.layer(layers = list(#"ambiguous" = ambiguous,
                                  "spliced" = spliced,
                                  "unspliced" = unspliced))

    pfile$close_all()
  }
