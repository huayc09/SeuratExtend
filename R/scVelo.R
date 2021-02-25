#' @title FUNCTION_TITLE
#' @description FUNCTION_DESCRIPTION
#' @param seu PARAM_DESCRIPTION
#' @param export.name PARAM_DESCRIPTION
#' @param velocyto.loom.path PARAM_DESCRIPTION
#' @param velocyto.loom.filenames PARAM_DESCRIPTION
#' @param fxn.convert.loomtoseurat.cellname PARAM_DESCRIPTION, Default: NULL
#' @param prefix PARAM_DESCRIPTION, Default: NULL
#' @param postfix PARAM_DESCRIPTION, Default: NULL
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
  function(seu, export.name, velocyto.loom.path, velocyto.loom.filenames,
           fxn.convert.loomtoseurat.cellname = NULL, prefix = NULL, postfix = NULL){
    library(Seurat)
    library(rlist)
    library(dplyr)
    library(loomR)

    check.metadata.na <- sapply(seu@meta.data, function(x) any(is.na(x)))
    if(any(check.metadata.na)) {
      message(paste0("NA(s) exist in meta data: ",
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

    genes <- rownames(seu)
    cells <- colnames(seu)
    seu.subset.genes <- F

    lfiles <- list()
    cellname <- list()
    # ambiguous <- list()
    spliced <- list()
    unspliced <- list()

    for (i in 1:length(velocyto.loom.filenames)) {
      lfiles[[i]] <- connect(filename = file.path(velocyto.loom.path, velocyto.loom.filenames[i]),
                             mode = "r+", skip.validate = T)

      if(is.null(fxn.convert.loomtoseurat.cellname)){
        cellname[[i]] <-
          strsplit(lfiles[[i]][["col_attrs"]][["CellID"]][], split = ":", fixed = T) %>%
          sapply(function(x) x[2])
        cellname[[i]] <- paste0(prefix[i],sub("x", "", cellname[[i]]),postfix[i])
      }else{
        cellname[[i]] <- fxn.convert.loomtoseurat.cellname(velocyto.loom.filenames = i,
                                                           loom.cellname = lfiles[[i]][["col_attrs"]][["CellID"]][])
      }
      if(!any(cellname[[i]] %in% cells)) stop(paste0(velocyto.loom.filenames[i],": wrong prefix/postfix?"))

      g <- lfiles[[i]][["row_attrs"]][["Gene"]][] %>% make.unique()
      n <- length(setdiff(genes, g))
      if(n>0){
        message("There are ", n, " genes in the Seurat object which cannot be found in the Velocyto output. ",
                "Did you use different version of reference genome?")
        genes <- intersect(genes, g)
        seu.subset.genes <- T
      }
      # ambiguous[[i]] <-
      #   lfiles[[i]][["layers"]][["ambiguous"]][,] %>%
      #   `colnames<-`(g) %>%
      #   `rownames<-`(cellname[[i]]) %>%
      #   .[rownames(.) %in% cells, genes]
      spliced[[i]] <-
        lfiles[[i]][["layers"]][["spliced"]][,] %>%
        `colnames<-`(g) %>%
        `rownames<-`(cellname[[i]]) %>%
        .[!is.na(rownames(.)),] %>%
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

    if(seu.subset.genes){
      seu.new <- CreateSeuratObject(GetAssayData(seu, slot = "counts")[genes,], meta.data = seu@meta.data)
      seu.new <- NormalizeData(seu.new)
      seu.new <- FindVariableFeatures(seu.new)
      seu.new[["RNA"]]@scale.data <- seu[[DefaultAssay(seu)]]@scale.data[genes,]
      for (i in names(seu@reductions)) {
        seu.new[[i]] <- CreateDimReducObject(embeddings = Embeddings(seu, reduction = i))
      }
      seu <- seu.new
    }
    pfile <- as.loom(x = seu, filename = export.name, overwrite = T)

    pfile$add.layer(layers = list(#"ambiguous" = ambiguous,
      "spliced" = spliced,
      "unspliced" = unspliced))

    pfile$close_all()
  }


scVelo.RunBasic <- function(loom, save.adata = "adata.obj"){
  library(reticulate)
  code <- paste0("import scvelo as scv\n",
                 "import scanpy as sc\n",
                 "import numpy as np\n",
                 "import matplotlib.cm as cm\n",
                 "scv.logging.print_version()\n",
                 "scv.settings.verbosity = 3\n",
                 "scv.settings.presenter_view = True\n",
                 "scv.set_figure_params('scvelo')\n",
                 "import loompy\n",
                 'ds = loompy.connect("', loom, '")\n',
                 'del ds.layers["norm_data"]\n',
                 'del ds.layers["scale_data"]\n',
                 'ds.close()\n',
                 'adata = scv.read("', loom, '")\n',
                 "adata\n",
                 "scv.pp.filter_and_normalize(adata, min_shared_counts=20, n_top_genes=2000)\n",
                 "scv.pp.moments(adata, n_pcs=30, n_neighbors=30)\n",
                 "scv.tl.velocity(adata)\n",
                 "scv.tl.velocity_graph(adata)\n",
                 "import pickle\n",
                 'pickle_out = open("', save.adata, '","wb")\n',
                 "pickle.dump(adata, pickle_out)\n",
                 "pickle_out.close()\n")
  py_run_string(code)
}
