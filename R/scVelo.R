#' @title Export loom for scVelo
#' @description Export Seurat object to loom for scVelo run, and merge
#' velocyto output (spliced/unspliced matrix) into the the loom
#' @param seu Seurat object
#' @param export.name Name of loom file
#' @param velocyto.loom.path Directory to velocyto output (loom file), Default: getwd()
#' @param velocyto.loom.filenames File names of velocyto output (loom file)
#' @param fxn.convert.loomtoseurat.cellname (Experimental) function to convert velocyto
#' cell names to Seurat cell names
#' @param prefix Seurat cell id prefix to barcode (usually sample name)
#' @param postfix Seurat cell id postfix to barcode (usually '-1')
#' @return a loom file
#' @details See above
#' @rdname scVelo.SeuratToLoom
#' @export

scVelo.SeuratToLoom <- function(
    seu,
    export.name,
    velocyto.loom.path = NULL,
    velocyto.loom.filenames,
    fxn.convert.loomtoseurat.cellname = NULL,
    cell.id.match.table = NULL,
    prefix = NULL,
    postfix = NULL
){
  library(Seurat)
  library(rlist)
  library(dplyr)
  library(rlang)
  import("loomR")
  import("SeuratDisk")

  check.metadata.na <- sapply(seu@meta.data, function(x) any(is.na(x)))
  if(any(check.metadata.na)) {
    message(paste0("NA(s) exist in meta data: ",
                   paste0(check.metadata.na %>% names(.)[.], collapse = ", "),
                   "\nColumn(s) removed"))
    seu@meta.data[,check.metadata.na] <- NULL
  }

  if(!is.null(velocyto.loom.path)) {
    loom.files <- file.path(velocyto.loom.path, velocyto.loom.filenames)
  } else {
    loom.files <- velocyto.loom.filenames
  }
  if(any(!file.exists(loom.files))) {
    stop(paste0("Loom file(s) not found: ",
                paste0(loom.files[!file.exists(loom.files)], collapse = ", ")
    )
    )
  }

  # orig.assay <- DefaultAssay(seu)
  DefaultAssay(seu) <- "RNA"
  genes <- rownames(seu)
  cells <- colnames(seu)
  seu.subset.genes <- F

  lfiles <- list()
  cellname <- list()
  # ambiguous <- list()
  spliced <- list()
  unspliced <- list()

  for (i in seq_along(velocyto.loom.filenames)) {
    lfiles[[i]] <- connect(filename = loom.files[i],
                           mode = "r+", skip.validate = T)

    if(is.null(fxn.convert.loomtoseurat.cellname)){
      cellname[[i]] <-
        strsplit(lfiles[[i]][["col_attrs"]][["CellID"]][], split = ":", fixed = T) %>%
        sapply(function(x) x[2])
      cellname[[i]] <- paste0(prefix[i],sub("x", "", cellname[[i]]),postfix[i])
    }else{
      cellname[[i]] <- fxn.convert.loomtoseurat.cellname(
        velocyto.loom.filenames = i,
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
  # spliced <- spliced %>%
  #   list.rbind() %>%
  #   .[cells,]
  # unspliced <- unspliced %>%
  #   list.rbind() %>%
  #   .[cells,]
  spliced <- list.rbind(spliced)
  unspliced <- list.rbind(unspliced)
  cell.looms <- rownames(spliced)
  nc <- length(setdiff(cells, cell.looms))
  if(nc>0){
    message("There are ", nc, " cells in the Seurat object which cannot be found in the Velocyto output. ")
    cells <- intersect(cells, cell.looms)
  }
  spliced <- spliced[cells,]
  unspliced <- unspliced[cells,]

  seu <- subset(seu, cells = cells, features = genes)
  seu.new <- CreateSeuratObject(GetAssayData(seu, slot = "counts"), meta.data = seu.new@meta.data)
  # if(add.reduction) {
  #   seu.new <- NormalizeData(seu.new)
  #   seu.new <- FindVariableFeatures(seu.new)
  #   seu.new <- ScaleData(seu.new)
  #   for (i in names(seu@reductions)) {
  #     seu.new[[i]] <- CreateDimReducObject(embeddings = Embeddings(seu, reduction = i))
  #   }
  # }
  pfile <- as.loom(x = seu.new, filename = export.name, overwrite = T)

  pfile$add_layer(x = t(spliced), name = "spliced")
  pfile$add_layer(x = t(unspliced), name = "unspliced")

  pfile$close_all()
}

#' @title FUNCTION_TITLE
#' @description FUNCTION_DESCRIPTION
#' @param loom PARAM_DESCRIPTION
#' @param save.adata PARAM_DESCRIPTION, Default: 'adata.obj'
#' @return OUTPUT_DESCRIPTION
#' @details DETAILS
#' @examples
#' \dontrun{
#' if(interactive()){
#'  #EXAMPLE1
#'  }
#' }
#' @rdname scVelo.RunBasic
#' @export

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

#' @title FUNCTION_TITLE
#' @description FUNCTION_DESCRIPTION
#' @param adata PARAM_DESCRIPTION
#' @param seu PARAM_DESCRIPTION
#' @param dr PARAM_DESCRIPTION, Default: NULL
#' @return OUTPUT_DESCRIPTION
#' @details DETAILS
#' @examples
#' \dontrun{
#' if(interactive()){
#'  #EXAMPLE1
#'  }
#' }
#' @rdname scVelo.AddDimReduc
#' @export

scVelo.AddDimReduc <- function(adata, seu, dr = NULL) {
  library(reticulate)
  library(rlang)
  dr_list <- Reductions(seu)
  dr <- dr %||% dr_list
  if(!all(dr %in% dr_list)){
    stop("Dimention reduction '", dr, "' not found.\n",
         "Available option(s): ", paste0(dr_list, collapse = ", "))
  }
  code <- paste0("import pickle\n",
                 'pickle_in = open("',adata,'","rb")\n',
                 "adata = pickle.load(pickle_in)\n",
                 "from numpy import genfromtxt")
  py_run_string(code)
  for (i in dr) {
    nc <- ncol(Embeddings(seu, reduction = i))
    export.dr(seu, i)
    code <- paste0(
      "adata.obsm['",i,"_cell_embeddings'] = genfromtxt('tmp/dr.csv', ",
      "delimiter=',', skip_header=1, usecols=range(1,", nc + 1 ,"))"
    )
    py_run_string(code)
  }
  code <- paste0(
    "import scvelo as scv\n",
    "scv.tl.velocity_graph(adata)\n",
    'pickle_out = open("', adata, '","wb")\n',
    "pickle.dump(adata, pickle_out)\n",
    "pickle_out.close()\n"
  )
  py_run_string(code)
}
