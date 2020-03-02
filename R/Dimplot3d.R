#' @title FUNCTION_TITLE
#' @description FUNCTION_DESCRIPTION
#' @param Seu PARAM_DESCRIPTION
#' @param group.by PARAM_DESCRIPTION, Default: 'seurat_clusters'
#' @param assay PARAM_DESCRIPTION, Default: NULL
#' @param reduction PARAM_DESCRIPTION, Default: 'pca'
#' @param dims PARAM_DESCRIPTION, Default: NULL
#' @param text_shift PARAM_DESCRIPTION, Default: 0.1
#' @param size PARAM_DESCRIPTION, Default: 16
#' @param width PARAM_DESCRIPTION, Default: 640
#' @param legend3d PARAM_DESCRIPTION, Default: T
#' @param writeWebGL PARAM_DESCRIPTION, Default: F
#' @param filename PARAM_DESCRIPTION, Default: NULL
#' @param slingPseudotime PARAM_DESCRIPTION, Default: NULL
#' @return OUTPUT_DESCRIPTION
#' @details DETAILS
#' @examples
#' \dontrun{
#' if(interactive()){
#'  #EXAMPLE1
#'  }
#' }
#' @rdname Dimplot3d
#' @export
Dimplot3d <- function(Seu, group.by = "seurat_clusters", assay = NULL, reduction = "pca",
                      dims = NULL, text_shift = 0.1, size = 16, width = 640,
                      legend3d = T, writeWebGL = F, filename = NULL, slingPseudotime = NULL){
  library(rgl)
  library(rlang)
  library(rlist)
  DefaultAssay(Seu) <- assay %||% DefaultAssay(Seu)
  Command <- paste("FindNeighbors",DefaultAssay(Seu),reduction, sep = ".")
  dims <- dims %||% Seu@commands[[Command]]$dims
  Seu <- RunUMAP(Seu, dims = dims, n.components = 3, reduction = reduction)
  f <- factor(Seu@meta.data[,group.by])
  Textdf <-
    Embeddings(Seu, reduction = "umap") %>%
    as.data.frame() %>%
    split(f) %>%
    lapply(function(x)apply(x, 2, mean)) %>%
    list.rbind()
  cols = gg_color_hue(nlevels(f))
  par3d(windowRect = 50 + c(0, 0, width, width))
  plot3d(
    Embeddings(Seu, reduction = "umap"),
    col = cols[f], type = 'p', alpha = 0.2)
  points3d(Textdf, size = size, alpha = 1, col = cols)
  text3d(Textdf + text_shift,
         texts = rownames(Textdf),
         adj = 0)
  check.slingPseudotime <-
    identical(Seu@meta.data[,group.by], Seu@misc[["slingshot"]][[toupper(reduction)]]$f)
  if(slingPseudotime %||% check.slingPseudotime){
    lapply(Seu@misc[["slingshot"]][[toupper(reduction)]]$SlingshotDataSet@lineages,
           function(x){lines3d(Textdf[x,], lwd = 8, alpha = 0.6)})
  }
  if(legend3d){
    legend3d("topright", legend = levels(f), pch = 16,
             col = cols, cex=1, inset=c(0.02))
  }
  if(writeWebGL){
    filename <- filename %||% paste0(paste(group.by, reduction, "dim", max(dims), sep = "_"), ".html")
    writeWebGL(dir = getwd(), filename = filename)
  }
}

# Dimplot3d(Seu, group.by = "cluster")
# Dimplot3d(Seu)
