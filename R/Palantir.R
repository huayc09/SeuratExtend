#' @title FUNCTION_TITLE
#' @description FUNCTION_DESCRIPTION
#' @param seu PARAM_DESCRIPTION
#' @param reduction PARAM_DESCRIPTION, Default: 'pca'
#' @param n_components PARAM_DESCRIPTION, Default: 20
#' @return OUTPUT_DESCRIPTION
#' @details DETAILS
#' @examples
#' \dontrun{
#' if(interactive()){
#'  #EXAMPLE1
#'  }
#' }
#' @rdname RunPalantirDiffusionMap
#' @export

RunPalantirDiffusionMap <- function(seu, reduction = "pca", n_components = 20) {
  library(Seurat)
  library(reticulate)

  if(!"tmp" %in% list.files()) dir.create("tmp")
  write.csv(Embeddings(seu, reduction = reduction), "tmp/dr.csv")

  code <- paste0("import palantir\n",
                 "import os\n",
                 "import matplotlib\n",
                 "import matplotlib.pyplot as plt\n",
                 "import pandas as pd\n",
                 "import numpy as np\n",
                 "import random\n\n",
                 "pca_projections = pd.read_csv('tmp/dr.csv', index_col=0)\n",
                 "dm_res = palantir.utils.run_diffusion_maps(pca_projections, n_components=",
                 n_components,
                 ")\n",
                 "dm_res['EigenVectors'].to_csv('tmp/dm_res.csv')\n",
                 "ms_data = palantir.utils.determine_multiscale_space(dm_res)\n",
                 "ms_data.to_csv('tmp/ms_data.csv')\n",
                 "tsne = palantir.utils.run_tsne(ms_data)\n",
                 "tsne.to_csv('tmp/tsne.csv')\n",
                 "quit")
  py_run_string(code)
  seu[["dm"]] <-
    read.csv("tmp/dm_res.csv", row.names = 1) %>%
    set_colnames(paste0("DM_", 1:ncol(.))) %>%
    as.matrix() %>%
    CreateDimReducObject(key = "DM_", assay = DefaultAssay(seu))
  seu[["ms"]] <-
    read.csv("tmp/ms_data.csv", row.names = 1) %>%
    set_colnames(paste0("MS_", 1:ncol(.))) %>%
    as.matrix() %>%
    CreateDimReducObject(key = "MS_", assay = DefaultAssay(seu))
  seu[["tsne"]] <-
    read.csv("tmp/tsne.csv", row.names = 1) %>%
    set_colnames(paste0("TSNE_", 1:ncol(.))) %>%
    as.matrix() %>%
    CreateDimReducObject(key = "TSNE_", assay = DefaultAssay(seu))
  return(seu)
}
