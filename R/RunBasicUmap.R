#' @title FUNCTION_TITLE
#' @description FUNCTION_DESCRIPTION
#' @param path PARAM_DESCRIPTION, Default: getwd()
#' @param project PARAM_DESCRIPTION
#' @param sub.path PARAM_DESCRIPTION, Default: '/outs/filtered_feature_bc_matrix'
#' @param saveRDS PARAM_DESCRIPTION, Default: T
#' @param save.path PARAM_DESCRIPTION, Default: 'rds'
#' @param spe PARAM_DESCRIPTION, Default: getOption("spe")
#' @param min.cells PARAM_DESCRIPTION, Default: 1
#' @param min.features PARAM_DESCRIPTION, Default: 200
#' @param dims PARAM_DESCRIPTION, Default: c(1:20)
#' @param resolution PARAM_DESCRIPTION, Default: 0.5
#' @return OUTPUT_DESCRIPTION
#' @details DETAILS
#' @examples
#' \dontrun{
#' if(interactive()){
#'  #EXAMPLE1
#'  }
#' }
#' @rdname RunBasicUmap
#' @export

RunBasicUmap <-
  function(path = getwd(), project, sub.path = "/outs/filtered_feature_bc_matrix", saveRDS = T, save.path = "rds",
           spe = getOption("spe"), min.cells = 1, min.features = 200, dims = c(1:20), resolution = 0.5){
    library(Seurat)
    library(dplyr)
    check_spe(spe)
    Seu<-Read10X(paste0(path,sub.path))
    Seu<-CreateSeuratObject(Seu, project = project, min.cells = min.cells, min.features = min.features)
    mt.pattern = c("human" = "^MT-", "mouse" = "^mt-")
    Seu[["percent.mt"]] <- PercentageFeatureSet(Seu, pattern = mt.pattern[spe])
    Seu <- Seu %>%
      subset(subset = percent.mt < 12.5) %>%
      NormalizeData() %>%
      FindVariableFeatures() %>%
      ScaleData() %>%
      RunPCA() %>%
      FindNeighbors(dims = dims) %>%
      FindClusters(resolution = resolution) %>%
      RunUMAP(dims = dims)
    if(saveRDS){
      if(!save.path %in% list.files()) dir.create(save.path)
      filename <- paste0(save.path, "/SeuratObj_", project, "_min.feature.", min.features, ".rds")
      saveRDS(Seu, filename)
    }
    return(Seu)
  }

# setwd("~/R documents/2020-3-15 AUID")
# list.files("CellRanger")
# # [1] "CAPS1_ZJB_NT_180065C"     "CAPS1_ZJB_TNFRi_190862A"  "CAPS2_NT_190862B"         "HC_180065B"
# # [5] "TRAPS1_WYH_NT_180065A"    "TRAPS1_WYH_TNFRi_180065D"
# path <- paste0(getwd(), "/CellRanger/",list.files("CellRanger"))[1]
# project = "project"
# sub.path = "/outs/filtered_feature_bc_matrix"
# saveRDS = T
# save.path = "rds"
# spe = "human"
# min.cells = 1
# min.features = 200
# percent.mt = 15
# dims = c(1:20)
# resolution = 0.5
