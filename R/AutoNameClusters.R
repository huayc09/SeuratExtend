#' @title FUNCTION_TITLE
#' @description FUNCTION_DESCRIPTION
#' @param Seu PARAM_DESCRIPTION
#' @param auto.cluster PARAM_DESCRIPTION, Default: 'EC'
#' @param species PARAM_DESCRIPTION, Default: 'Mouse'
#' @param name.orig PARAM_DESCRIPTION, Default: 'seurat_clusters'
#' @param name.new PARAM_DESCRIPTION, Default: 'auto_cluster'
#' @return OUTPUT_DESCRIPTION
#' @details DETAILS
#' @examples
#' \dontrun{
#' if(interactive()){
#'  #EXAMPLE1
#'  }
#' }
#' @rdname AutoNameClusters
#' @export
AutoNameClusters <- function(Seu, auto.cluster = "EC", species = "Mouse",
                             name.orig = "seurat_clusters", name.new = "auto_cluster"){
  if(species == "Mouse") ref <- AutoCluster_Mouse[[auto.cluster]]
  f <- factor(Seu@meta.data[[name.orig]])
  Exp <-
    CalcScoreGeneral_v3(Seu, features = ref$Gene, group.by = name.orig, method = "zscore") %>%
    as.data.frame %>%
    mutate(cluster = ref$Cluster, Max.exp = colnames(.)[apply(., 1, which.max)])
  Seu@meta.data[[name.new]] <-
    split(Exp$cluster, Exp$Max.exp) %>%
    lapply(function(x) paste0(x, collapse = ", ")) %>%
    unlist %>%
    c(setdiff(levels(f), Exp$Max.exp) %>% setNames(.,.)) %>%
    .[as.vector(f)]
  return(Seu)
}
