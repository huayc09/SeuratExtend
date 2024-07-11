package_data <- list(
  "loomR" = list(
    website = "https://github.com/mojaveazure/loomR",
    tutorial = "https://satijalab.org/loomR/loomR_tutorial.html",
    install_info = 'install.packages("hdf5r")\nremotes::install_github(repo = "mojaveazure/loomR")',
    install = function(){
      remotes::install_github(repo = 'mojaveazure/loomR', ref = 'develop')
    }
  ),
  "nichenetr" = list(
    website = "https://github.com/saeyslab/nichenetr",
    tutorial = "https://github.com/saeyslab/nichenetr",
    install_info = 'remotes::install_github("saeyslab/nichenetr")',
    install = function(){
      remotes::install_github("saeyslab/nichenetr")
    }
  ),
  "biomaRt" = list(
    website = "https://bioconductor.org/packages/release/bioc/html/biomaRt.html",
    tutorial = "https://bioconductor.org/packages/release/bioc/vignettes/biomaRt/inst/doc/biomaRt.html",
    install_info = 'BiocManager::install("biomaRt")',
    install = function(){
      BiocManager::install("biomaRt")
    }
  ),
  "mgsa" = list(
    website = "https://www.bioconductor.org/packages/release/bioc/html/mgsa.html",
    tutorial = "https://www.bioconductor.org/packages/release/bioc/vignettes/mgsa/inst/doc/mgsa.pdf",
    install_info = 'BiocManager::install("mgsa")',
    install = function(){
      BiocManager::install("mgsa")
    }
  ),
  "AUCell" = list(
    website = "https://bioconductor.org/packages/release/bioc/html/AUCell.html",
    tutorial = "https://bioconductor.org/packages/release/bioc/vignettes/AUCell/inst/doc/AUCell.html",
    install_info = 'remotes::install_github("aertslab/AUCell")',
    install = function(){
      remotes::install_github("aertslab/AUCell")
    }
  ),
  "ComplexHeatmap" = list(
    website = "https://bioconductor.org/packages/release/bioc/html/ComplexHeatmap.html",
    tutorial = "https://jokergoo.github.io/ComplexHeatmap-reference/book/",
    install_info = 'BiocManager::install("ComplexHeatmap")',
    install = function(){
      BiocManager::install("ComplexHeatmap")
    }
  ),
    "SeuratDisk" = list(
      website = "https://github.com/mojaveazure/seurat-disk",
      tutorial = "https://github.com/mojaveazure/seurat-disk",
      install_info = 'remotes::install_github("mojaveazure/seurat-disk")',
      install = function(){
        remotes::install_github("mojaveazure/seurat-disk")
      }
  ),
  "DelayedArray" = list(
    website = "https://bioconductor.org/packages/release/bioc/html/DelayedArray.html",
    tutorial = "https://bioconductor.org/packages/release/bioc/html/DelayedArray.html",
    install_info = 'BiocManager::install("DelayedArray")',
    install = function(){
      BiocManager::install("DelayedArray")
    }
  ),
  "DelayedMatrixStats" = list(
    website = "https://bioconductor.org/packages/release/bioc/html/DelayedMatrixStats.html",
    tutorial = "https://bioconductor.org/packages/release/bioc/html/DelayedMatrixStats.html",
    install_info = 'BiocManager::install("DelayedMatrixStats")',
    install = function(){
      BiocManager::install("DelayedMatrixStats")
    }
  ),
  "slingshot" = list(
    website = "https://bioconductor.org/packages/release/bioc/html/slingshot.html",
    tutorial = "https://bioconductor.org/packages/release/bioc/html/slingshot.html",
    install = function(){
      BiocManager::install("slingshot")
    }
  )
)

import <- function(package, hpc.mode = F, detach = F){
  if(!require(package, character.only = T)) {
    message("Required package '", package, "' not installed\n")
    if(package %in% names(package_data)) {
      message("Website: ", package_data[[package]]$website)
      message("Tutorial: ", package_data[[package]]$tutorial)
      message("Recommended installation:\n", package_data[[package]]$install_info,"\n")
    }
    if(!hpc.mode){
      input <- readline(prompt="Try install? y/[n] ")
    } else input <- "y"
    if(input %in% c("y","yes","Y","Yes")){
      if(package %in% names(package_data))
        package_data[[package]]$install() else
          install.packages(package)
    }
  }
  library(package, character.only = T)
  if(detach) detach(paste0("package:",package), unload = TRUE, character.only = T)
}

