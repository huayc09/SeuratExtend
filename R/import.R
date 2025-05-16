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

import <- function(package, hpc.mode = FALSE, detach = FALSE){
  # Suppress warnings from require - we'll handle messaging ourselves
  is_installed <- suppressWarnings(require(package, character.only = TRUE, quietly = TRUE))
  
  if(!is_installed) {
    message("Required package '", package, "' not installed")
    if(package %in% names(package_data)) {
      message("Website: ", package_data[[package]]$website)
      message("Tutorial: ", package_data[[package]]$tutorial)
      message("Recommended installation:\n", package_data[[package]]$install_info,"\n")
    }
    
    # Ask if user wants to install
    if(!hpc.mode){
      input <- readline(prompt="Try install? y/[n] ")
    } else input <- "y"
    
    if(input %in% c("y","yes","Y","Yes")){
      message("Attempting to install package '", package, "'...")
      
      # Try to install the package
      tryCatch({
        if(package %in% names(package_data)) {
          package_data[[package]]$install()
        } else {
          install.packages(package)
        }
        
        # Check if installation was successful
        if(suppressWarnings(require(package, character.only = TRUE, quietly = TRUE))) {
          message("Package '", package, "' installed and loaded successfully")
          
          # If detach was requested, detach after loading
          if(detach) {
            detach(paste0("package:", package), unload = TRUE, character.only = TRUE)
            return(TRUE)
          }
          return(TRUE)
        } else {
          message("Package installation appeared to complete but the package could not be loaded")
          message("Please try installing manually using the command above")
          return(FALSE)
        }
      }, error = function(e) {
        message("Error during installation: ", e$message)
        message("Please try installing manually using the command above")
        return(FALSE)
      })
    } else {
      message("Package installation skipped")
      return(FALSE)
    }
  } else {
    # Package is already installed and loaded
    if(detach) {
      detach(paste0("package:", package), unload = TRUE, character.only = TRUE)
    }
    return(TRUE)
  }
}

