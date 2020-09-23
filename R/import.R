package_data <- list(
  "loomR" = list(
    website = "https://github.com/mojaveazure/loomR",
    tutorial = "https://satijalab.org/loomR/loomR_tutorial.html",
    install_info = 'install.packages("hdf5r")\ndevtools::install_github(repo = "mojaveazure/loomR")',
    install = function(){
      install.packages("hdf5r")
      devtools::install_github(repo = 'mojaveazure/loomR', ref = 'develop')
    }
  ),
  "nichenetr" = list(
    website = "https://github.com/saeyslab/nichenetr",
    tutorial = "https://github.com/saeyslab/nichenetr",
    install_info = 'devtools::install_github("saeyslab/nichenetr")',
    install = function(){
      devtools::install_github("saeyslab/nichenetr")
    }
  ),
  "biomaRt" = list(
    website = "https://bioconductor.org/packages/release/bioc/html/biomaRt.html",
    tutorial = "https://bioconductor.org/packages/release/bioc/vignettes/biomaRt/inst/doc/biomaRt.html",
    install_info = 'BiocManager::install("biomaRt")',
    install = function(){
      BiocManager::install("biomaRt")
    }
  )
)

#' @title FUNCTION_TITLE
#' @description FUNCTION_DESCRIPTION
#' @param package PARAM_DESCRIPTION
#' @param hpc.mode PARAM_DESCRIPTION, Default: F
#' @return OUTPUT_DESCRIPTION
#' @details DETAILS
#' @examples
#' \dontrun{
#' if(interactive()){
#'  #EXAMPLE1
#'  }
#' }
#' @rdname import
#' @export

import <- function(package, hpc.mode = F){
  if(!require(package, character.only = T)) {
    message(paste0("Required package '", package, "' not installed\n"))
    message(paste0("Website: ", package_data[[package]]$website))
    message(paste0("Tutorial: ", package_data[[package]]$tutorial))
    message(paste0("Recommended installation:\n", package_data[[package]]$install_info),"\n")
    if(!hpc.mode){
      input <- readline(prompt="Try install? y/[n] ")
    } else input <- "y"
    if(input %in% c("y","yes","Y","Yes")) package_data[[package]]$install()
  }
  library(package, character.only = T)
}

