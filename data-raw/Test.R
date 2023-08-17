library(SeuratExtend)
library(Seurat)
library(dplyr)
library(mosaic)
library(rlist)
library(purrr)
library(roxygen2)
library(sinew)
options(max.print = 100, spe = "human", nCores = 12)

usethis::use_data(PanglaoDB_data, overwrite = TRUE)
usethis::use_gpl_license(version = 3, include_future = TRUE)
rmarkdown::render("vignettes/CalcStats.Rmd", output_format = "md_document")
makeOxygen(VlnPlot2.default)
roxygenize()


