library(SeuratExtend)
library(Seurat)
library(dplyr)
library(mosaic)
library(rlist)
library(purrr)
library(roxygen2)
library(sinew)
options(max.print = 50, spe = "mouse", nCores = 12)

usethis::use_data(mouse_human_genesymbols, overwrite = TRUE)
makeOxygen(write.py.fun)
roxygenize()


