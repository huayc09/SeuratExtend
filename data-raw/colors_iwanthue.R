library(openxlsx)
library(tidyr)
library(dplyr)
file <- "d:/Documents/R documents/SeuratExtend_databases/2021-5-21 preset colors/2021-8-12 preset colors I want hue.xlsx"
col.space <- c("default","intense","pastel","all","all_hard")
color_presets <- list()
for (i in 1:5) {
  col_list <- read.xlsx(file, sheet = i, colNames = F) %>%
    as.list() %>%
    lapply(function(x) x[!is.na(x)]) %>%
    split(f = lengths(.)) %>%
    lapply(function(x) {
      names(x) <- seq_along(x)
      x})
  color_presets[[col.space[i]]] <- col_list
}

usethis::use_data(color_presets, overwrite = TRUE)
# rename
presets_color_iwh <- color_presets
usethis::use_data(presets_color_iwh, overwrite = TRUE)
