setwd("~/Documents/SeuratExtend-database/2024-4-5 color pro from i want hue")

files <- list.files(pattern = ".csv")
files <- files[!grepl("_original",files)]
theme <- c("default","light","red","yellow","green","blue","purple","bright")

lines <- readLines(files[1])
color_palettes <- strsplit(lines, ",")
color_palettes <- split(color_palettes, f = lengths(color_palettes))

color_pro_presets <- list()
for (i in theme) {
  file_diff <- paste0("diffSortedPalettes_pro_",i,".csv")
  file_hue <- paste0("hueSortedPalettes_pro_",i,".csv")
  file <- c(diff = file_diff, hue = file_hue)
  for (j in c("hue","diff")) {
    lines <- readLines(file[j])
    color_palettes <- strsplit(lines, ",")
    color_pro_presets[[i]][[j]] <- split(color_palettes, f = lengths(color_palettes))
  }
}

presets_color_pro <- color_pro_presets
usethis::use_data(presets_color_pro, overwrite = TRUE)
