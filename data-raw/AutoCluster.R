## code to prepare `mydataset` dataset goes here

library(readxl)
AutoCluster_Mouse <- list()
AutoCluster_Mouse[["EC"]] <-
  read_excel("data-raw/AutoCluster_Mouse.xlsx", sheet = "EC") %>%
  as.data.frame()

usethis::use_data(AutoCluster_Mouse, overwrite = TRUE)




