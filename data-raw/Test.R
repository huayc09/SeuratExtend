# seuratObj <- readRDS("~/Documents/scRNA/HEV-seurat3/Robjects/seuratObj.rds")
# SeuObjList <- readRDS("~/Documents/scRNA/2020-2-10 EC PyMT and E0771/rds/SeuObjList.rds")

library(SeuratExtend)
library(Seurat)
library(dplyr)
library(mosaic)
library(rlist)
library(purrr)
library(roxygen2)
library(sinew)
options(max.print = 50, spe = "mouse", nCores = 12)
seu <- readRDS("~/R documents/2020-2-10 EC PyMT and E0771/rds/PyMTEC_old.rds")

seu <- seuratObj
matr <- as.matrix(GetAssayData(seuratObj))[1:20,]
meta <- seuratObj@meta.data
f <- seuratObj@meta.data$cluster
group.by = "cluster"
method <- "tscore"
assay = "RNA"
features = c("Sele","Selp","asdfsdf", "asdfsdff")
species = "Mouse"
auto.cluster = "EC"
name.orig = "seurat_clusters"
Seu <- SeuObjList[[1]]
name.new = "auto_cluster"
dataset = "GO:1901342"
root = "BP"
spe = "mouse"
ratio = 0.2
n.min = 0
n.max = Inf
only.end.terms = T
slot = "counts"
assay = "RNA"
nCores = 4
makeOxygen(HumanToMouseGenesymbol)
roxygenize()
CalcScoreGeneral_v3(Seu, features, group.by, "zscore")

