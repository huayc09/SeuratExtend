library(Seurat)
pbmc <- Read10X("pbmc3k_10X/outs/filtered_feature_bc_matrix/")
pbmc <- CreateSeuratObject(counts = pbmc, project = "pbmc3k", min.cells = 3, min.features = 200)
pbmc[["percent.mt"]] <- PercentageFeatureSet(pbmc, pattern = "^MT-")
pbmc <- subset(pbmc, subset = nFeature_RNA > 200 & nFeature_RNA < 2500 & percent.mt < 5)
pbmc <- NormalizeData(pbmc)
pbmc <- FindVariableFeatures(pbmc)
pbmc <- ScaleData(pbmc)
pbmc <- RunPCA(pbmc)
pbmc <- FindNeighbors(pbmc, dims = 1:10)
pbmc <- FindClusters(pbmc, resolution = 0.5)
pbmc <- RunUMAP(pbmc, dims = 1:10)
DimPlot(pbmc, reduction = "umap")
new.cluster.ids <- c("CD4 T Naive", "Mono CD14", "CD4 T Memory", "B cell", "CD8 T cell", "Mono FCGR3A", "NK cell", "DC", "Platelet")
pbmc$cluster <- new.cluster.ids[pbmc$seurat_clusters]
DimPlot(pbmc, group.by = "cluster", label = T)
# remove
cells <- CellSelector(DimPlot(pbmc))
cell.rm1 <- intersect(cells, colnames(pbmc)[pbmc$seurat_clusters == "0"])
cells <- CellSelector(DimPlot(pbmc))
cell.rm2 <- setdiff(cells, colnames(pbmc)[pbmc$seurat_clusters == "8"])
cell.rm <- c(cell.rm1, cell.rm2)
# keep enough number of DC and platelet
table(pbmc$seurat_clusters)
cell.kp <- colnames(pbmc)[pbmc$seurat_clusters %in% c("8")]
cell.kp <- c(cell.kp, sample(colnames(pbmc)[pbmc$seurat_clusters %in% c("7")], 12))
DimPlot(pbmc, cells.highlight = cell.rm)
DimPlot(pbmc, cells.highlight = cell.kp)

cells <- sample(setdiff(colnames(pbmc), c(cell.rm, cell.kp)), 500 - length(cell.kp))
cells <- c(cells, cell.kp)
DimPlot(pbmc, cells.highlight = cells)

pbmc_sub <- subset(pbmc, cells = cells)
pbmc_sub <- subset(pbmc_sub, features = rownames(pbmc_sub)[rowSums(GetAssayData(pbmc_sub)) >= 1] )
DimPlot(pbmc_sub)
table(pbmc_sub$cluster)
pbmc_sub@meta.data$orig.ident <- sample(c("sample1","sample2","sample2"), size = 500, replace = T)
pbmc_sub@meta.data$orig.ident <- factor(pbmc_sub$orig.ident)
pbmc_sub$cluster <- factor(pbmc_sub$cluster)
Idents(pbmc_sub) <- 'cluster'
DimPlot(pbmc_sub, label = T)
DimPlot(pbmc_sub, group.by = "orig.ident")
pbmc <- pbmc_sub

set.seed(42)

pbmc$condition <- gsub("sample","condition",pbmc$orig.ident)
pbmc$sample_id <- paste0(
  pbmc$condition, "_rep",
  sample(1:3, ncol(pbmc), replace = TRUE)
)

usethis::use_data(pbmc, overwrite = TRUE)

