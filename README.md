# SeuratExtend: An Enhanced Toolkit for scRNA-seq Analysis

## Overview

`SeuratExtend` is an R package designed to provide an improved and easy-to-use toolkit for scRNA-seq analysis and visualization, built upon the Seurat object. While `Seurat` is a widely-used tool in the R community that offers a foundational framework for scRNA-seq analysis, it has limitations when it comes to more advanced analysis and customized visualization. `SeuratExtend` expands upon `Seurat` by offering an array of enhanced visualization tools, an integrated functional and pathway analysis pipeline, seamless integration with popular Python tools, and a suite of utility functions for data manipulation and presentation. Designed to be user-friendly even for beginners, the package retains a level of professionalism that ensures rigorous analysis.

**Key Features**:

- **Enhanced Data Visualization**: Includes heatmaps, violin plots, dimensional reduction (UMAP) plots, waterfall plots, dot plots, proportion bars, volcano plots, and GSEA plots.
- **Integrated Functional and Pathway Analysis**: Supports GO and Reactome databases, with the option to use custom databases.
- **Python Tool Integration**: Easily apply tools like scVelo, SCENIC, and Palantir within R using the Seurat object.
- **Utility Functions**: Assorted functions for calculations and color selections to streamline your scRNA-seq analysis.

## Resources

- **GitHub Repository**: Access the source code and contribute to SeuratExtend on [GitHub](https://github.com/huayc09/SeuratExtend).
- **Online Tutorial**: For a comprehensive guide on using SeuratExtend, visit our [tutorial website](https://huayc09.github.io/SeuratExtend/).
- **SeuratExtend Chatbot**: Try our AI-powered assistant (beta version, powered by ChatGPT) for help with scRNA-seq analysis: [scRNA-seq Assistant](https://chatgpt.com/g/g-8scQjmzkd-scrna-seq-assistant).

## Citation

If you use SeuratExtend in your research, please cite:

Hua, Y., Weng, L., Zhao, F., and Rambow, F. (2025). SeuratExtend: streamlining single-cell RNA-seq analysis through an integrated and intuitive framework. Gigascience 14, giaf076. https://doi.org/10.1093/gigascience/giaf076.

## Installation

Install `SeuratExtend` directly from GitHub:

```R
if (!requireNamespace("remotes", quietly = TRUE)) {
    install.packages("remotes")
}
remotes::install_github("huayc09/SeuratExtend")
```

## Vignettes and Tutorials

### [Quick Start-Up Guide](#quick-start-up-guide-1)

### [What's New in v1.2.0](vignettes/News.md)

### [Enhanced Visualization](vignettes/Visualization.md)
- [Create an Enhanced Dimensional Reduction Plot](vignettes/Visualization.md#create-an-enhanced-dimensional-reduction-plot) `DimPlot2` `FeaturePlot3` `FeaturePlot3.grid` `theme_umap_arrows`
- [Generate a Heatmap Plot](vignettes/Visualization.md#generate-a-heatmap-plot) `Heatmap`
- [Create Enhanced Dot Plots](vignettes/Visualization.md#create-enhanced-dot-plots-new-in-v110) **(New in v1.1.0)** `DotPlot2`
- [Create an Enhanced Violin Plot](vignettes/Visualization.md#create-an-enhanced-violin-plot) `VlnPlot2`
- [Visualize Cluster Distribution in Samples](vignettes/Visualization.md#visualize-cluster-distribution-in-samples) `ClusterDistrBar`
- [Generate a Waterfall Plot](vignettes/Visualization.md#generate-a-waterfall-plot) `WaterfallPlot`
- [Create Volcano Plots](vignettes/Visualization.md#create-volcano-plots-new-in-v110) **(New in v1.1.0)** `VolcanoPlot`
- [Explore Color Functions](vignettes/Visualization.md#explore-color-functions) `color_pro` `color_iwh` `ryb2rgb` `save_colors`

### [Geneset Enrichment Analysis (GSEA)](vignettes/GSEA.md)
- [Conduct GSEA using the GO or Reactome database](vignettes/GSEA.md#conduct-gsea-using-the-go-or-reactome-database) `GeneSetAnalysisGO` `GeneSetAnalysisReactome`
- [Perform GSEA using customized genesets](vignettes/GSEA.md#perform-gsea-using-customized-genesets) `GeneSetAnalysis`
- [Find pathways in the GO/Reactome database or customized genesets](vignettes/GSEA.md#find-pathways-in-the-goreactome-database-or-customized-genesets) `SearchDatabase` `SearchPathways`
- [Convert GO/Reactome pathway IDs to pathway names](vignettes/GSEA.md#convert-goreactome-pathway-ids-to-pathway-names) `RenameGO` `RenameReactome`
- [Filter the GO/Reactome pathway list based on certain criteria](vignettes/GSEA.md#filter-the-goreactome-pathway-list-based-on-certain-criteria) `FilterGOTerms` `FilterReactomeTerms`
- [Create a GSEA plot emulating the Broad Institute analysis](vignettes/GSEA.md#create-a-gsea-plot-emulating-the-broad-institute-analysis) `GSEAplot`

### [Trajectory and Pseudotime Analysis](vignettes/Trajectory.md)
-  [scVelo Tutorial for Trajectory Analysis](vignettes/Trajectory.md#analyzing-single-cell-trajectories-with-scvelo) `scVelo.SeuratToAnndata` `scVelo.Plot` 
-  [Palantir Tutorial for Trajectory and Pseudotime Analysis](vignettes/Trajectory.md#palantir-tutorial-for-trajectory-and-pseudotime-analysis) `Palantir.RunDM` `Palantir.Pseudotime`
-  [MAGIC for Denoising and Smoothing Gene Expression](vignettes/Trajectory.md#magic-for-denoising-and-smoothing-gene-expression) `Palantir.Magic`
-  [CellRank Tutorial for Trajectory Analysis](vignettes/Trajectory.md#cellrank-tutorial-for-trajectory-analysis) `Cellrank.Compute` `Cellrank.Plot`
-  [Gene Expression Dynamics Along Differentiation Trajectories](vignettes/Trajectory.md#gene-expression-dynamics-along-differentiation-trajectories) `GeneTrendCurve.Palantir` `GeneTrendHeatmap.Palantir` `GeneTrendCurve.Slingshot` `GeneTrendHeatmap.Slingshot`
-  [Slingshot Tutorial for Pseudotime Analysis](vignettes/Trajectory.md#slingshot-tutorial-for-pseudotime-analysis) `RunSlingshot` 
-  [Integration of Seurat with Python Tools](vignettes/Trajectory.md#integration-of-seurat-with-python-tools) `create_condaenv_seuratextend` `Seu2Adata` `Seu2Loom` `adata.LoadLoom` `adata.AddDR` `adata.AddMetadata` `adata.Save` `adata.Load`

### [SCENIC for Gene Regulatory Networks Analysis](vignettes/SCENIC.md)
- [Importing SCENIC Loom Files into Seurat](vignettes/SCENIC.md#importing-scenic-loom-files-into-seurat) `ImportPyscenicLoom`
- [Visualizing SCENIC Results](vignettes/SCENIC.md#visualizing-scenic-results) 

### [Utility Tools and Functions](vignettes/Utilities.md)
- [Facilitate Gene Naming Conversions](vignettes/Utilities.md#facilitate-gene-naming-conversions) `HumanToMouseGenesymbol` `MouseToHumanGenesymbol` `EnsemblToGenesymbol` `GenesymbolToEnsembl` `UniprotToGenesymbol`
- [Compute Statistics Grouped by Clusters](vignettes/Utilities.md#compute-statistics-grouped-by-clusters) `CalcStats`
- [Assess Proportion of Positive Cells in Clusters](vignettes/Utilities.md#assess-proportion-of-positive-cells-in-clusters) `feature_percent`
- [Run Standard Seurat Pipeline](vignettes/Utilities.md#run-standard-seurat-pipeline) `RunBasicSeurat`

### [Single-Cell RNA-seq Analysis Course](#single-cell-rna-seq-analysis-course-new-in-v110-1) **(New in v1.1.0)**

### [FAQ](vignettes/FAQ.md)

## Quick Start-Up Guide

This quick start-up guide provides an overview of the most frequently
used functions in single-cell RNA sequencing (scRNA-seq) analysis. After
running the standard Seurat pipeline (refer to this [Seurat pbmc3k
tutorial](https://satijalab.org/seurat/articles/pbmc3k_tutorial)), you
should have a Seurat object ready for further analysis. Below, we
illustrate the use of a subset of the pbmc dataset as an example to
demonstrate various functionalities of the `SeuratExtend` package.

### Visualizing Clusters

``` r
library(Seurat)
library(SeuratExtend)

# Visualizing cell clusters using DimPlot2
DimPlot2(pbmc, theme = theme_umap_arrows())
```

![](vignettes/quick_start_files/figure-gfm/unnamed-chunk-1-1.png)<!-- -->

### Analyzing Cluster Distribution

To check the percentage of each cluster within different samples:

``` r
# Cluster distribution bar plot
ClusterDistrBar(pbmc$orig.ident, pbmc$cluster)
```

![](vignettes/quick_start_files/figure-gfm/unnamed-chunk-2-1.png)<!-- -->

### Marker Gene Analysis with Heatmap

To examine the marker genes of each cluster and visualize them using a
heatmap:

``` r
# Calculating z-scores for variable features
genes.zscore <- CalcStats(
  pbmc,
  features = VariableFeatures(pbmc),
  group.by = "cluster",
  order = "p",
  n = 4)
  
# Displaying heatmap
Heatmap(genes.zscore, lab_fill = "zscore")
```

![](vignettes/quick_start_files/figure-gfm/unnamed-chunk-3-1.png)<!-- -->

### Enhanced Dot Plots (New in v1.1.0)

``` r
# Create grouped features
grouped_features <- list(
  "B_cell_markers" = c("MS4A1", "CD79A"),
  "T_cell_markers" = c("CD3D", "CD8A", "IL7R"),
  "Myeloid_markers" = c("CD14", "FCGR3A", "S100A8")
)

DotPlot2(pbmc, features = grouped_features)
```

![](vignettes/quick_start_files/figure-gfm/unnamed-chunk-4-1.png)<!-- -->

### Enhanced Visualization of Marker Genes

For visualizing specific markers via a violin plot that incorporates box
plots, median lines, and performs statistical testing:

``` r
# Specifying genes and cells of interest
genes <- c("CD3D", "CD14", "CD79A")
cells <- WhichCells(pbmc, idents = c("B cell", "CD8 T cell", "Mono CD14"))

# Violin plot with statistical analysis
VlnPlot2(
  pbmc,
  features = genes,
  group.by = "cluster",
  cells = cells,
  stat.method = "wilcox.test")
```

![](vignettes/quick_start_files/figure-gfm/unnamed-chunk-5-1.png)<!-- -->

### Visualizing Multiple Markers on UMAP

Displaying three markers on a single UMAP, using RYB coloring for each
marker:

``` r
FeaturePlot3(pbmc, feature.1 = "CD3D", feature.2 = "CD14", feature.3 = "CD79A", pt.size = 1)
```

![](vignettes/quick_start_files/figure-gfm/unnamed-chunk-6-1.png)<!-- -->

### Create Volcano Plots (New in v1.1.0)

Create a basic volcano plot comparing two cell types:

``` r
VolcanoPlot(pbmc, 
            ident.1 = "B cell",
            ident.2 = "CD8 T cell")
```

![](vignettes/quick_start_files/figure-gfm/unnamed-chunk-7-1.png)<!-- -->

### Conducting Geneset Enrichment Analysis (GSEA)

Examining all the pathways of the immune process in the Gene Ontology
(GO) database, and visualizing by a heatmap that displays the top
pathways of each cluster across multiple cell types:

``` r
options(spe = "human")
pbmc <- GeneSetAnalysisGO(pbmc, parent = "immune_system_process", n.min = 5)
matr <- RenameGO(pbmc@misc$AUCell$GO$immune_system_process)
go_zscore <- CalcStats(
  matr,
  f = pbmc$cluster,
  order = "p",
  n = 3)
Heatmap(go_zscore, lab_fill = "zscore")
```

![](vignettes/quick_start_files/figure-gfm/unnamed-chunk-8-1.png)<!-- -->

### Detailed Comparison of Two Cell Types

Using a GSEA plot to focus on a specific pathway for deeper comparative
analysis:

``` r
GSEAplot(
  pbmc,
  ident.1 = "B cell",
  ident.2 = "CD8 T cell",
  title = "GO:0042113 B cell activation (335g)",
  geneset = GO_Data$human$GO2Gene[["GO:0042113"]])
```

![](vignettes/quick_start_files/figure-gfm/unnamed-chunk-9-1.png)<!-- -->

### Importing and Visualizing SCENIC Analysis

After conducting Gene Regulatory Networks Analysis using pySCENIC,
import the output and visualize various aspects within Seurat:

``` r
# Downloading a pre-computed SCENIC loom file
scenic_loom_path <- file.path(tempdir(), "pyscenic_integrated-output.loom")
download.file("https://zenodo.org/records/10944066/files/pbmc3k_small_pyscenic_integrated-output.loom", scenic_loom_path, mode = "wb")

# Importing SCENIC Loom Files into Seurat
pbmc <- ImportPyscenicLoom(scenic_loom_path, seu = pbmc)

# Visualizing variables such as cluster, gene expression, and SCENIC regulon activity with customized colors
DimPlot2(
  pbmc,
  features = c("cluster", "orig.ident", "CEBPA", "tf_CEBPA"),
  cols = list("tf_CEBPA" = "OrRd"),
  theme = NoAxes()
) + theme_umap_arrows()
```

![](vignettes/quick_start_files/figure-gfm/unnamed-chunk-10-1.png)<!-- -->

``` r
# Creating a waterfall plot to compare regulon activity between cell types
DefaultAssay(pbmc) <- "TF"
WaterfallPlot(
  pbmc,
  features = rownames(pbmc),
  ident.1 = "Mono CD14",
  ident.2 = "CD8 T cell",
  exp.transform = FALSE,
  top.n = 20)
```

![](vignettes/quick_start_files/figure-gfm/unnamed-chunk-11-1.png)<!-- -->

### Trajectory Analysis with Palantir in R

Trajectory analysis helps identify developmental pathways and
transitions between different cell states. In this section, we
demonstrate how to perform trajectory analysis using the Palantir
algorithm on a subset of myeloid cells, integrating everything within
the R environment.

#### Download and Prepare the Data

First, we download a small subset of myeloid cells to illustrate the
analysis:

``` r
# Download the example Seurat Object with myeloid cells
mye_small <- readRDS(url("https://zenodo.org/records/10944066/files/pbmc10k_mye_small_velocyto.rds", "rb"))
```

#### Diffusion Map Calculation

Palantir uses diffusion maps for dimensionality reduction to infer
trajectories. Here's how to compute and visualize them:

``` r
# Compute diffusion map
mye_small <- Palantir.RunDM(mye_small)
```

    ## Determing nearest neighbor graph...

``` r
# Visualize the first two diffusion map dimensions
DimPlot2(mye_small, reduction = "ms")
```

![](vignettes/quick_start_files/figure-gfm/unnamed-chunk-13-1.png)<!-- -->

#### Pseudotime Calculation

Pseudotime ordering assigns each cell a time point in a trajectory,
indicating its progression along a developmental path:

``` r
# Calculate pseudotime with a specified start cell
mye_small <- Palantir.Pseudotime(mye_small, start_cell = "sample1_GAGAGGTAGCAGTACG-1")
```

    ## Sampling and flocking waypoints...
    ## Time for determining waypoints: 0.00112607479095459 minutes
    ## Determining pseudotime...
    ## Shortest path distances using 30-nearest neighbor graph...
    ## Time for shortest paths: 0.014574062824249268 minutes
    ## Iteratively refining the pseudotime...
    ## Correlation at iteration 1: 1.0000
    ## Entropy and branch probabilities...
    ## Markov chain construction...
    ## Identification of terminal states...
    ## Computing fundamental matrix and absorption probabilities...
    ## Project results to all cells...

``` r
# Store pseudotime results in meta.data for easy plotting
ps <- mye_small@misc$Palantir$Pseudotime
colnames(ps)[3:4] <- c("fate1", "fate2")
mye_small@meta.data[,colnames(ps)] <- ps

# Visualize pseudotime and cell fates
DimPlot2(
  mye_small,
  features = colnames(ps),
  reduction = "ms",
  cols = list(continuous = "A", Entropy = "D"),
  theme = NoAxes())
```

![](vignettes/quick_start_files/figure-gfm/unnamed-chunk-14-1.png)<!-- -->

#### Visualization Along Trajectories

Visualizing gene expression or regulon activity along calculated
trajectories can provide insights into dynamic changes:

``` r
# Create smoothed gene expression curves along trajectory
GeneTrendCurve.Palantir(
  mye_small,
  pseudotime.data = ps,
  features = c("CD14", "FCGR3A")
)
```

![](vignettes/quick_start_files/figure-gfm/unnamed-chunk-15-1.png)<!-- -->

``` r
# Create a gene trend heatmap for different fates
GeneTrendHeatmap.Palantir(
  mye_small,
  features = VariableFeatures(mye_small)[1:10],
  pseudotime.data = ps,
  lineage = "fate1"
)
```

![](vignettes/quick_start_files/figure-gfm/unnamed-chunk-16-1.png)<!-- -->

### scVelo Analysis

scVelo is a Python tool used for RNA velocity analysis. We demonstrate
how to integrate and analyze velocyto-generated data within the Seurat
workflow using scVelo.

#### Preparing for scVelo

First, download the pre-calculated velocyto loom file:

``` r
# Download velocyto loom file
loom_path <- file.path(tempdir(), "pbmc10k_mye_small.loom")
download.file("https://zenodo.org/records/10944066/files/pbmc10k_mye_small.loom", 
              loom_path,
              mode = "wb")  # Use binary mode for Windows compatibility

# Set up the path for saving the AnnData object in the HDF5 (h5ad) format
if (.Platform$OS.type == "windows") {
    adata_path <- normalizePath(file.path(tempdir(), "mye_small.h5ad"), winslash = "/")
} else {
    adata_path <- file.path(tempdir(), "mye_small.h5ad")
}

# Integrate Seurat Object and velocyto loom into an AnnData object
scVelo.SeuratToAnndata(
  mye_small,
  filename = adata_path,
  velocyto.loompath = loom_path,
  prefix = "sample1_",
  postfix = "-1"
)
```

    ## scVelo version: 0.3.0
    ## Filtered out 10891 genes that are detected 20 counts (shared).
    ## Normalized count data: X, spliced, unspliced.
    ## Extracted 2000 highly variable genes.
    ## Logarithmized X.
    ## computing neighbors
    ##     finished (0:00:00) --> added 
    ##     'distances' and 'connectivities', weighted adjacency matrices (adata.obsp)
    ## computing moments based on connectivities
    ##     finished (0:00:00) --> added 
    ##     'Ms' and 'Mu', moments of un/spliced abundances (adata.layers)
    ## computing velocities
    ##     finished (0:00:00) --> added 
    ##     'velocity', velocity vectors for each individual cell (adata.layers)
    ## computing velocity graph (using 1/28 cores)
    ## WARNING: Unable to create progress bar. Consider installing `tqdm` as `pip install tqdm` and `ipywidgets` as `pip install ipywidgets`,
    ## or disable the progress bar using `show_progress_bar=False`.
    ##     finished (0:00:01) --> added 
    ##     'velocity_graph', sparse matrix with cosine correlations (adata.uns)

    ## NULL

#### Plotting scVelo Results

Once the data is processed, visualize the RNA velocity:

``` r
# Plot RNA velocity
scVelo.Plot(color = "cluster", basis = "ms_cell_embeddings", 
            save = "quick_start_scvelo.png", figsize = c(5,4))
```

<img src="vignettes/figures/scvelo_quick_start_scvelo.png" width="700" />

For detailed usage of the functions and more advanced analysis, please refer to the vignettes and tutorials.

## Single-Cell RNA-seq Analysis Course (New in v1.1.0)

A comprehensive 6-lesson course originally presented at the [Institute for AI in Medicine (IKIM)](https://www.ikim.uk-essen.de/institute), University Hospital Essen on October 8, 2024, organized by the [Department of Applied Computational Cancer Research](https://www.ikim.uk-essen.de/groups/accr). The course materials have been updated for SeuratExtend v1.1.0 and are now freely available online. Starting with fundamentals of R and Seurat, the course progressively builds to cover enhanced visualization, functional analysis, quality control, and cutting-edge methods including trajectory analysis, regulatory networks, and cell-cell communication. Perfect for beginners while providing depth needed for advanced applications.

### [Lesson 1: Introduction to R Programming](https://huayc09.github.io/SeuratExtend/articles/single-cell-course/1.basic-R.html)
Essential R programming fundamentals tailored for scRNA-seq analysis. Covers basic data types, data structures (vectors, matrices, data frames), file operations, and package management. Perfect for beginners starting their journey in bioinformatics.

### [Lesson 2: Basic Single-Cell Analysis with Seurat](https://huayc09.github.io/SeuratExtend/articles/single-cell-course/2.Seurat.html)
Comprehensive walkthrough of the standard Seurat workflow, from raw count matrix to cell type annotation. Learn about data normalization, dimensionality reduction, clustering, and visualization through hands-on analysis of PBMC data.

### [Lesson 3: Advanced Visualization with SeuratExtend](https://huayc09.github.io/SeuratExtend/articles/single-cell-course/3.Visualization.html)
Master advanced visualization techniques using SeuratExtend's enhanced plotting functions. Explore DimPlot2, FeaturePlot3, Heatmap, and other tools to create publication-ready figures. Includes practical examples of customizing plots and color schemes.

### [Lesson 4: Gene Set Enrichment Analysis and Utilities](https://huayc09.github.io/SeuratExtend/articles/single-cell-course/4.GSEA.html)
Master functional enrichment analysis using GO and Reactome databases through SeuratExtend's integrated GSEA pipeline. Learn to perform custom gene set analysis, interpret enrichment scores, and utilize helpful utility functions for gene naming conversions and cell proportions.

### [Lesson 5: Core Workflow Enhancements](https://huayc09.github.io/SeuratExtend/articles/single-cell-course/5.Core-Enhancement.html)
Elevate your scRNA-seq analysis with advanced quality control, doublet removal, data integration using Harmony, cell cycle analysis, and alternative normalization methods like SCTransform. Understand key considerations for processing and analyzing multi-sample datasets.

### [Lesson 6: Advanced Analytical Methods (Part 1)](https://huayc09.github.io/SeuratExtend/articles/single-cell-course/6.Advanced.html) [(Part 2)](https://huayc09.github.io/SeuratExtend/articles/single-cell-course/6.Advanced-2.html)
Explore cutting-edge techniques including trajectory analysis with scVelo/Palantir, cell-cell communication using CellChat/NicheNet, regulatory network inference with SCENIC, and specialized analyses for TCR/BCR data and copy number variations.

## License

The SeuratExtend R package code is licensed under GPL-3.0.

The data files (*.rda files in the 'data' folder) are released under CC0 1.0 Universal (CC0 1.0) Public Domain Dedication, meaning they are in the public domain and can be used without any restrictions.

## Publications Using SeuratExtend

1. Hua, Y., Vella, G., Rambow, F., et al. (2022). Cancer immunotherapies transition endothelial cells into HEVs that generate TCF1+ T lymphocyte niches through a feed-forward loop. **Cancer Cell** 40, 1600-1618. https://doi.org/10.1016/j.ccell.2022.11.002
2. Hua, Y., Wu, N., Miao, J., Shen, M. (2023). Single-cell transcriptomic analysis in two patients with rare systemic autoinflammatory diseases treated with anti-TNF therapy. **Front. Immunol.** 14. https://doi.org/10.3389/fimmu.2023.1091336
3. Verhoeven, J., Jacobs, K.A., Rizzollo, F., Lodi, F., Hua, Y., Poźniak, J., Narayanan Srinivasan, A., Houbaert, D., Shankar, G., More, S., et al. (2023). Tumor endothelial cell autophagy is a key vascular-immune checkpoint in melanoma. **EMBO Mol. Med.** 15, e18028. https://doi.org/10.15252/emmm.202318028
4. Dobersalske, C., Rauschenbach, L., Hua, Y., Berliner, C., Steinbach, A., Grüneboom, A., Kokkaliaris, K.D., Heiland, D.H., Berger, P., Langer, S., et al. (2024). Cranioencephalic functional lymphoid units in glioblastoma. **Nat. Med.** https://doi.org/10.1038/s41591-024-03152-x
