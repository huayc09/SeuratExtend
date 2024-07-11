---
title: "SCENIC for Gene Regulatory Networks Analysis"
author: "Yichao Hua"
date: "2024-5-1"
version: "SeuratExtend v1.0.0"
---

## Table of Contents

1.  [Importing SCENIC Loom Files into Seurat](#importing-scenic-loom-files-into-seurat)
2.  [Visualizing SCENIC Results](#visualizing-scenic-results)

## Importing SCENIC Loom Files into Seurat

[SCENIC](https://www.nature.com/articles/nmeth.4463) (Single-Cell
Regulatory Network Inference and Clustering) is a computational method
that provides deep insights into the regulatory networks governing gene
expression in single cells. It is highly recommended to use the Nextflow
pipeline to run SCENIC, which can be found
[here](https://github.com/aertslab/SCENICprotocol/). This process
requires a loom file as input, which can be generated directly using the
`Seu2Loom()` function. Currently, the `SeuratExtend` package does not
integrate the `RunScenic` functionality directly (requiring the use of
Nextflow command line), but if there is user demand, this could be
considered for future updates. The workflow results in a file named
“pyscenic\_integrated-output.loom,” which includes a list of
transcription factors (TFs) and their regulated genes, as well as a
TF-cell matrix representing the AUCell values. The AUCell score
represents the enrichment score of all genes regulated by a TF,
indicating regulon activity. The `ImportPyscenicLoom()` function allows
for the import of SCENIC-generated loom files into Seurat objects,
facilitating further analysis and visualization within the Seurat
framework.

As an example, we use a pre-computed SCENIC loom file which can be
downloaded as follows:

```{r}
library(SeuratExtend)
scenic_loom_path <- file.path(tempdir(), "pyscenic_integrated-output.loom")
download.file("https://zenodo.org/records/10944066/files/pbmc3k_small_pyscenic_integrated-output.loom", scenic_loom_path)

# Importing SCENIC Loom Files into Seurat
pbmc <- ImportPyscenicLoom(scenic_loom_path, seu = pbmc)
```

If you prefer to import SCENIC results without specifying a Seurat
object:

```{r}
# Importing SCENIC results without an existing Seurat object
scenic_output <- ImportPyscenicLoom(scenic_loom_path)
```

### Examining SCENIC Outputs

SCENIC results are stored in `seu@misc$SCENIC`, which includes
`seu@misc$SCENIC$Regulons`, a list where each element’s name is a TF
name and the value is a list of genes it regulates. Additionally,
`seu@misc$SCENIC$RegulonsAUC`, a TF-cell AUCell matrix, is also loaded
into the “TF” assay of the Seurat object, making the manipulation of TF
regulon activity as straightforward as handling gene expression data.

```{r}
# Viewing the outputs
tf_auc <- pbmc@misc$SCENIC$RegulonsAUC
head(tf_auc, 4:5)
```

    ##                         AHR     ARID3A        ARNT     ARNTL        ATF1
    ## CTATAAGATCGTTT-1 0.01406902 0.03347861 0.000000000 0.1144568 0.016730526
    ## GTGATTCTGGTTCA-1 0.00000000 0.00000000 0.000000000 0.1730939 0.017191920
    ## ACGTTGGACCGTAA-1 0.00000000 0.03189382 0.000000000 0.1463286 0.003671087
    ## GGATACTGCAGCTA-1 0.02483910 0.01347068 0.006867406 0.0945589 0.016670344

```{r}
tf_gene_list <- pbmc@misc$SCENIC$Regulons
head(tf_gene_list, 5)
```

    ## $AHR
    ##  [1] "NFYC-AS1" "MYSM1"    "ZZZ3"     "FUBP1"    "WDR77"    "TMEM183A" "TGOLN2"   "SEC22C"   "ATXN7"   
    ## [10] "RBPJ"     "C5orf24"  "CREBRF"   "POLR1C"   "YIPF3"    "RUNX2"    "BCLAF1"   "HOXA9"    "HOXA10"  
    ## [19] "CREB5"    "STAG3"   
    ##  [ reached getOption("max.print") -- omitted 29 entries ]
    ## 
    ## $ARID3A
    ##  [1] "CDC7"      "IL6R"      "LAMC1"     "LINC01136" "CDC42EP3"  "HNMT"      "GPBAR1"    "OARD1"    
    ##  [9] "GNB2"      "ATP6V1F"   "YWHAZ"     "XPA"       "ANKRD22"   "APLP2"     "PSMB5"     "LGALS3"   
    ## [17] "RPH3AL"    "EMILIN2"   "ME2"       "MAFB"     
    ##  [ reached getOption("max.print") -- omitted 4 entries ]
    ## 
    ## $ARNT
    ## [1] "EPHB3"   "SMIM14"  "GPR68"   "ANXA2"   "MORF4L1" "ZFP3"   
    ## 
    ## $ARNTL
    ## [1] "SLC4A10" "AP2M1"   "EMC2"    "PHF20L1" "SLC43A1" "CFL1"    "SPG21"   "CRB3"    "RPS28"  
    ## 
    ## $ATF1
    ##  [1] "CAMK2N1"  "PAFAH2"   "BCAS2"    "S100A13"  "ZNF281"   "LIN9"     "THUMPD2"  "CALM2"    "ACVR1"   
    ## [10] "NAB1"     "ZDBF2"    "SNRK"     "IFRD2"    "DUSP7"    "MITF"     "FOXP1"    "DHX36"    "KIAA0232"
    ## [19] "TEC"      "G3BP2"   
    ##  [ reached getOption("max.print") -- omitted 59 entries ]

## Visualizing SCENIC Results

Once SCENIC data is integrated into a Seurat object, users can leverage
a variety of visualization tools provided in the [Enhanced
Visualization](#enhanced-visualization) section to explore and interpret
these regulatory networks. Both the extracted `tf_auc` matrix or the
Seurat object itself can be used as inputs. Here are some practical
examples:

### Identifying Top Activated TFs in Each Cluster

```{r}
tf_zscore <- CalcStats(tf_auc, f = pbmc$cluster, order = "p", n = 4, t = TRUE)
Heatmap(tf_zscore, lab_fill = "zscore")
```

![](SCENIC_files/figure-markdown_strict/unnamed-chunk-4-1.png)

### Comparing TF Gene Expression Levels and Regulon Activity (AUCell)

Since we have imported SCENIC results into the “TF” assay, we can easily
access the corresponding AUCell values by prefixing “tf\_” to the TF
name:

```{r}
DimPlot2(
  pbmc,
  features = c("ETS1", "ATF3", "tf_ETS1", "tf_ATF3"),
  cols = list("tf_ETS1" = "D", "tf_ATF3" = "D")
)
```

![](SCENIC_files/figure-markdown_strict/unnamed-chunk-5-1.png)

### Simplifying Regulon Activity Access by Setting Default Assay

If you find manually adding “tf\_” to each transcription factor
cumbersome, you can set the default assay to “TF”, which simplifies
operations involving regulon activity. For example, to create a
waterfall plot that compares the regulon activity between two cell
types, you can do the following:

```{r}
# Setting the default assay to "TF" for easier access to regulon activity
DefaultAssay(pbmc) <- "TF"

# Creating a waterfall plot to compare regulon activity between monocytes and CD8 T cells
WaterfallPlot(
  pbmc,
  features = rownames(pbmc),  # Using all available TFs in the "TF" assay
  ident.1 = "Mono CD14",      # First group of cells
  ident.2 = "CD8 T cell",     # Second group of cells
  exp.transform = FALSE,      # Disable transformation of expression data
  top.n = 20                  # Display the top 20 most differentially active TFs
)
```

![](SCENIC_files/figure-markdown_strict/unnamed-chunk-6-1.png)

These examples illustrate how to integrate and utilize SCENIC analysis
within the Seurat framework, providing a comprehensive approach to
understanding gene regulatory mechanisms at the single-cell level.
