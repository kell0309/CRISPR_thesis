# Achilles CRISPR Screen Analysis

## Overview
Analysis of the DepMap 20Q2 Achilles genome wide CRISPR knockout screen to identify essential genes in colorectal carcinoma cell lines. The pipeline processes raw sgRNA count data from 9 colorectal cancer cell lines across 18 replicates, applying TMM normalisation, batch-specific LFC calculation and edgeR statistical testing to identify 33 common essential genes. MAGeCK MLE was additionally applied as a complementary analytical approach.

Data was obtained from the Broad Institute Project Achilles https://figshare.com/articles/dataset/DepMap_20Q2_Public/12280541

## Contents
**Achilles_scripts/**
- `formatting_data.R`; data loading, QC filtering and replicate metadata preparation
- `Statistical_test_ver2.R`; TMM normalisation, LFC calculation and edgeR GLM statistical testing
- `mageck_format.R`; count matrix and design matrix preparation for MAGeCK MLE
- `plot.R` ;all visualisation plots including QC, volcano, heatmap, UMAP and Venn diagram

**results/plots/**
- All figures generated from the analysis pipeline

**results/**
- `mageck_counts.txt`; sgRNA count matrix formatted for MAGeCK MLE input
- `mageck_design.txt`;design matrix encoding batch structure for MAGeCK MLE

## How to Run
Download DepMap 20Q2 from https://figshare.com/articles/dataset/DepMap_20Q2_Public/12280541 then run scripts in order:

set the working directory;
```r
setwd("path/to/CRISPR_thesis")
```
and then the scripts

```r
source("Achilles_scripts/formatting_data.R")
source("Achilles_scripts/Statistical_test_ver2.R")
source("Achilles_scripts/mageck_format.R")
source("Achilles_scripts/plot.R")
```

## Requirements
R 4.1.0 or higher with the following packages:
`edgeR`, `ggplot2`, `ggrepel`, `pheatmap`, `umap`, `ggvenn`, `readr`, `dplyr`


```r
install.packages("tidyverse")
install.packages("BiocManager")
BiocManager::install("edgeR")
install.packages("ggplot2")
install.packages("ggrepel")
install.packages("pheatmap")
install.packages("umap")
install.packages("ggvenn")
install.packages("readr")
install.packages("dplyr")
```
MAGeCK MLE was run via Galaxy EU (https://usegalaxy.eu/) version 0.5.9.2.1
