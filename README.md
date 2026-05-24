# Achilles CRISPR Screen Analysis

**Author:** Khalif Hussein  
**Institution:** University of Skövde  
**Project:** Bachelor Degree Thesis in Bioinformatics, 2026



## Overview

This folder contains the bioinformatics pipeline developed for the analysis of the  Achilles genome-wide CRISPR knockout screen dataset. The analysis was performed as part of a comparative study of CRISPR screening analytical methods applied to colorectal carcinoma cell lines.

The raw sgRNA count data was obtained from the Broad Institute Project Achilles dataset.



## Contents

### Achilles_scripts/
The R scripts used for data processing, statistical analysis and visualisation:

- `formatting_data.R`
- `Statistical_test_ver2.R`
- `mageck_format.R`
- `plot.R`

### results/plots/
Figures generated from the analysis pipeline including quality control plots, volcano plot, heatmap, bar chart, UMAP and Venn diagram.

### results/
MAGeCK MLE input files used for the Galaxy EU analysis:

- `mageck_counts.txt`
- `mageck_design.txt`



## Data Source
**Dataset:** https://figshare.com/articles/dataset/DepMap_20Q2_Public/12280541
**Source:** Broad Institute Project Achilles  


## Dependencies
Analysis was performed in R using the following packages:
- `edgeR`
- `ggplot2`
- `ggrepel`
- `pheatmap`
- `umap`
- `ggvenn`
- `readr`
- `dplyr`

MAGeCK MLE analysis was performed using Galaxy EU (https://usegalaxy.eu/) version 0.5.9.2.1.
