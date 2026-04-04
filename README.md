# CRISPR Screen Analysis Pipeline

## Overview
This pipeline analyses  CRISPR knockout screen data from the Achilles dataset.
It performs data loading, TMM normalisation and statistical testing to identify
essential genes across multiple cancer cell lines.

## Pipeline steps
- `01_load_data.R` — loads all cell line count files
- `02_normalisation.R` — TMM normalisation via edgeR
- `03_fold_change.R` — statistical testing via glmQLFTest / glmLRT

## Requirements

### R version
R 4.1.0 or higher

### Packages
Install the following before running:

install.packages("tidyverse")
install.packages("BiocManager")
BiocManager::install("edgeR")

## How to run

### Step 1: Clone the repository
git clone https://github.com/kell0309/CRISPR_thesis.git
cd CRISPR_thesis

### Step 2: Open RStudio
Open RStudio and set your working directory to the cloned folder:
setwd("path/to/CRISPR_thesis")

### Step 3:  Run the scripts in order

#### Run script 1: loading
("Achilles_scripts/01_load_data.R")
 You will be prompted to enter the path to your data folder
 Example: C:/Users/username/Documents/data/achilles

#### Run script 2:  normalisation
("Achilles_scripts/02_normalization.R")

#### Run script 3: statistical testing
("Achilles_scripts/03_fold_change.R")

## Data format
Each input file should be a tab separated .txt file.

## Notes
- Cell lines with a single replicate are analysed using a fixed biological
  coefficient of variation of 0.4
- All other cell lines use dispersion estimated directly from the replicates
- FDR threshold of 0.05 is recommended for significant hits
