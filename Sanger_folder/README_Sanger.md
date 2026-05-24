Sanger Pipeline for CRISPR Cancer Screening data essential genes

This is the Sanger processing showing the code used for the 4 steps: importing/preprocessing, statistical testing, exporting, and Processing via plots.
Biomanager and Cran packages both required for this code to function.

The original dataset comes from the Sanger Depmap Achilles dataset.
The end result of the exporting is vising under Sanger_Results, no plots were saved in the git repository.

All code is currated and reading through it is expected for proper use. Only the first three scripts were intended to be run blindly.

Initially an additional folder with the unprocessed dataset was uploaded in this repository however due to size limitations it was removed after the results were acquired

The importation and pre processing stage combines the raw sanger Count data into a single masterframe that will be processed
The Statistical testing using edgeR functions normalises the datasets and creates a DGE object in addition to constructing pseudo-essential genes
The exporting was created to compare the results of the original run to future versions as well as for any contrast with Achilles.
The processing attempts to isolate data from the main results as well as construct plots. It is advised to run in individual chunks

Author: Vasileios Theocharis
For background information read our Bachelor's thesis project or contact us.
