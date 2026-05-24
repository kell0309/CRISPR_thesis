Sanger Pipeline for CRISPR Cancer Screening data essential genes

This is the Sanger processing showing the code used for the 4 steps: 
-import datasets
statistical tests
export_function
Processing datasets.
The following packages are required for the code to function as well as their dependencies: 
dplyr
pheatmap
ggVenDiagram
corrplots
ggplots2
ggrepel
as well as biomanager packages:edgerR, limma, clusterprofiler, org.Hs.eg.db, enrichplot.

Initially an additional folder with the unprocessed data was uploaded in this repository however due to size limitations it was removed after the results were acquired.
The original dataset comes from the Sanger CRISPR Dependencies 2024.
The end result of the exporting is vising under Sanger_Results, no plots were saved in the git repository.

All code is currated and reading through it is expected for proper use. Only the first three scripts were intended to be run all at once.

Import datasets: The importation and pre processing stage combines the raw sanger Count data into a single masterframe with merged replicates that will be processed. It is intended to be given an input file path.
The Sanger CRISPR screening raw counts are expected to be renamed after their cell line and replicate identifiers as present in their zip.
This code includes silence categorisation functions to use if needed.

statistical tests: The Statistical testing using edgeR functions normalise the datasets and creates a DGE object based on statistical tests. in addition it constructing pseudo-essential genes used in processing
The statistical tests include TMM normalisation, and statistical tests to measure for accuracy of the data. Lastly it also measures Logarithmic Fold Change as the main measure
This code include silenced functions for more complex datasets

export_function: The exporting was created to compare the results of the original run to future versions as well as for any contrast with Achilles. It is only meant to be used for external data

processing datasets: The processing attempts to isolate data from the main results as well as construct plots. It is advised to run in individual chunks and carefully curating the code.
This code includes pearson's correlation, Spearman's correlation, Cross dataset comparisons, Enrichment and Pathway analysis. As well as plots in order to represent said results.

Author: Vasileios Theocharis
For background information read our Bachelor's thesis project or contact us.
