# 01_load_data.R
# Loading libraries
library(data.table)
library(tidyverse)
# Ask user for folder path
file_path <- readline(prompt = "Enter the full path to the file: ")

# Find textfiles
txt_files <- list.files(file_path, pattern = "\\.txt$", full.names = TRUE)

if (length(txt_files) == 0) {
  stop("No .txt files found in the folder.")
}

cat("Found", length(txt_files), "files:\n")
print(basename(txt_files))

# Get files without extensions
file_names <- basename(txt_files)
file_names <- gsub("\\.txt$","",file_names)
counts_list <- lapply(txt_files, function(x) read_tsv(x, show_col_types = FALSE))

# Name each element in the list by cell line
names(counts_list) <- file_mames

# identify rep columns for each file
counts_list <- lapply(counts_list, function(df) {
  
  rep_cols <- colnames(select(df, contains("_Rep")))
  
  return(df)
})

# Sanity check
cat("\nDone! Files loaded:", length(counts_list), "\n")
cat("first file rep mean (first 5 rows):\n")
print(head(counts_list[[1]]))
