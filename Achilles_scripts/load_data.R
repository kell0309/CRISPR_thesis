# 01_load_data.R
# Loading libraries
library(tidyverse)
library(data.table)

# Ask user for folder path
file_path <- readline(prompt = "Enter the full path to the file: ")

# Find textfiles
txt_files <- list.files(file_path, pattern = "\\.txt$", full.names = TRUE)

if (length(txt_files) == 0) {
  stop("No .txt files found in the folder.")
}

cat("Found", length(txt_files), "files:\n")
print(basename(txt_files))

# Read all files into a named list
counts_list <- txt_files |>                                          
  set_names(tools::file_path_sans_ext(basename(txt_files))) |>
  map(~ read_tsv(.x, show_col_types = FALSE))

# Add calculated mean of reps per sgRNA
counts_list <- counts_list |>                                        
  map(function(df) {
    
    # identify rep columns (e.g. C2BBE1_Rep1, Rep2, Rep3)
    rep_cols <- df |>
      select(contains("_Rep")) |>
      colnames()
    
    # calculate mean across reps for each sgRNA row
    df <- df |>                                                      
      mutate(rep_mean_calculated = rowMeans(across(all_of(rep_cols)),
                                            na.rm = TRUE))
    return(df)
  })

# Sanity check
cat("\nDone! Files loaded:", length(counts_list), "\n")
cat("Example - first file rep mean (first 5 rows):\n")
print(counts_list[[1]] |> select(sgRNA, gene, rep_mean_calculated) |> head(5))