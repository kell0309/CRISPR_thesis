
# Format Achilles Data into Sanger Format
# Output: one file per colorectal cell line
# Format: sgRNA | gene | rep1 | rep2 | ... | plasmid


library(readr)
library(readxl)
library(dplyr)


# STEP 1: Load reference files

guide_map <- read_csv("Achilles_guide_map.csv") %>%
  select(sgrna, gene) %>%
  mutate(gene = gsub(" \\(.*\\)", "", gene))  # Remove Entrez IDs

sample_info <- read_csv("sample_info.csv") %>%
  select(DepMap_ID, stripped_cell_line_name, lineage)

replicate_map <- read_csv("Achilles_replicate_map.csv") %>%
  filter(passes_QC == "True")


# STEP 2: Identify colorectal cell lines and their replicates


colorectal_samples <- sample_info %>%
  filter(lineage == "colorectal") %>%
  inner_join(replicate_map, by = "DepMap_ID")

cat("Colorectal cell lines:", length(unique(colorectal_samples$stripped_cell_line_name)), "\n")


# STEP 3: Load Achilles count matrix


cat("Loading Achilles count matrix (may take a moment)...\n")
achilles_raw <- read_excel("Achilles_colorectal_with_pDNA.xlsx")
colnames(achilles_raw)[1] <- "barcode"

all_cols <- colnames(achilles_raw)

# Define pDNA columns per batch
pdna_batch2 <- grep("batch2", grep("pDNA|pDNA", all_cols, value = TRUE), value = TRUE)
pdna_batch3 <- grep("batch3", grep("pDNA|pDNA", all_cols, value = TRUE), value = TRUE)
pdna_batch4 <- grep("batch4", grep("pDNA|pDNA", all_cols, value = TRUE), value = TRUE)

cat("pDNA batch2 columns:", length(pdna_batch2), "\n")
cat("pDNA batch3 columns:", length(pdna_batch3), "\n")
cat("pDNA batch4 columns:", length(pdna_batch4), "\n")

# Map barcodes to sgRNA names and genes
achilles_mapped <- achilles_raw %>%
  inner_join(guide_map, by = c("barcode" = "sgrna")) %>%
  relocate(barcode, gene)

cat("sgRNAs after gene mapping:", nrow(achilles_mapped), "\n")


# STEP 4: Create output directory


output_dir <- "achilles_per_cellline"
if (!dir.exists(output_dir)) dir.create(output_dir)


# STEP 5: Generate one file per cell line


cell_lines <- unique(colorectal_samples$stripped_cell_line_name)
cat("\nGenerating files for", length(cell_lines), "cell lines...\n")

for (cl in cell_lines) {
  
  # Get replicates and batch for this cell line
  cl_info <- colorectal_samples %>% filter(stripped_cell_line_name == cl)
  rep_ids  <- cl_info$replicate_ID
  batch    <- unique(cl_info$pDNA_batch)
  
  # Match replicate column names in the Achilles file
  rep_cols <- all_cols[all_cols %in% rep_ids]
  
  if (length(rep_cols) == 0) {
    cat("  Skipping", cl, "- no matching columns found\n")
    next
  }
  
  # Select appropriate pDNA columns based on batch
  pdna_cols <- switch(as.character(batch[1]),
                      "2" = pdna_batch2,
                      "3" = pdna_batch3,
                      "4" = pdna_batch4,
                      pdna_batch3  
  )
  
  if (length(pdna_cols) == 0) {
    cat("  Skipping", cl, "- no pDNA columns for batch", batch[1], "\n")
    next
  }
  
  # Build output: sgRNA | gene | replicates | plasmid
  cl_data <- achilles_mapped %>%
    select(barcode, gene, all_of(rep_cols), all_of(pdna_cols)) %>%
    mutate(plasmid = rowMeans(across(all_of(pdna_cols)))) %>%
    select(sgRNA = barcode, gene, all_of(rep_cols), plasmid)
  
  # Rename replicate columns to clean names
  rep_new_names <- paste0(cl, "_Rep", seq_along(rep_cols))
  colnames(cl_data)[3:(2 + length(rep_cols))] <- rep_new_names
  
  # Save file
  filename <- file.path(output_dir, paste0(cl, "_counts.txt"))
  write_tsv(cl_data, filename)
  
  cat("  Saved:", cl, "->", length(rep_cols), "replicates,", nrow(cl_data), "sgRNAs\n")
}

cat("\nDone! All files saved in:", output_dir, "\n")
cat("Each file has columns: sgRNA | gene | Rep1 | Rep2 | ... | plasmid\n")
cat("These are ready to combine with your Sanger file for analysis.\n")