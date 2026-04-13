# load the libraries
library(readr)
library(dplyr)
library(stringr)
library(data.table)

# load sample info
file_path_name <- "C:/Users/khali/Downloads/12280541/sample_info.csv"
full_sample_info <- read_csv(file_path_name, show_col_types = FALSE)

# extract relevant columns
columns <- c("DepMap_ID", "sample_collection_site", "primary_or_metastasis",
             "primary_disease", "Subtype", "sex", "age")
sample_info <- full_sample_info[, columns]

# filter colon colorectal cell lines
colon_lines <- sample_info %>%
  filter(
    str_detect(sample_collection_site, regex("colon", ignore_case = TRUE)) &
      str_detect(primary_disease, regex("colorectal|colon|rectal", ignore_case = TRUE))
  )
cat(sprintf("Colorectal cell lines found: %d\n", nrow(colon_lines)))

# load replicate map
file_path_name_rep <- "C:/Users/khali/Downloads/12280541/Achilles_replicate_map.csv"
full_sample_info_rep <- read_csv(file_path_name_rep, show_col_types = FALSE)

# match colon DepMap IDs to replicate IDs
colon_replicates <- full_sample_info_rep %>%
  filter(DepMap_ID %in% colon_lines$DepMap_ID)

# keep only QC passing replicates
colon_replicates_qc <- colon_replicates %>%
  filter(passes_QC == TRUE)

cat(sprintf("Replicates before QC filter: %d\n", nrow(colon_replicates)))
cat(sprintf("Replicates after QC filter:  %d\n", nrow(colon_replicates_qc)))

# merge cell line metadata onto replicates
colon_replicates_metadata <- colon_replicates_qc %>%
  left_join(colon_lines, by = "DepMap_ID")
cat(sprintf("Replicate metadata: %d rows x %d cols\n",
            nrow(colon_replicates_metadata), ncol(colon_replicates_metadata)))

# load raw counts
file_path_raw <- "C:/Users/khali/Downloads/12280541/Achilles_raw_readcounts.csv"
raw_header <- fread(file_path_raw, nrows = 0)

# get matched replicate columns
colon_replicate_ids <- colon_replicates_qc$replicate_ID
matched_cols <- intersect(colon_replicate_ids, colnames(raw_header))
cat(sprintf("Replicate columns matched: %d\n", length(matched_cols)))

# get pDNA columns
all_cols <- colnames(raw_header)
pdna_cols <- all_cols[str_detect(all_cols, regex("pdna", ignore_case = TRUE))]
cat(sprintf("pDNA columns found: %d\n", length(pdna_cols)))

# combine sgRNA + pDNA + colon replicate columns and load
sgrna_col <- all_cols[1]
keep_cols <- c(sgrna_col, pdna_cols, matched_cols)
colon_raw_counts <- as.data.frame(fread(file_path_raw, select = keep_cols))

# rename sgRNA column
colnames(colon_raw_counts)[1] <- "Construct_Barcode"
cat(sprintf("Raw counts loaded: %d sgRNAs x %d cols\n",
            nrow(colon_raw_counts), ncol(colon_raw_counts)))

# load guide map
file_path_guide <- "C:/Users/khali/Downloads/12280541/Achilles_guide_map.csv"
guide_map <- read_csv(file_path_guide, show_col_types = FALSE)
cat(sprintf("Guide map loaded: %d rows x %d cols\n", nrow(guide_map), ncol(guide_map)))

# load common essentials
Achilles_common_path <- "C:/Users/khali/Downloads/12280541/common_essentials.csv"
Achilles_common <- read_csv(Achilles_common_path, show_col_types = FALSE)
cat(sprintf("Common essential genes: %d\n", nrow(Achilles_common)))

# load nonessentials
Achilles_noncommon_path <- "C:/Users/khali/Downloads/12280541/nonessentials.csv"
Achilles_noncommon <- read_csv(Achilles_noncommon_path, show_col_types = FALSE)
cat(sprintf("Nonessential genes: %d\n", nrow(Achilles_noncommon)))

# flag essential and nonessential genes in guide map
guide_map <- guide_map %>%
  mutate(
    is_essential    = gene %in% Achilles_common$gene,
    is_nonessential = gene %in% Achilles_noncommon$gene
  )
cat(sprintf("Essential guides flagged:    %d\n", sum(guide_map$is_essential)))
cat(sprintf("Nonessential guides flagged: %d\n", sum(guide_map$is_nonessential)))
cat(sprintf("Test genes:                  %d\n",
            sum(!guide_map$is_essential & !guide_map$is_nonessential)))

# merge guide map onto master counts
colon_raw_counts <- merge(colon_raw_counts,
                          guide_map[, c("sgrna", "gene", "is_essential", "is_nonessential")],
                          by.x = "Construct_Barcode",
                          by.y = "sgrna",
                          all.x = TRUE)

# move gene column to front
colon_raw_counts <- colon_raw_counts %>%
  relocate(gene, is_essential, is_nonessential, .after = 1)

cat(sprintf("Master matrix after merge: %d rows x %d cols\n",
            nrow(colon_raw_counts), ncol(colon_raw_counts)))

# load dropped guides
dropped_guides_path <- "C:/Users/khali/Downloads/12280541/Achilles_dropped_guides.csv"
dropped_guides <- read_csv(dropped_guides_path, show_col_types = FALSE)
colnames(dropped_guides)[1] <- "Construct_Barcode"
cat(sprintf("Dropped guides loaded: %d\n", nrow(dropped_guides)))

# remove dropped guides from master matrix
before <- nrow(colon_raw_counts)
colon_raw_counts <- colon_raw_counts %>%
  filter(!Construct_Barcode %in% dropped_guides$Construct_Barcode)
after <- nrow(colon_raw_counts)

cat(sprintf("Guides before removal: %d\n", before))
cat(sprintf("Guides removed:        %d\n", before - after))
cat(sprintf("Guides remaining:      %d\n", after))
cat(sprintf("\nFinal master matrix: %d rows x %d cols\n",
            nrow(colon_raw_counts), ncol(colon_raw_counts)))