# load libraries
library(edgeR)

# load master matrix
colon_raw_counts <- readRDS("C:/Users/khali/CRISPR_thesis/colon_raw_counts.rds")
cat(sprintf("Loaded: %d sgRNAs x %d cols\n", nrow(colon_raw_counts), ncol(colon_raw_counts)))

# separate annotation and count columns
annotation_cols <- c("Construct_Barcode", "gene", "is_essential", "is_nonessential")
count_cols <- colnames(colon_raw_counts)[!colnames(colon_raw_counts) %in% annotation_cols]

# separate pDNA and cell line columns
pdna_cols <- count_cols[grepl("pdna", count_cols, ignore.case = TRUE)]
cellline_cols <- count_cols[!grepl("pdna", count_cols, ignore.case = TRUE)]

cat(sprintf("pDNA columns: %d\n", length(pdna_cols)))
cat(sprintf("Cell line columns: %d\n", length(cellline_cols)))

# extract count matrix only (no annotation columns)
count_matrix <- as.matrix(colon_raw_counts[, count_cols])
rownames(count_matrix) <- make.unique(colon_raw_counts$Construct_Barcode)

# TMM normalisation using edgeR
dge <- DGEList(counts = count_matrix)
dge <- calcNormFactors(dge, method = "TMM")
norm_counts <- cpm(dge, normalized.lib.sizes = TRUE, log = FALSE)

cat(sprintf("Normalised matrix: %d rows x %d cols\n", nrow(norm_counts), ncol(norm_counts)))

# check medians across samples
cat("\nMedian counts per slfc ample after TMM normalisation:\n")
print(round(apply(norm_counts, 2, median), 2))

# define which pDNA columns belong to which batch
pdna_batch2 <- pdna_cols[grepl("batch2", pdna_cols)]
pdna_batch3 <- pdna_cols[grepl("batch3", pdna_cols)]

cat(sprintf("pDNA batch2 columns: %d\n", length(pdna_batch2)))
cat(sprintf("pDNA batch3 columns: %d\n", length(pdna_batch3)))

# calculate mean pDNA per batch as T0 reference
pdna_batch2_mean <- rowMeans(norm_counts[, pdna_batch2] + 1)
pdna_batch3_mean <- rowMeans(norm_counts[, pdna_batch3] + 1)

# calculate LFC per cell line vs its correct pDNA batch
lfc_matrix <- matrix(NA, nrow = nrow(norm_counts), ncol = length(cellline_cols))
rownames(lfc_matrix) <- rownames(norm_counts)
colnames(lfc_matrix) <- cellline_cols

for (col in cellline_cols) {
  if (grepl("batch2", col)) {
    lfc_matrix[, col] <- log2((norm_counts[, col] + 1) / pdna_batch2_mean)
  } else {
    lfc_matrix[, col] <- log2((norm_counts[, col] + 1) / pdna_batch3_mean)
  }
}

cat(sprintf("LFC matrix: %d sgRNAs x %d cell lines\n", 
            nrow(lfc_matrix), ncol(lfc_matrix)))

# quick validation - essential genes should have negative LFC
essential_guides <- colon_raw_counts$is_essential
nonessential_guides <- colon_raw_counts$is_nonessential

cat(sprintf("\nMean LFC of essential guides:    %.3f\n", 
            mean(lfc_matrix[essential_guides, ])))
cat(sprintf("Mean LFC of nonessential guides: %.3f\n", 
            mean(lfc_matrix[nonessential_guides, ])))

# aggregate guide LFCs to gene level using median
gene_names <- colon_raw_counts$gene

# get unique genes
unique_genes <- unique(gene_names)
cat(sprintf("Unique genes: %d\n", length(unique_genes)))

# create gene level LFC matrix
gene_lfc_matrix <- matrix(NA, nrow = length(unique_genes), ncol = length(cellline_cols))
rownames(gene_lfc_matrix) <- unique_genes
colnames(gene_lfc_matrix) <- cellline_cols

for (gene in unique_genes) {
  gene_guides <- which(gene_names == gene)
  if (length(gene_guides) == 1) {
    gene_lfc_matrix[gene, ] <- lfc_matrix[gene_guides, ]
  } else {
    gene_lfc_matrix[gene, ] <- apply(lfc_matrix[gene_guides, ], 2, median)
  }
}

cat(sprintf("Gene level LFC matrix: %d genes x %d cell lines\n",
            nrow(gene_lfc_matrix), ncol(gene_lfc_matrix)))

# check top depleted genes across all cell lines
mean_lfc <- rowMeans(gene_lfc_matrix)
top20 <- sort(mean_lfc)[1:20]
cat("\nTop 20 most depleted genes across all colorectal cell lines:\n")
print(round(top20, 3))

# edgeR statistical testing per cell line
results_list <- list()

for (cell_line in unique(colon_replicates_metadata$DepMap_ID)) {
  
  # get cell line name
  cell_name <- cell_line
  
  # get replicates for this cell line
  cell_reps <- colon_replicates_metadata$replicate_ID[colon_replicates_metadata$DepMap_ID == cell_line]
  cell_reps <- intersect(cell_reps, cellline_cols)
  
  # get pDNA batch for this cell line
  batch <- unique(colon_replicates_metadata$pDNA_batch[colon_replicates_metadata$DepMap_ID == cell_line])
  pdna_reps <- pdna_cols[grepl(paste0("batch", batch), pdna_cols)]
  
  # build count matrix for this cell line
  cell_counts <- count_matrix[, c(pdna_reps, cell_reps)]
  
  # build group factor
  group <- factor(c(rep("pDNA", length(pdna_reps)), rep("cellline", length(cell_reps))))
  
  # create DGEList
  dge_cell <- DGEList(counts = cell_counts, group = group, genes = colon_raw_counts$gene)
  dge_cell <- calcNormFactors(dge_cell, method = "TMM")
  
  # design matrix
  design <- model.matrix(~group)
  
  # estimate dispersion
  dge_cell <- estimateDisp(dge_cell,design)
  
  # fit model
  fit <- glmFit(dge_cell, design)
  lrt <- glmLRT(fit)
  
  # extract results
  res <- topTags(lrt, n = Inf, sort.by = "none")$table
  res$cell_line <- cell_name
  res$DepMap_ID <- cell_line
  
  results_list[[cell_line]] <- res
  cat(sprintf("Processed: %s\n", cell_name))
}

# combine all cell line results into one dataframe
all_results <- do.call(rbind, results_list)
cat(sprintf("Total results: %d rows x %d cols\n", nrow(all_results), ncol(all_results)))

# filter essential genes FDR < 0.05 and LFC < -1
essential_hits <- all_results[all_results$FDR < 0.05 & all_results$logFC < -1, ]
cat(sprintf("Essential gene hits: %d\n", nrow(essential_hits)))

# check how many unique genes and cell lines
cat(sprintf("Unique essential genes: %d\n", length(unique(essential_hits$genes))))
cat(sprintf("Cell lines with hits: %d\n", length(unique(essential_hits$cell_line))))

# find common essentials - essential in >= 50% of cell lines (>= 5 out of 9)
n_cell_lines <- length(unique(colon_replicates_metadata$DepMap_ID))
threshold <- ceiling(n_cell_lines * 0.5)
gene_counts <- table(essential_hits$genes)
common_essentials <- names(gene_counts[gene_counts >= threshold])

cat(sprintf("\nCell lines tested: %d\n", n_cell_lines))
cat(sprintf("Threshold: essential in >= %d cell lines\n", threshold))
cat(sprintf("Common essential genes: %d\n", length(common_essentials)))

# save common essentials as dataframe
common_essentials_df <- data.frame(
  gene = common_essentials,
  n_cell_lines_essential = as.integer(gene_counts[common_essentials])
)
common_essentials_df <- common_essentials_df[order(-common_essentials_df$n_cell_lines_essential), ]

cat("\nTop 20 common essential genes:\n")
print(head(common_essentials_df, 20))

# save results
write.csv(all_results, 
          "C:/Users/khali/CRISPR_thesis/edgeR_all_results.csv",
          row.names = FALSE)
write.csv(common_essentials_df,
          "C:/Users/khali/CRISPR_thesis/common_essentials.csv",
          row.names = FALSE)
cat("\nResults saved!\n")

