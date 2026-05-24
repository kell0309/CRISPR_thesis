# load libraries
library(edgeR)
library(dplyr)

# load master matrix
colon_raw_counts <- readRDS("C:/Users/khali/CRISPR_thesis/results/colon_raw_counts.rds")
dim(colon_raw_counts)

# separate annotation and count columns
annotation_cols <- c("Construct_Barcode", "gene", "is_essential", "is_nonessential")
count_cols <- colnames(colon_raw_counts)[!colnames(colon_raw_counts) %in% annotation_cols]

# separate pDNA and cell line columns
pdna_cols <- count_cols[grepl("pdna", count_cols, ignore.case = TRUE)]
cellline_cols <- count_cols[!grepl("pdna", count_cols, ignore.case = TRUE)]

# view lenghts of both pdna and cellline 
length(pdna_cols)
length(cellline_cols)

# extract count matrix
count_matrix <- as.matrix(colon_raw_counts[, count_cols])
rownames(count_matrix) <- make.unique(colon_raw_counts$Construct_Barcode)

# TMM normalisation
dge <- DGEList(counts = count_matrix)
dge <- calcNormFactors(dge, method = "TMM")
norm_counts <- cpm(dge, normalized.lib.sizes = TRUE, log = FALSE)

#dimensions of the TMM normalisation
dim(norm_counts)

# checking the median counts per sample after the normalisation
print(round(apply(norm_counts, 2, median), 2))

# define pDNA batches and split them
pdna_batch2 <- pdna_cols[grepl("batch2", pdna_cols)]
pdna_batch3 <- pdna_cols[grepl("batch3", pdna_cols)]

length(pdna_batch2)
length(pdna_batch3)

# calculate mean pDNA per batch
pdna_batch2_mean <- rowMeans(norm_counts[, pdna_batch2] + 1)
pdna_batch3_mean <- rowMeans(norm_counts[, pdna_batch3] + 1)

# calculate LFC per cell line vs correct pDNA batch
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

# dimensionality of the logfold matrix
dim(lfc_matrix)

# validation 
essential_guides <- colon_raw_counts$is_essential
nonessential_guides <- colon_raw_counts$is_nonessential


# aggregate guide LFCs to gene level using median
gene_names <- colon_raw_counts$gene
unique_genes <- unique(gene_names)
length(unique_genes)

# creating lfc matrix
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

# dimensions check for the gene lfc matrix 
dim(gene_lfc_matrix)

# top 20 depleted genes
mean_lfc <- rowMeans(gene_lfc_matrix)
top20 <- sort(mean_lfc)[1:20]
cat("\nTop 20 most depleted genes across all colorectal cell lines:\n")
print(round(top20, 3))

# edgeR statistical testing per cell line
results_list <- list()

for (cell_line in unique(colon_replicates_metadata$DepMap_ID)) {
  
  # gett cell line and their replicate columns
  cell_name <- cell_line
  cell_reps <- colon_replicates_metadata$replicate_ID[colon_replicates_metadata$DepMap_ID == cell_line]
  
  # only keep the replicates that are present
  cell_reps <- intersect(cell_reps, cellline_cols)
  
  # assign the pdna batches belongs to which
  batch <- unique(colon_replicates_metadata$pDNA_batch[colon_replicates_metadata$DepMap_ID == cell_line])
  pdna_reps <- pdna_cols[grepl(paste0("batch", batch), pdna_cols)]
  
  # label the groups are pDNA and cell line
  cell_counts <- count_matrix[, c(pdna_reps, cell_reps)]
  group <- factor(c(rep("pDNA", length(pdna_reps)), rep("cellline", length(cell_reps))))
  
  # build the DGElist and normalise
  dge_cell <- DGEList(counts = cell_counts, group = group, genes = colon_raw_counts$gene)
  dge_cell <- calcNormFactors(dge_cell, method = "TMM")
  
  # build the deisgn matrix and also calculate the dispersion
  design <- model.matrix(~group)
  dge_cell <- estimateDisp(dge_cell, design)
  
  
  #fit the model and run the test
  fit <- glmFit(dge_cell, design)
  lrt <- glmLRT(fit)
  
  # extract all results 
  res <- topTags(lrt, n = Inf, sort.by = "none")$table
  res$cell_line <- cell_name
  res$DepMap_ID <- cell_line
  results_list[[cell_line]] <- res
}

# obtain dispersion and BCV
dge_cell$common.dispersion
sqrt(dge_cell$common.dispersion)

# combine results
all_results <- do.call(rbind, results_list)
dim(all_results)


# filter essential genes
essential_hits <- all_results[all_results$FDR < 0.05 & all_results$logFC < -1, ]
nrow(essential_hits)
length(unique(essential_hits$genes))
length(unique(essential_hits$cell_line))

# find common essentials - count unique cell lines per gene
n_cell_lines <- length(unique(colon_replicates_metadata$DepMap_ID))
threshold <- ceiling(n_cell_lines * 0.5)

gene_cellline_counts <- essential_hits %>%
  group_by(genes) %>%
  summarise(n_cell_lines = n_distinct(cell_line)) %>%
  as.data.frame()

common_essentials_df <- gene_cellline_counts[gene_cellline_counts$n_cell_lines >= threshold, ]
colnames(common_essentials_df) <- c("gene", "n_cell_lines_essential")
common_essentials_df <- common_essentials_df[order(-common_essentials_df$n_cell_lines_essential), ]

# total number of cell lines testes
n_cell_lines  

# minimum cell lines to qualify as common essential
threshold 
nrow(common_essentials_df)
max(common_essentials_df$n_cell_lines_essential)

# view the top 20 common essential gene
print(head(common_essentials_df, 20))

# save results
write.csv(all_results,
          "C:/Users/khali/CRISPR_thesis/results/edgeR_all_results.csv",
          row.names = FALSE)
write.csv(common_essentials_df,
          "C:/Users/khali/CRISPR_thesis/results/common_essentials.csv",
          row.names = FALSE)