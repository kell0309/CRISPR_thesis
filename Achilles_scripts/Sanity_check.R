# presentation_summary.R
library(tidyverse)

# ---- 1. number of cell lines and files ----
cat("=== DATASET OVERVIEW ===\n")
cat("Number of cell lines:", length(counts_list), "\n")
cat("Cell lines:\n")
print(names(counts_list))

# ---- 2. gene and sgrna counts per cell line ----
cat("\n=== GENES AND sgRNAs PER CELL LINE ===\n")
for (name in names(counts_list)) {
  cat("\n", name, "\n")
  cat("  Total sgRNAs :", nrow(counts_list[[name]]), "\n")
  cat("  Unique genes :", n_distinct(counts_list[[name]]$gene), "\n")
}

# ---- 3. genes appearing in most cell lines ----
cat("\n=== TOP 20 GENES ACROSS MOST CELL LINES ===\n")

# get all genes from all files into one dataframe
all_genes <- data.frame()
for (name in names(counts_list)) {
  temp <- data.frame(
    cell_line = name,
    gene      = unique(counts_list[[name]]$gene)
  )
  all_genes <- rbind(all_genes, temp)
}

# count how many cell lines each gene appears in
gene_counts <- as.data.frame(table(all_genes$gene))
colnames(gene_counts) <- c("gene", "n_cell_lines")
gene_counts <- gene_counts[order(-gene_counts$n_cell_lines), ]

print(head(gene_counts, 20))

# ---- 4. save summary to csv ----
write.csv(gene_counts, "gene_frequency_across_celllines.csv", row.names = FALSE)
cat("\nSaved gene frequency table to gene_frequency_across_celllines.csv\n")

# ---- 5. plot top 20 most common genes ----
pdf("presentation_plots.pdf", width = 12, height = 8)

# plot 1 - top 20 genes across cell lines
top20 <- head(gene_counts, 20)
barplot(top20$n_cell_lines,
        names.arg = top20$gene,
        las       = 2,
        col       = "steelblue",
        main      = "Top 20 genes appearing across most cell lines",
        ylab      = "Number of cell lines",
        cex.names = 0.7)

# plot 2 - sgRNA count distribution before normalisation
n <- length(counts_list)
par(mfrow = c(ceiling(n/3), 3))
for (name in names(counts_list)) {
  hist(log2(counts_list[[name]]$rep_mean_calculated + 1),
       main   = name,
       xlab   = "log2(mean count + 1)",
       col    = "steelblue",
       breaks = 50)
}
par(mfrow = c(1,1))

# plot 3 - before vs after normalisation for first file
first <- names(counts_list)[1]
par(mfrow = c(1, 2))

hist(log2(counts_list[[first]]$rep_mean_calculated + 1),
     main   = paste(first, "- before"),
     xlab   = "log2(mean count + 1)",
     col    = "steelblue",
     breaks = 50)

rep_cols_norm <- colnames(select(counts_list_norm[[first]], contains("_Rep")))
norm_mean <- rowMeans(select(counts_list_norm[[first]], all_of(rep_cols_norm)))
hist(norm_mean,
     main   = paste(first, "- after"),
     xlab   = "log2 CPM",
     col    = "coral",
     breaks = 50)

par(mfrow = c(1,1))
dev.off()

cat("\nAll plots saved to presentation_plots.pdf\n")