# Load libraries
library(ggplot2)
library(ggrepel)
library(pheatmap)
library(umap)
library(readr)
library(ggvenn)


# source the scripts to load all variables
source("C:/Users/khali/CRISPR_thesis/Achilles_scripts/formatting_data.R")
source("C:/Users/khali/CRISPR_thesis/Achilles_scripts/Statistical_test_ver2.R")

 # output directory
output_dir <- "C:/Users/khali/CRISPR_thesis/results/plots"
dir.create(output_dir, showWarnings = FALSE)



# plot 1: LFC distribution of essential and nonessentials
# Calcuatinng the mean LFC across all the cell lines for each sgRNA guide group
lfc_values <- data.frame(
  LFC   = c(rowMeans(lfc_matrix[colon_raw_counts$is_essential, ]),
            rowMeans(lfc_matrix[colon_raw_counts$is_nonessential, ])),
  group = c(rep("Essential", sum(colon_raw_counts$is_essential)),
            rep("Nonessential", sum(colon_raw_counts$is_nonessential)))
)

# code for the plot density 
p1 <- ggplot(lfc_values, aes(x = LFC, fill = group)) +
  geom_density(alpha = 0.5) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "black") +
  labs(title = "Essential vs Nonessential Guide LFC Distribution",
       x = "Mean LFC", y = "Density", fill = "Guide Type") +
  theme_minimal()

# save the file in the directed output directory
ggsave(file.path(output_dir, "QC1_lfc_distribution.png"), p1, width = 8, height = 5)



# pLot 2: reads per depth boxpot
# Shape the count matrix to a long format for ggplot
counts_long <- data.frame(
  sample = rep(colnames(count_matrix), each = nrow(count_matrix)),
  counts = as.vector(count_matrix),
  type   = ifelse(grepl("pdna", rep(colnames(count_matrix), each = nrow(count_matrix)),
                        ignore.case = TRUE), "pDNA", "Cell Line")
)

# clean sample names for all cell lines remove the batch and protocol information for clean axis labelling
counts_long$sample_clean <- gsub("-311Cas9|_311Cas9|311Cas9|-Avana 4 |_Avana 4 ", "", counts_long$sample)
counts_long$sample_clean <- gsub("_p[0-9]_batch[0-9]|_batch[0-9]|Avana4pDNA20160601|Avana_4_Hu_pDNA_", "", counts_long$sample_clean)


# code for boxplot using ggplot
p2 <- ggplot(counts_long, aes(x = sample_clean, y = log2(counts + 1), fill = type)) +
  geom_boxplot(outlier.size = 0.3, outlier.alpha = 0.3) +
  scale_fill_manual(values = c("Cell Line" = "#2166ac", "pDNA" = "#d6604d")) +
  labs(title = "Guide Count Distribution per Sample",
       x = "", y = "log2(counts + 1)", fill = "Type") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, size = 10),
        legend.position = "top")

ggsave(file.path(output_dir, "QC2_reads_per_sample.png"), p2, width = 14, height = 6)



# plot 3: replicate s correlation heatmap
# pearson correlation on all cell line replicates
cor_matrix <- cor(norm_counts[, cellline_cols], method = "pearson")

# reshape to the long format for ggplot and title plot
cor_df <- as.data.frame(as.table(cor_matrix))
colnames(cor_df) <- c("Sample1", "Sample2", "Correlation")

# plot ggplot for correlation heatmap and set midpint to 0.95
p3 <- ggplot(cor_df, aes(x = Sample1, y = Sample2, fill = Correlation)) +
  geom_tile() +
  scale_fill_gradient2(low = "blue", mid = "white", high = "red", midpoint = 0.95) +
  labs(title = "Replicate Correlation Heatmap", x = "", y = "") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, size = 6),
        axis.text.y = element_text(size = 6))

ggsave(file.path(output_dir, "QC3_replicate_correlation.png"), p3, width = 10, height = 8)



#  Plot 4: Volcano Plot 
# flag genes as essential if FDR < 0.05 and logFC < -1
volcano_df <- all_results
volcano_df$sig <- ifelse(volcano_df$FDR < 0.05 & volcano_df$logFC < -1, 
                         "Essential", "Not Essential")
volcano_df$neg_log10_fdr <- -log10(volcano_df$FDR + 1e-300)

# get top 10 essential genes to label and label them 
top_labels <- volcano_df[volcano_df$sig == "Essential", ]
top_labels <- top_labels[order(top_labels$logFC), ][1:10, ]

p4 <- ggplot(volcano_df, aes(x = logFC, y = neg_log10_fdr, colour = sig)) +
  geom_point(alpha = 0.4, size = 0.8) +
  geom_text_repel(data = top_labels, aes(label = genes), size = 3, 
                  max.overlaps = 20) +
  scale_colour_manual(values = c("Essential" = "#d6604d", "Not Essential" = "grey70")) +
  geom_vline(xintercept = -1, linetype = "dashed", colour = "black") +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", colour = "black") +
  labs(title = "Volcano Plot  Colorectal Cancer Essential Genes",
       x = "Log2 Fold Change", y = "-log10(FDR)", colour = "") +
  theme_minimal()

ggsave(file.path(output_dir, "Plot4_volcano.png"), p4, width = 10, height = 7)


#  Plot 5: Heatmap  all common essential genes 
# get gene lfc matrix for common essential genes only
heatmap_genes <- common_essentials_df$gene
heatmap_matrix <- gene_lfc_matrix[rownames(gene_lfc_matrix) %in% heatmap_genes, ]

# clean column names for display away with protocol and  batch info
colnames_clean <- gsub("-311Cas9|_311Cas9|311Cas9|-Avana.*|_Avana.*", "", colnames(heatmap_matrix))
colnames_clean <- gsub("_p[0-9]_batch[0-9]|_batch[0-9]", "", colnames_clean)
colnames(heatmap_matrix) <- colnames_clean

#pheatmap that saves directly to file 
pheatmap(heatmap_matrix,
         color = colorRampPalette(c("#2166ac", "white", "#d6604d"))(100),
         cluster_rows = TRUE,
         cluster_cols = TRUE,
         show_rownames = TRUE,
         show_colnames = TRUE,
         fontsize_row = 6,
         fontsize_col = 7,
         main = "Heatmap  Common Essential Genes",
         filename = file.path(output_dir, "Plot5_heatmap.png"),
         width = 12, height = 14)


#  Plot 6: Bar chart top 30 common essential genes
top30 <- head(common_essentials_df, 30)

p6 <- ggplot(top30, aes(x = reorder(gene, n_cell_lines_essential), 
                        y = n_cell_lines_essential, fill = n_cell_lines_essential)) +
  geom_bar(stat = "identity") +
  scale_fill_gradient(low = "#2166ac", high = "#d6604d") +
  coord_flip() +
  labs(title = "Top 30 Common Essential Genes",
       x = "Gene", y = "Number of Cell Lines Essential In", fill = "") +
  theme_minimal() +
  theme(legend.position = "none",
        axis.text.y = element_text(size = 8))

ggsave(file.path(output_dir, "Plot6_common_essentials_bar.png"), p6, width = 8, height = 10)

# Plot 7: UMAP  cell lines coloured by batch
# transpose so rows = cell lines, cols = genes
umap_input <- t(gene_lfc_matrix)

# run umap
# setting see for reproducibility
set.seed(42)
umap_result <- umap(umap_input)

# build plot dataframe
umap_df <- data.frame(
  UMAP1      = umap_result$layout[, 1],
  UMAP2      = umap_result$layout[, 2],
  sample     = rownames(umap_input),
  batch      = ifelse(grepl("batch2", rownames(umap_input)), "Batch 2", "Batch 3")
)

# add clean cell line names
umap_df$cell_line <- gsub("-311Cas9|_311Cas9|311Cas9|-Avana.*|_Avana.*", "", umap_df$sample)
umap_df$cell_line <- gsub("_p[0-9]_batch[0-9]|_batch[0-9]|_Rep.*", "", umap_df$cell_line)


# code for plotting the umap 
p7 <- ggplot(umap_df, aes(x = UMAP1, y = UMAP2, colour = batch, label = cell_line)) +
  geom_point(size = 4, alpha = 0.8) +
  geom_text_repel(size = 3, max.overlaps = 20) +
  scale_colour_manual(values = c("Batch 2" = "#2166ac", "Batch 3" = "#d6604d")) +
  labs(title = "UMAP  Cell Lines Coloured by Batch",
       x = "UMAP1", y = "UMAP2", colour = "Batch") +
  theme_minimal()

ggsave(file.path(output_dir, "Plot7_umap_batch.png"), p7, width = 8, height = 7)


# Plot 8: MLE and common essential genes venn diagram
# load MAGeCK results
mageck_path <- "C:/Users/khali/CRISPR_thesis/Galaxy46-[MAGeCK mle on dataset 44 and 45_ Gene Summary (MLE)] (2).tabular"
mageck_results <- read_tsv(mageck_path, show_col_types = FALSE)
dim(mageck_results)


# load edgeR common essentials
edger_path <- "C:/Users/khali/CRISPR_thesis/results/common_essentials.csv"
common_essentials_df <- read_csv(edger_path, show_col_types = FALSE)
dim(mageck_results)

# print top 3 hits per cell line
cell_lines_reliable <- c("RKO", "SW403", "LS180", "GP5d", "C2BBE1", "HT115")

for (cl in cell_lines_reliable) {
  beta_col <- paste0(cl, "|beta")
  fdr_col  <- paste0(cl, "|fdr")
  top3 <- mageck_results[order(mageck_results[[beta_col]]), c("Gene", beta_col, fdr_col)]
  top3 <- head(top3[top3[[fdr_col]] < 0.05, ], 3)
  print(top3)
}

# clean edgeR gene names to match MAGeCK format
common_essentials_df$gene_clean <- gsub(" \\(.*\\)", "", common_essentials_df$gene)

# find MAGeCK essential genes per cell line
cell_lines <- c("RKO", "SW403", "LS180", "CL11", "GP5d", "C2BBE1", "SW48", "SW948", "HT115")

mageck_essential_list <- lapply(cell_lines, function(cl) {
  beta_col <- paste0(cl, "beta")
  fdr_col  <- paste0(cl, "fdr")
  mageck_results$Gene[mageck_results[[beta_col]] < 0 &
                        mageck_results[[fdr_col]] < 0.05]
})
names(mageck_essential_list) <- cell_lines

# check how many essential genes MAGeCK found per cell line
sapply(mageck_essential_list, length)


# find common essentials using threshold of 3 cell lines
mageck_gene_counts <- table(unlist(mageck_essential_list))
mageck_common <- names(mageck_gene_counts[mageck_gene_counts >= 5])

# build venn diagram comparing mageck and edger
venn_list <- list(
  edgeR  = common_essentials_df$gene_clean,
  MAGeCK = mageck_common
)

p8 <- ggvenn(venn_list,
             fill_color = c("#2166ac", "#d6604d"),
             stroke_size = 0.5,
             set_name_size = 5) +
  labs(title = "Overlap between edgeR and MAGeCK MLE Essential Genes")

# save plot
output_dir <- "C:/Users/khali/CRISPR_thesis/results/plots"
ggsave(file.path(output_dir, "Plot8_venn_diagram.png"), p8, width = 8, height = 7)


# show overlapping genes
overlap <- intersect(common_essentials_df$gene_clean, mageck_common)
length(overlap)
print(overlap)

# save overlap results
overlap_df <- data.frame(gene = overlap)
write.csv(overlap_df,
          "C:/Users/khali/CRISPR_thesis/results/edgeR_MAGeCK_overlap.csv",
          row.names = FALSE)
