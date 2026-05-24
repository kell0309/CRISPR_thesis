# The plots can be ran in sequences but no variable is printed directly
# This is meant to be used as an example for how our data was processed
# and printed and is only advised as a recommendation
# There is repetitive parts of the function more as part of the guiding
# If you run this code, please run small chunks at a time and ensure
# you only print was you want to visualise
# B next to a library means it needs biomanager to be installed
# use this to install biomanager and then the libraries if you
# do not have them already: install.packages("BiocManager")
# BiocManager::install()
library("pheatmap")
library("ggVennDiagram")
library("clusterProfiler") #B
library("corrplot")
library("org.Hs.eg.db") #B
library("enrichplot") #B
library("ggplot2")
library("patchwork")
library("ggrepel")
library("dplyr")

#artificial essential matrix by using the overlap of statistically significant genes
# and the pseudo essentials.
matx<- norplres[norplres$sgRNA %in% plasmid_essen$sgRNA, ]



# simple function for the ven diagrams used for the two datasets 2 versions
# one for essential genes, one for statistically important genes
# these functions are used to give a significant gene
# For basic testing sample and plasmid controls were used

#this function is meant to be used on only the essentials of Achilles comparing
#the gene columns
#And whatever input you want. I used all of my data in these comparisons
vennextr <- function(sang_data, column1, column2) {
  #we remove any shared sgRNAs with the same gene
  plas <- unique(sang_data[[column1]])
  #Here we need to remove the numeric values next to the genes
  samp <- sub(" \\(.*\\)", "", common_essentials[[column2]])
  #defining the two columns
  venno <- list( "Sanger" = plas, "Achilles" = samp)
  ggVennDiagram(venno, )
}
#same process as the function above but instead for all results
venncomplete <- function(sang_data, column1, column2) {
  #we remove any shared sgRNAs with the same gene
  plas <- unique(sang_data[[column1]])
  #Here we need to remove the numberic values next to the genes
  samp <- sub(" \\(.*\\)", "", edgeR_all_results[[column2]])
  #defining the two columns
  venno <- list( "Sanger" = plas, "Achilles" = samp)
  ggVennDiagram(venno)
}

# Here is a basic preparation of the data from achilles to remove the numeric value
# in addition we find the overlap and lastly we used the single gene function
# to the two different genes variables norplres (sanger) and plasmid_essen
# this chunk was made after the venn diagrams below but placed to top
achilles_uniq<-sub(" \\(.*\\)", "", edgeR_all_results$genes)
# the intersect between the Achilles Results and the Sanger results
common_genes<- intersect(achilles_uniq, unique(norplres$gene))
#really not needed this is how MCM6 was identified
  #single_gene<- intersect(sub(" \\(.*\\)", "", common_essentials$gene), unique(norplres$gene))
#common_genes
#single_gene


# shared essentials with statistical significance 0
vennextr(matx,"gene", "gene")

#shared essentials 0
vennextr(plasmid_essen, "gene", "gene")
# essential vs total number overlap
#sanger unique vs complete dataset achilles
venncomplete(matx , "gene", "genes")
# only one is not overlapping
# Achilles essentials vs statistically significant sanger
vennextr(norplres, "gene", "gene")
# 1 shared MCM6

#comparison of all genes
venncomplete(norpl$genes, "gene", "genes")
#results achilles: sanger 849, achilles 1364, shared 17160



# removing all sgRNAs that have a more highly expressed sgRNA with the same gene
gene_level <- norplres %>%
  #making groups by genes
  group_by(gene) %>%
  #taking only the highest log fold change
  slice_max(order_by = abs(logFC), n = 1) %>%
  #returning results without the groupings
  ungroup()
# removing all the insigificant genes of the achilles results
gene_level_edgeR <- edgeR_all_results %>%
  #removing the arithmetic brachets present in EdgeR
  mutate(genes = sub(" \\(.*\\)", "", genes)) %>%
  # filtering according to FDR
  filter(FDR < 0.05) %>%
  # same process as gene_level
  group_by(genes) %>%
  slice_max(order_by = abs(logFC), n = 1) %>%
  ungroup()


#preparing the data for the heatmap using LFC from gene_level
Heat_mat <- gene_level[, "logFC", drop = FALSE]
#putting back the rownames
rownames(Heat_mat) <- gene_level$gene
# turning heatrow into a matrix
Heat_mat <- as.matrix(Heat_mat)

#taking only the top values for the heatmap
topheat<- Heat_mat[order(abs(Heat_mat[,1]), decreasing = TRUE)[1:50], , drop = FALSE]
#AI was used to help guide through the function including colors and guide through
#The structure of the graph
pheatmap(
  # taking the top 50 genes
  topheat,
  # no clustering
  cluster_cols   = FALSE,
  # the two colours for positive and negative LFCs and the 
  color          = colorRampPalette(c("#2166AC", "white", "#D6604D"))(100),
  # font size
  fontsize_row   = 7,
  #title
  main = "Sanger top 50 Essential gene by LFC"
  )

# Below is a function used for the positive log fold changes but nothing was used
#Heat_mat_pos <- Heat_mat[Heat_mat[,1] > 0, , drop = FALSE]

# positive heatmap
#pheatmap(
#  Heat_mat_pos,
#  cluster_cols = FALSE, color = colorRampPalette(c("white", "#D6604D"))(100),
#  fontsize_row = 7, main = "Positive fold changes")

# the Volcano plot is based on the Volcano plot made for the Achilles code
# making a new variable identical to the Sanger dataset
volcano_df <- norplres
# taking information by a statistically significant FDR and LFC
volcano_df$sig <- ifelse(volcano_df$FDR < 0.05 & volcano_df$logFC < -1, "Negative LFC", "Positive LFC")
# Small number added to prevent an FDR of 0 from causing errors
volcano_df$neg_log10_fdr <- -log10(volcano_df$FDR + 1e-300)

#Labels and extaction of th top 20 genes, this is used in the volcano
#however it was not important in the final results.
top_labels <- volcano_df[volcano_df$sig == "Negative LFC", ]
top_labels <- top_labels[order(top_labels$logFC), ][1:20, ]

#defining the plot and which data each axis used
pheat <- ggplot(volcano_df, aes(x = logFC, y = neg_log10_fdr, colour = sig)) +
  #Defining size
  geom_point(alpha = 0.4, size = 0.8) +
  # Limiting overlap and defining data to be used
  geom_text_repel(data = top_labels, aes(label = gene), size = 3, 
                  max.overlaps = 20) +
  # defining the two groups and how they are named between the positive
  # and negative lfc groups
  scale_colour_manual(values = c("Negative LFC" = "#d6604d", "Positive LFC" = "grey70")) +
  geom_vline(xintercept = -1, linetype = "dashed", colour = "black") +
  #defining the intercepts and the lines
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", colour = "black") +
  # title x label and y label
  labs(title = "Sanger Volcano Plot LFC distribution",
       x = "Log2 Fold Change", y = "-log10(FDR)", colour = "") +
  theme_minimal()
#calling upon it
pheat



#Enrichment Analysis function was not used as much as originally planned
enrichpaths<- function(typeset) {
  # dataset and function
  # basic method similar to the pathway analysis but requiring
  # keytype to be explicitly named
  enrichGO(gene = typeset,
            # data structure and database
            OrgDb = org.Hs.eg.db, keyType = "ENTREZID", ont = "BP",
            pAdjustMethod = "BH", pvalueCutoff = 0.05)
}
# this was used earlier to extract the plasmid essentials,
# unique is redundant and no information was found ti was experimented 
# with insignificantly
as.character(unique(plasmid_essen$gene))
# this is a dataframe identical to the Sanger genes, with Entrez IDs
# same process as in the achilles dataset above
gene.df <- bitr( gene_level$gene,
                 fromType = "SYMBOL", toType = "ENTREZID",
                 OrgDb = org.Hs.eg.db)


# these were made for plot experimentation essent_enrich
# however became the enrichment analysis results later for the Sanger dataset
Essent_enrich<-enrichpaths(as.character(unique(gene.df$ENTREZID)))
enrich_simple<-enrichpaths(as.character(unique(norplres$gene[norplres$logFC > 0.8])))

# check for essent_enrichment
#Essent_enrich
#unique(norplres$gene[norplres$logFC > 0.9])
#unique(norplres$gene[norplres$logFC < -4])



# experimentation plots, not used
#dotplot(Essent_enrich, showCategory = 10, title = "")

#cnetplot(Essent_enrich, showCategory = 20, title = "")

#dotplot(enrich_simple, showCategory = 20, title = "")


#this part had a lot of experimentation, however these were the functions ran
#for the research paper so no clean up was made outside the comments
#Correlation matrices used for the matrix first by extracting the only the names
essential_Id <- norpl$genes$gene %in% plasmid_essen$gene
# matching the essentials and non essentials
norpl_essential <- norpl[essential_Id, ]
norpl_noness    <- norpl[!essential_Id, ]
# using a counts per million for both datasets
cormat_ess    <- cpm(norpl_essential, log = TRUE)
cormat_noness <- cpm(norpl_noness, log = TRUE)
# non essentials defines the moethod in the correlation plot
cormat <- cpm(norpl, log = TRUE)
#essentials uses pearson's here while
cor_ess <- cor(cormat_ess[, -1], method = "pearson")

#setting it so two plots are visualised
par(mfrow = c(1, 2))
# All Sanger cell lines compared using pearson's correlation
# only the top function was used
cor_mat <- cor(cormat[, -1], method = "pearson")
corrplot(cor_mat, method = "color", type = "upper",
         tl.col = "black", tl.srt = 45, tl.cex = 0.7,
         addCoef.col = "black", number.cex = 0.6,
         col = colorRampPalette(c("#2166AC", "white", "#D6604D"))(100),
         mar = c(2, 0, 3, 0),
         title = "total Cell line correlations", line = 0)

corrplot(cor_ess, method = "color", type = "upper",
         tl.col = "black", tl.srt = 45, tl.cex = 0.7,
         addCoef.col = "black", number.cex = 0.6,
         col = colorRampPalette(c("#2166AC", "white", "#D6604D"))(100),
         mar = c(2, 0, 3, 0),
         title = "Essential genes cell line correlations", line = -2)
# we set it one to 1 plot to be visualised         
par(mfrow = c(1, 1))

# Non-essential not significant only 2 genes
cor_noness <- cor(cormat_noness[, -1], method = "pearson")
corrplot(cor_noness, method = "color", type = "upper",
         tl.col = "black", tl.srt = 45, tl.cex = 0.7,
         addCoef.col = "black", number.cex = 0.6,
         col = colorRampPalette(c("#2166AC", "white", "#D6604D"))(100),
         main = "Non-essential genes - cell line correlations")






# intersect of the shared genes between the two datasets
shared_genes <- intersect(gene_level$gene, gene_level_edgeR$genes)
#merging the datasets into one dataframe to the LFC and the genes
merged <- merge(
  #extracting the columns from both datasets with genes and LFC
  gene_level[gene_level$gene %in% shared_genes, c("gene", "logFC")],
  gene_level_edgeR[gene_level_edgeR$genes %in% shared_genes, c("genes", "logFC")],
  #matching them based on the gene columns
  by.x = "gene", by.y = "genes")
#column names set as the datasets
colnames(merged) <- c("gene", "Sanger", "Achilles")

# filtering so that only the concordant genes are preserved
merged_filt <- merged %>%
  mutate(concordant = sign(Sanger) == sign(Achilles)) %>%
  #note the 0 is not important as specifically 0 it is a remnant
  #of experimentation
  filter(concordant, abs(Sanger) > 0, abs(Achilles) > 0)

#separation for the pathway enrichment no pathways were found
#merged_paths <- merged_filt[, -4]
#merged_paths <- bitr(as.character(merged_paths$gene),
#       fromType = "SYMBOL", toType = "ENTREZID",
#       OrgDb    = org.Hs.eg.db)

heat_mat_shared <- as.matrix(merged_filt[, c("Sanger", "Achilles")])
rownames(heat_mat_shared) <- merged_filt$gene
# This is the shared genes LFC, it was made after the spearman's
# below to extract the significant overlapping genes
pheatmap(heat_mat_shared,
         cluster_cols = FALSE,
         color        = colorRampPalette(c("#2166AC", "white", "#D6604D"))(100),
         fontsize_row = 7,
         angle_col = 0,
         main         = "Shared genes - LFC comparison")


#ín addition after measuring the raw similarity here is a spearman's test
spearman_result <- cor.test(merged$Sanger, merged$Achilles, method = "spearman")


#spearman's plot and defining datasets
ggplot(merged, aes(x = Sanger, y = Achilles)) +
  #defining points
  geom_point(alpha = 0.4, size = 1) +
  # making the line red and smooth
  geom_smooth(method = "loess", color = "red") +
  # defining the two axis y and x with a centre of 0
  # originally it was going to be modified but it was apt for comparison
  geom_hline(yintercept = 0, linetype = "dashed") +
  geom_vline(xintercept = 0, linetype = "dashed") +
  #defimimg the two texts
  annotate("text", x = Inf, y = Inf, 
           #extracting from the spearman's results the ρ value and the
           #p-value
           label = paste0("ρ = ", round(spearman_result$estimate, 3),
                          "\np = ", signif(spearman_result$p.value, 3)),
           # the size would have been adjusted but it was a good location
           hjust = 1.1, vjust = 1.5, size = 4) +
  # defining the title and the X and Y label
  labs(title = "Spearman's Correlation: Sanger vs Achilles",
       x = "Sanger LFC", y = "Achilles LFC") +
  theme_bw()



#taking the shared significant genes and going through the process to find the 
# common term and pathway enrichment
entrez_shared <- bitr(as.character(merged$gene),
                      fromType = "SYMBOL", toType = "ENTREZID",
                      OrgDb    = org.Hs.eg.db)

enrich_shared <- enrichpaths(as.character(entrez_shared$ENTREZID))
#enrichment analysis of shared datasets, nothing significant
kegg_shared <- enrichKEGG(gene          = entrez_shared$ENTREZID,
                          organism      = "hsa",
                          pAdjustMethod = "BH",
                          pvalueCutoff  = 0.05)




# this is renaming the Symbols for the genes to entrez ids from the
# human library in the Achilles dataset
gene_entrez_edgeR <- bitr(as.character(gene_level_edgeR$genes),
                          fromType = "SYMBOL",
                          toType   = "ENTREZID",
                          OrgDb    = org.Hs.eg.db)
# simple definition of a pathway analysis data and column of the IDs
# followed by humans, padjust method in this case BH and the normal
# cuttoff value of 0.05 to ensure significance
# this is used at the very end for the 4 ways plot
pathhhh_edgeR <- enrichKEGG(gene          = gene_entrez_edgeR$ENTREZID,
                            organism      = "hsa",
                            pAdjustMethod = "BH",
                            pvalueCutoff  = 0.05)

# Run KEGG pathway enrichment for the Sanger dataset same as the achilles
# the dataframe with entrez ids was made earlier for experimentation
pathhhh <- enrichKEGG(
  gene = gene.df$ENTREZID,
  organism = "hsa",
  pAdjustMethod = "BH",
  pvalueCutoff = 0.05
)
# merging the sanger and the achilles pathway enrichment into variables used for
# the final 4 way plot into 2 identical boxplots
achilpat<- barplot(pathhhh_edgeR,  title = "Pathway analysis Achiles")+
  theme(axis.text.x = element_text(size = 10))
sangerpat<- barplot(pathhhh, title = "Pathway analysis Sanger")+
  theme(axis.text.x = element_text(size = 10))

# this was a 2 way plot constructed to ensure it works
# achilpat / sangerpat + plot_annotation(tag_levels = list(c("A)","B)")))


# below is the enrichment shared and the 
enr_shared<- barplot(enrich_shared, title = "Enrichment shared genes")+
  theme(axis.text.x = element_text(size = 10))
path_shared<- barplot(kegg_shared, title = "Pathway analysis shared genes") +
  theme(axis.text.x = element_text(size = 10))
#these were the two original plots made for th one above, for confirmation
#enr_shared / path_shared + plot_annotation(tag_levels = list(c("A)","B)")))

# making the 4 way plot set up in the study based on location
(enr_shared | path_shared) /
  (achilpat | sangerpat) +
  #this defines each plot with the letter
  plot_annotation(tag_levels = "A")
