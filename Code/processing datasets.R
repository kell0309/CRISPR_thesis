library("pheatmap")
library("ggVennDiagram")
library("clusterProfiler")
library("corrplot")
library("org.Hs.eg.db")

# simple function for the ven diagrams used for the two datasets 2 versions
# one for essential genes, one for statistically important genes
# these functions are used to give a significant gene
# For basic testing sample and plasmid controls were used

vennextr <- function(column) {
  #we remove any shared sgRNAs with the same gene
  plas <- unique(plasmid_essen[[column]])
  #Here we need to remove the numberic values next to the genes
  samp <- common_essentials[[column]] <- sub(" \\(.*\\)", "", common_essentials[[column]])
  venno <- list( "Sanger" = plas, "Achilles" = samp)
  ggVennDiagram(venno)
}

venncomplete <- function(sang_data, column1, column2) {
  #we remove any shared sgRNAs with the same gene
  plas <- unique(sang_data[[column1]])
  #Here we need to remove the numberic values next to the genes
  samp <- edgeR_all_results[[column2]] <- sub(" \\(.*\\)", "", edgeR_all_results[[column2]])
  venno <- list( "Sanger" = plas, "Achilles" = samp)
  ggVennDiagram(venno)
}

achilles_uniq<-sub(" \\(.*\\)", "", edgeR_all_results$genes)
common_genes<- intersect(achilles_uniq, unique(plasmid_essen$gene))
single_gene<- intersect(sub(" \\(.*\\)", "", common_essentials$gene), unique(plasmid_essen$gene))

common_genes


pheatmap(
  Heat_mat,
  #annotation_row = row_anno,
  cluster_cols   = FALSE,
  color          = colorRampPalette(c("#2166AC", "white", "#D6604D"))(100),
  fontsize_row   = 7,
  main = "Essential gene hits - LFC"
)



Heat_mat <- norplres[, "logFC", drop = FALSE]
rownames(Heat_mat) <- norplres$sgRNA  
Heat_mat <- as.matrix(Heat_mat)
row_anno <- as.data.frame(sample_essen)
row_anno <- row_anno[row_anno$sgRNA %in% rownames(Heat_mat), ]
row_anno <- row_anno[!duplicated(row_anno$sgRNA), ]
rownames(row_anno) <- row_anno$sgRNA
row_anno$sgRNA <- NULL

pheatmap(
  Heat_mat,
  #annotation_row = row_anno,
  cluster_cols   = FALSE,
  color          = colorRampPalette(c("#2166AC", "white", "#D6604D"))(100),
  fontsize_row   = 7,
  main = "Essential gene hits - LFC"
  )






#Enrichment Analysis

enrichpaths<- function(typeset) {
  # dataset and function
  enrichGO(gene = typeset,
           # data structure and database
           OrgDb = org.Hs.eg.db, keytype = "SYMBOL", ont = "BP",
           pAdjustMethod = "BH", pvalueCutoff = 0.05)
}
as.character(unique(plasmid_essen$gene))

Essent_enrich<-enrichpaths(as.character(unique(plasmid_essen$gene)))
enrich_simple<-enrichpaths(norcsres$gene[norcsres$logFC > 2])

unique(norcsres$gene[norcsres$logFC > 0.9])
unique(norcsres$gene[norcsres$logFC < -5])



dotplot(ress_enr, showCategory = 20, title = "")



#Correlation matrices
cormat <- cpm(norpl, log = TRUE)

# All cell lines
cor_mat <- cor(cormat[, -1], method = "pearson")
corrplot(cor_mat, method = "color", type = "upper",
         tl.col = "black", tl.srt = 45,
         addCoef.col = "black", number.cex = 0.6,
         col = colorRampPalette(c("#2166AC", "white", "#D6604D"))(100),
         main = "Cell line correlations")

# cell line  plasmid
cor_plasmid <- cor(cormat[, "plasmid"], cormat[, -1], method = "pearson")
corrplot(cor_plasmid, method = "color",
         tl.col = "black", tl.srt = 45,
         addCoef.col = "black", number.cex = 0.8,
         col = colorRampPalette(c("#2166AC", "white", "#D6604D"))(100),
         main = "Cell line vs plasmid correlations")



essential_Id <- norpl$genes$gene %in% plasmid_essen$gene

norpl_essential <- norpl[essential_Id, ]
norpl_noness    <- norpl[!essential_Id, ]

cormat_ess    <- cpm(norpl_essential, log = TRUE)
cormat_noness <- cpm(norpl_noness, log = TRUE)

# Essential
cor_ess <- cor(cormat_ess[, -1], method = "pearson")
corrplot(cor_ess, method = "color", type = "upper",
         tl.col = "black", tl.srt = 45, tl.cex = 0.7,
         addCoef.col = "black", number.cex = 0.6,
         col = colorRampPalette(c("#2166AC", "white", "#D6604D"))(100),
         main = "Essential genes - cell line correlations")

# Non-essential
cor_noness <- cor(cormat_noness[, -1], method = "pearson")
corrplot(cor_noness, method = "color", type = "upper",
         tl.col = "black", tl.srt = 45, tl.cex = 0.7,
         addCoef.col = "black", number.cex = 0.6,
         col = colorRampPalette(c("#2166AC", "white", "#D6604D"))(100),
         main = "Non-essential genes - cell line correlations")


plasmid_essen
#two results from the ven diagrams
vennextr("gene")
venncomplete(plasmid_essen, "gene", "genes")
venncomplete(norplres, "gene", "genes")
