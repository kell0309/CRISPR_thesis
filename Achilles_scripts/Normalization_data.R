# 02_normalisation.R
library(tidyverse)
library(edgeR)
 """
 1.motivate why use TMM¨
 2. why we use CPM
 
 """
tmm_normalise <- function(df) {
  
  # step 1: find rep columns only (no plasmid, no mean)
  rep_cols <- colnames(select(df, contains("_Rep")))
  
  # step 2: build count matrix (rows = sgRNAs, cols = replicates)
  count_matrix <- as.matrix(select(df, all_of(rep_cols)))
  rownames(count_matrix) <- df$sgRNA
  
  # step 3: create edgeR object and calculate TMM factors
  dge <- DGEList(counts = count_matrix)
  dge <- calcNormFactors(dge, method = "TMM")
  
  # step 4: extract normalised counts (log2 CPM)
  norm_matrix <- cpm(dge, log = TRUE, prior.count = 1)
  norm_df <- as.data.frame(norm_matrix)
  norm_df$sgRNA <- rownames(norm_matrix)
  
  # step 5: join back with gene and plasmid columns
  result <- merge(
    select(df, sgRNA, gene, plasmid,),
    norm_df,
    by = "sgRNA"
  )
  
  return(result)
}

counts_list_norm <- lapply(counts_list, tmm_normalise)

cat("\nFirst file after normalisation:\n")
print(head(counts_list_norm[[1]]))
