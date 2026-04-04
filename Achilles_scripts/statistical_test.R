# 03_logFold_R

# loading the library
library(edgeR)

# Creating a function for the statistical test
fold_change <- function(df, cellline_name) {
  
  rep_cols <- colnames(df)[grepl("_Rep", colnames(df))]
  count_matrix <- as.matrix(df[, c("plasmid", rep_cols)])
  rownames(count_matrix) <- df$sgRNA
  group <- factor(c("plasmid", rep("screen", length(rep_cols))))
  
  dge <- DGEList(counts = count_matrix, group = group)
  dge <- calcNormFactors(dge, method = "TMM")
  design <- model.matrix(~group)
  
  if (length(rep_cols) < 2) {
    
    # single replicate  use fixed BCV with glmLRT
    cat("\n", cellline_name, "— single replicate, using fixed BCV = 0.4\n")
    bcv <- 0.4
    fit <- glmFit(dge, design, dispersion = bcv^2)
    test <- glmLRT(fit)
    
  } else {
    
    # multiple replicates  estimate dispersion and use glmQLFTest
    cat("\n", cellline_name, "— estimating dispersion from replicates\n")
    dge <- estimateDisp(dge, design)
    fit <- glmQLFit(dge, design)
    test <- glmQLFTest(fit)
    
  }
  
  # extract results
  results <- topTags(test, n = Inf, sort.by = "PValue")$table
  results$sgRNA <- rownames(results)
  
  # merge with gene names
  final <- merge(
    df[, c("sgRNA", "gene")],
    results,
    by = "sgRNA"
  )
  
  # flag which method was used
  final$dispersion_method <- ifelse(length(rep_cols) < 2,
                                    "fixed BCV glmLRT",
                                    "estimated glmQLFTest")
  return(final)
}

counts_list_fc <- mapply(
  function(df, name) fold_change(df, name),
  counts_list,
  names(counts_list),
  SIMPLIFY = FALSE
)

cat("\nDone! Fold change calculated for", length(counts_list_fc), "cell lines\n")
print(head(counts_list_fc[[1]]))
