library(edgeR)
library(dplyr)

#below I just made a single function to normalise by TMM all the data plas_def being the data set
normalisation<- function(plas_def) {
  con_DGE<-DGEList(counts = plas_def, group = factor(ifelse(grepl("plasmid", 
  #Sapply to make it so the columns have to be numeric, this removes the non count columns
  colnames(plas_def[sapply(plas_def, is.numeric)])), "plasmid", "sample")))
  #finally running the data through TMM
  norfact<-calcNormFactors(con_DGE, method = "TMM", refColumn = 1)
  #cpm(norfact, normalized.lib.sizes = TRUE, log =FALSE)
}
#For the time based cell lines they will each be separated by their time value before they are normalised
time_seps<-function(x) {
  cols_tim<-(which(colnames(x) == "plasmid") + 1):ncol(x)
  sam_tim<-x[, cols_tim]
  facto_tim<-sub("\\d+$", "", colnames(sam_tim))
  
  split_tim <- split.default(sam_tim, facto_tim)
  lapply(split_tim, function(group) cbind(x[, 1:which(colnames(x) == "plasmid")], group))
}
#normalisaiton of the non time based data
norcs<-normalisation(non_CRISPR_C6596666.sample)
norpl<-normalisation(non_ERS717283.plasmid)
#nestled in here both the normalisaiton and the time based data, cuts on lines
#they will not be used but a function for them was created
#poscs<-lapply(time_seps(pos_CRISPR_C6596666.sample), normalisation)
#pospl<-lapply(time_seps(pos_ERS717283.plasmid), normalisation)

#masterframe is self explanatory the DGElist
folds_data<-function(masterframe) {
  #design matrix for the estimatedistribution
  #Instead of rewriting group, just pull it from your DGEList
  design <- model.matrix(~ masterframe$samples$group)
  #estimate distribution
  con_sar <- estimateDisp(masterframe, design)

  fit <- glmFit(con_sar, design)
  lrt <- glmLRT(fit, coef = 2)

  results <- topTags(lrt, n = Inf, p.value = 0.05) |>
    as.data.frame() |>
    arrange(FDR)
}
#printing out the statistical tests
norcsres<-folds_data(norcs)
norplres<-folds_data(norpl)

#Extracting the number of cell lines where each sgRNA is less thatn one in expression as well
#As a list of which cell lines they are in
essential_genes<- function(masterframe, threshhold) {
  #we are stanrding after the 3rd and 4th columns to only take numberic values
  masterframe <- masterframe[masterframe[, 3] >= 10, ]
  colstart<- masterframe[, 4:ncol(masterframe)]
  #we turn rowSums into a true or false
  Essential_Number <- rowSums(colstart < 1, na.rm = TRUE)
  #each of the cell lines that has a 0 count on the gene
  Essential_Names<- apply(colstart < 1, 1, function(row) names(which(row)))
  #taking the first 3 columns from the master frame and binding the new columns
  binding<- cbind(masterframe[, 1:3], Essential_Number, Essential_Names = I(Essential_Names))
  #we only care for keeping the genes with atleast a single cell line
  binding[binding$Essential_Number >= ncol(colstart) * threshhold, ]
}
# the two datasets of the essential genes holding the information 
sample_essen<-essential_genes(non_CRISPR_C6596666.sample, 0.1)
# the reasion for the threshhold being variable is to take into account
# the difference not only in size but diversity as sample_essen is made up of
# cell lines from the same origins
plasmid_essen<-essential_genes(non_ERS717283.plasmid, 0.5)

