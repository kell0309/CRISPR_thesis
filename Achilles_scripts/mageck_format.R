# load libraries
library(dplyr)
library(stringr)

# load master matrix
colon_raw_counts <- readRDS("C:/Users/khali/CRISPR_thesis/colon_raw_counts.rds")
dim(colon_raw_counts)

# prepare MAGeCK counts matrix
mageck_input <- colon_raw_counts %>%
  select(-is_essential, -is_nonessential) %>%
  select(-contains("batch4"))

# rename Construct_Barcode to sgRNA
colnames(mageck_input)[1] <- "sgRNA"

# clean gene column 
mageck_input$gene <- gsub(" \\(.*\\)", "", mageck_input$gene)

# clean all column names - remove spaces and special characters
colnames(mageck_input) <- colnames(mageck_input) %>%
  str_replace_all(" ", "_") %>%
  str_replace_all("[^A-Za-z0-9_]", "") %>%
  str_replace_all("__", "_")

# fix duplicate pDNA batch3 column name
colnames(mageck_input)[colnames(mageck_input) == "Avana_4_Hu_pDNA_MAA40_93015_02pguL_batch3"][2] <- "Avana_4_Hu_pDNA_MAA40_93015_02pguL_v2_batch3"

# check final structure
dim(mageck_input)
colnames(mageck_input)          
head(mageck_input[, 1:5])

# save counts matrix
write.table(mageck_input,
            "C:/Users/khali/CRISPR_thesis/Achilles_scripts/mageck_counts.txt",
            sep = "\t", row.names = FALSE, quote = FALSE)

# build design matrix
design_matrix <- data.frame(
  Samples = c(
    "Avana4pDNA20160601311cas9_RepG09_batch2",
    "Avana4pDNA20160601311cas9_RepG10_batch2",
    "Avana4pDNA20160601311cas9_RepG11_batch2",
    "Avana4pDNA20160601311cas9_RepG12_batch2",
    "Avana_4_Hu_pDNA_MAA40_93015_02pguL_batch3",
    "Avana_4_Hu_pDNA_MAA40_93015_02pguL_v2_batch3",
    "Avana_4_Hu_pDNA_MAA40_93015_batch3",
    "Avana_4_Hu_pDNA_MAF34_112718_02pguL_batch3",
    "RKO311Cas9_RepA_p6_batch2",
    "RKO311Cas9_RepB_p6_batch2",
    "RKO311Cas9_RepC_p6_batch2",
    "SW403311Cas9_RepB_p5_batch2",
    "LS_180311Cas9_Rep_A_p4_batch2",
    "LS_180311Cas9_Rep_B_p4_batch2",
    "CL11311Cas9_RepB_p6_batch3",
    "GP5d311Cas9_RepA_p6_batch3",
    "GP5d311Cas9_RepB_p6_batch3",
    "C2BBE1311Cas9_Rep_A_p5_batch3",
    "C2BBE1311Cas9_Rep_B_p5_batch3",
    "C2BBE1311Cas9_Rep_C_p5_batch3",
    "SW48311Cas9_RepA_p6_batch3",
    "SW48311Cas9_RepB_p6_batch3",
    "SW948311Cas9_RepB_p6_batch3",
    "HT115Avana_4_Rep_A_p6_batch3",
    "HT115Avana_4_Rep_B_p6_batch3",
    "HT115Avana_4_Rep_C_p6_batch3"
  ),
  baseline = c(1,1,1,1, 1,1,1,1, 0,0,0,0,0,0, 0,0,0,0,0,0,0,0,0,0,0,0),
  RKO    =   c(1,1,1,1, 0,0,0,0, 1,1,1,0,0,0, 0,0,0,0,0,0,0,0,0,0,0,0),
  SW403  =   c(1,1,1,1, 0,0,0,0, 0,0,0,1,0,0, 0,0,0,0,0,0,0,0,0,0,0,0),
  LS180  =   c(1,1,1,1, 0,0,0,0, 0,0,0,0,1,1, 0,0,0,0,0,0,0,0,0,0,0,0),
  CL11   =   c(0,0,0,0, 1,1,1,1, 0,0,0,0,0,0, 1,0,0,0,0,0,0,0,0,0,0,0),
  GP5d   =   c(0,0,0,0, 1,1,1,1, 0,0,0,0,0,0, 0,1,1,0,0,0,0,0,0,0,0,0),
  C2BBE1 =   c(0,0,0,0, 1,1,1,1, 0,0,0,0,0,0, 0,0,0,1,1,1,0,0,0,0,0,0),
  SW48   =   c(0,0,0,0, 1,1,1,1, 0,0,0,0,0,0, 0,0,0,0,0,0,1,1,0,0,0,0),
  SW948  =   c(0,0,0,0, 1,1,1,1, 0,0,0,0,0,0, 0,0,0,0,0,0,0,0,1,0,0,0),
  HT115  =   c(0,0,0,0, 1,1,1,1, 0,0,0,0,0,0, 0,0,0,0,0,0,0,0,0,1,1,1)
)

# verify samples match between counts and design matrix
matched <- intersect(colnames(mageck_input)[3:ncol(mageck_input)], design_matrix$Samples)
length(matched)       


# preview design matrix
head(design_matrix)
dim(design_matrix)

# save design matrix
write.table(design_matrix,
            "C:/Users/khali/CRISPR_thesis/mageck_design.txt",
            sep = "\t", row.names = FALSE, quote = FALSE)
