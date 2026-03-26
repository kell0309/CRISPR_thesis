setwd("C:/Users/khali/Documents/Achilles/achilles_per_cellline")

getwd()  # confirm it changed

files <- list.files(pattern = "*.txt")
files    # confirm files are found

# Now merge
df_all <- Reduce(function(x, y) merge(x, y, by = c("sgRNA", "gene","plasmid")), df_list)

