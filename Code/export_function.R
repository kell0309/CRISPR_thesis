# this is a tiny amount of code only meant for exporting the results based on the script directory
results_dir <- file.path(dirname(dirname(rstudioapi::getSourceEditorContext()$path)), "Sanger_Results")
# creating a small function that will go through all my code
save_csv <- function(data, fname) {
  write.csv(data, file.path(results_dir, fname), row.names = FALSE)
}
write.csv(plasmid_essen, file.path(results_dir, "Essential_ERS717283.csv"), row.names = FALSE)
write.csv(norcsres, file.path(results_dir, "edgeR_ERS717283.csv"), row.names = FALSE)
write.csv(norcs, file.path(results_dir, "DGE_ERS717283.csv"), row.names = FALSE)

write.csv(sample_essen, file.path(results_dir, "Essential_CRISPR_C6596666.csv"), row.names = FALSE)
write.csv(norplres, file.path(results_dir, "edgeR_C6596666.csv"), row.names = FALSE)
write.csv(norpl, file.path(results_dir, "DGE_C6596666.csv"), row.names = FALSE)
