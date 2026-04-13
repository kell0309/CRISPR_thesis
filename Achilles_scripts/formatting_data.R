# load the libraries in
library(readxl)

# load the file path
file_path_name <- readline(prompt = "Please enter the file path to your files: ")

# Read the file to  handle both Excel and CSV
if (str_ends(file_path_name, "\\.xlsx|\\.xls")) {
  data <- read_excel(file_path_name)
} else if (str_ends(file_path_name, "\\.csv")) {
  data <- read_csv(file_path_name)
} else {
  stop("File must be .xlsx, .xls or .csv")
}


# Extract specific columns from raw data
cat("\n next ste extracting columns needed in sample info")

# Define the column

columns <- c(
  "DepMap_ID",
  "sample_collection_site",
  "primary_or_metastasis",
  "primary_disease",
  "Subtype"
)

sample_info <- raw_data[, cols_needed]