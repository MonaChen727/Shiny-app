# Load necessary libraries
library(data.table) # for fread()
library(here)       # for handling paths

here::i_am("scripts/gene_num.R")


# List all .gct.gz files starting with "gene_tpm_v10_" in the specified directory
files <- list.files(
  path = here::here("gtex_v10_shiny/data/raw_data"),
  pattern = "^gene_tpm_v10_.*\\.gct\\.gz$", 
  full.names = TRUE
)


#test all the expression files contain the same gene
# Function to get nrow and ncol of a .gct.gz file
get_nrow <- function(file) {
  data <- fread(file, skip = 2) # Skip the first two metadata lines in .gct files
  return(nrow(data))
}
get_ncol <- function(file) {
  data <- fread(file, skip = 2) # Skip the first two metadata lines in .gct files
  return(ncol(data))
}

# Apply the function to each file and store results
file_rows <- sapply(files, get_nrow)
file_cols <- sapply(files, get_ncol)

file_rows[file_rows != 59033]
## all the output is 59033

gene_n <- function(file) {
  data <- read.table(gzfile(file), sep = "\t", skip = 2, header = TRUE)
  return(length(unique(exp[,2])))  # Count unique gene names in the second column
}

# Apply the function to each file in the 'files' vector and store results
gene_numbers <- sapply(files, gene_n)
gene_numbers[gene_numbers != 57853]
## all the ouptut is 57853



