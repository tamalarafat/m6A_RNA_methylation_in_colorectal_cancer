# Specify the projects dir containing all the functions of scRNA-seq analysis pipeline
source("~/Documents/Projects/my_projects/m6A_RNA_methylation_in_colorectal_cancer/library_and_packages/library_handler.R", .GlobalEnv)

# Functions to perfrom analysis tasks
funcset <- list.files("~/Documents/Projects/my_projects/nanopore_RNA_DNA_analysis_functions", pattern = "*.R$", full.names = TRUE)
sapply(funcset, source, .GlobalEnv)

# Directory to store the results
res_dir = "~/Documents/Projects/m6A_RNA_modification_project/modification_analysis_output/"

store_folder = "02_preprocessed_modification_sites"

# Data directroy
data_dir <- "~/Documents/Projects/m6A_RNA_modification_project/modification_analysis_output/01_converted_modified_base_count/"

file_paths <- list.files(data_dir, pattern = "*.RData", full.names = TRUE)

for (i in 1:length(file_paths)){
  
  sample_name <- sub("\\..*", "", basename(file_paths[i]))
  
  # Run the function with provided arguments
  filter_modification_sites(file_paths[i], retain_strand = "+", retain_Ndiff_percent = 25, retain_Ndel_percent = 25, append_data = FALSE, res_dir, store_folder)
}








