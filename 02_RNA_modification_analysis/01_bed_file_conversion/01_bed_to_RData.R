# Specify the projects dir containing all the functions of scRNA-seq analysis pipeline
source("~/Documents/Projects/my_projects/m6A_RNA_methylation_in_colorectal_cancer/library_and_packages/library_handler.R", .GlobalEnv)

# Functions to perfrom analysis tasks
funcset <- list.files("~/Documents/Projects/my_projects/nanopore_RNA_DNA_analysis_functions", pattern = "*.R$", full.names = TRUE)
sapply(funcset, source, .GlobalEnv)

# Directory to store the results
res_dir = "~/Documents/Projects/m6A_RNA_modification_project/modification_analysis_output/"

store_folder = "01_converted_modified_base_count"

if(!dir.exists(str_c(res_dir))){
  dir.create(str_c(res_dir), showWarnings = TRUE, recursive = FALSE, mode = "0777")
}

# Data directroy
data_dir <- "~/Documents/Projects/m6A_RNA_modification_project/samples_bed_files/"
# Example data directory
# data_dir <- "~/Documents/Projects/m6A_RNA_modification_project/example_bed_files/"

file_paths <- list.files(data_dir, pattern = "*.bed", full.names = TRUE)

for (i in 1:length(file_paths)){
  
  sample_name <- sub("\\..*", "", basename(file_paths[i]))
  
  # Read and convert the ".bed" file 
  bed_conversion(file_paths[i], res_dir, store_folder, sample_name)
}
