# Specify the projects dir containing all the functions of scRNA-seq analysis pipeline
source("~/Documents/Projects/my_projects/m6A_RNA_methylation_in_colorectal_cancer/library_and_packages/library_handler.R", .GlobalEnv)

# Functions to perfrom analysis tasks
funcset <- list.files("~/Documents/Projects/my_projects/nanopore_RNA_DNA_analysis_functions", pattern = "*.R$", full.names = TRUE)
sapply(funcset, source, .GlobalEnv)

# Directory to store the results
res_dir = "~/Documents/Projects/m6A_RNA_modification_project/modification_analysis_output/"

store_folder = "05_modification_sites_test_statistics"

# Create the output directory if it doesn't exist
if (!dir.exists(file.path(res_dir, store_folder))) {
  dir.create(file.path(res_dir, store_folder), recursive = TRUE, mode = "0777")
}

store_dir = file.path(res_dir, store_folder)

# Load conditions from config.yaml
config_file_path <- "~/Documents/Projects/my_projects/m6A_RNA_methylation_in_colorectal_cancer/library_and_packages/configuration.yaml"

# Data directroy
data_dir <- "~/Documents/Projects/m6A_RNA_modification_project/modification_analysis_output/04_filtered_modification_sites/"

file_paths <- list.files(data_dir, pattern = "*.RData", full.names = TRUE)

# Run Z-statistics function
significant_mod_sites <- compute_Zstatistics(file_paths, config_file_path, filter_by_pval = 0.05)

# Save the table with z-statistics data
save(significant_mod_sites, file = file.path(store_dir, "modification_sites_Z_statistics.RData"))
