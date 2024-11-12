# Specify the projects dir containing all the functions of scRNA-seq analysis pipeline
source("~/Documents/Projects/my_projects/m6A_RNA_methylation_in_colorectal_cancer/library_and_packages/library_handler.R", .GlobalEnv)

# Functions to perfrom analysis tasks
funcset <- list.files("~/Documents/Projects/my_projects/nanopore_RNA_DNA_analysis_functions", pattern = "*.R$", full.names = TRUE)
sapply(funcset, source, .GlobalEnv)

# Directory to store the results
res_dir = "~/Documents/Projects/m6A_RNA_modification_project/modification_analysis_output/"

store_folder = "03_merged_samples_modification_sites"

# Create the output directory if it doesn't exist
if (!dir.exists(file.path(res_dir, store_folder))) {
  dir.create(file.path(res_dir, store_folder), recursive = TRUE, mode = "0777")
}

store_dir = file.path(res_dir, store_folder)

# Data directroy
data_dir <- "~/Documents/Projects/m6A_RNA_modification_project/modification_analysis_output/02_preprocessed_modification_sites/"

file_paths <- list.files(data_dir, pattern = "*.RData", full.names = TRUE)

# Run merging function
merged_data <- merge_all_samples(file_paths, c("ref_seq_id", "start_pos", "strand"))

# Save the final merged data
save(merged_data, file = file.path(store_dir, "joint_modification_sites.RData"))
