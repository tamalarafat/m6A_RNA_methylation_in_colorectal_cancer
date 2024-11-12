# Specify the projects dir containing all the functions of scRNA-seq analysis pipeline
source("~/Documents/Projects/my_projects/m6A_RNA_methylation_in_colorectal_cancer/library_and_packages/library_handler.R", .GlobalEnv)

# Functions to perfrom analysis tasks
funcset <- list.files("~/Documents/Projects/my_projects/nanopore_RNA_DNA_analysis_functions", pattern = "*.R$", full.names = TRUE)
sapply(funcset, source, .GlobalEnv)

# Directory to store the results
res_dir = "~/Documents/Projects/m6A_RNA_modification_project/modification_analysis_output/"

store_folder = "09_consensus_motif_plots"

# Create the output directory if it doesn't exist
if (!dir.exists(file.path(res_dir, store_folder))) {
  dir.create(file.path(res_dir, store_folder), recursive = TRUE, mode = "0777")
}

store_dir = file.path(res_dir, store_folder)


# Load conditions from config.yaml
config_file_path <- "~/Documents/Projects/my_projects/m6A_RNA_methylation_in_colorectal_cancer/library_and_packages/configuration.yaml"

# Data directroy
file_paths <- "~/Documents/Projects/m6A_RNA_modification_project/modification_analysis_output/07_modification_sites_and_motifs/significant_mod_sites_signature_motifs.RData"

motif_data_list <- motif_data_per_condition(file_paths, config_file_path)

motif_data <- motif_data_list[1]

temp_motif_cols = motif_data_list[[2]]

analyte_type = motif_data_list[[3]]

# Get fasta subset and gene transcript table
generate_consensus_motif_logo(data_list = motif_data, motif_col_id = temp_motif_cols, store_dir = store_dir)

motif_base_proportion_plot(data_list = motif_data, motif_col_id = temp_motif_cols, analyte_type = analyte_type, store_dir = store_dir)
