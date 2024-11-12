# Specify the projects dir containing all the functions of scRNA-seq analysis pipeline
source("~/Documents/Projects/my_projects/m6A_RNA_methylation_in_colorectal_cancer/library_and_packages/library_handler.R", .GlobalEnv)

# Functions to perfrom analysis tasks
funcset <- list.files("~/Documents/Projects/my_projects/nanopore_RNA_DNA_analysis_functions", pattern = "*.R$", full.names = TRUE)
sapply(funcset, source, .GlobalEnv)

# Directory to store the results
res_dir = "~/Documents/Projects/m6A_RNA_modification_project/modification_analysis_output/"

store_folder = "07_modification_sites_and_motifs"

# Create the output directory if it doesn't exist
if (!dir.exists(file.path(res_dir, store_folder))) {
  dir.create(file.path(res_dir, store_folder), recursive = TRUE, mode = "0777")
}

store_dir = file.path(res_dir, store_folder)

# Load conditions from config.yaml
config_file_path <- "~/Documents/Projects/my_projects/m6A_RNA_methylation_in_colorectal_cancer/library_and_packages/configuration.yaml"

# Data directroy
file_paths <- "~/Documents/Projects/m6A_RNA_modification_project/modification_analysis_output/05_modification_sites_test_statistics/modification_sites_Z_statistics.RData"

# Directory of the fasta file
fasta_dir = "~/Documents/Projects/m6A_RNA_modification_project/modification_analysis_output/06_fasta_and_gene_transcript_table/fasta_subset.fasta"

# Detect and store motif sequences to modification table
temp_df = find_motif(file_paths = file_paths, config_file_path = config_file_path, fasta_path = fasta_dir)

# Save the gene transcript table
save(temp_df, file = file.path(store_dir, "significant_mod_sites_motifs.RData"))
