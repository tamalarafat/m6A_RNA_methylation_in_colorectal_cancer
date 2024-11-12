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
file_paths <- "~/Documents/Projects/m6A_RNA_modification_project/modification_analysis_output/07_modification_sites_and_motifs/significant_mod_sites_motifs.RData"

# Process the input file to detect signature motifs
mod_sites_signature_motifs <- detect_signature_motifs(file_paths, config_file_path)

# Save the modified data frame with signature motifs
save(mod_sites_signature_motifs, file = file.path(store_dir, "significant_mod_sites_signature_motifs.RData"))

# Load conditions and effect direction from config.yaml
config <- yaml::read_yaml(config_file_path)

# Load configuration again to access signature motifs (already loaded above, but repeated for modularity)
signature_motif_pattern <- config$signature_motif_pattern

# Save separate tables for each signature motif with "YES" in the respective column
for (i in seq_along(signature_motif_pattern)) {
  # Get the signature name
  signature_name <- names(signature_motif_pattern)[i]
  
  # Filter rows where the signature motif is marked as "YES"
  temp_df <- mod_sites_signature_motifs[mod_sites_signature_motifs[[signature_name]] == "YES", ]
  
  # Save the filtered data frame for the current signature
  save(temp_df, file = file.path(store_dir, paste0(signature_name, "_significant_mod_sites.RData")))
}
