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
file_paths <- "~/Documents/Projects/m6A_RNA_modification_project/modification_analysis_output/07_modification_sites_and_motifs/significant_mod_sites_signature_motifs.RData"

# Gene Transcript symbol table
gene_transcript_symbols <- "~/Documents/Projects/m6A_RNA_modification_project/modification_analysis_output/06_fasta_and_gene_transcript_table/gene_transcript_symbol.RData"

# Annotate the transcript ids with gene ID symbols
annotated_mod_table <- annotate_modification_sites(file_paths, config_file_path, gene_transcript_symbols)

temp_df <- annotated_mod_table[[2]]

save(temp_df, file = file.path(store_dir, "annotated_mod_sites.Rdata"))
