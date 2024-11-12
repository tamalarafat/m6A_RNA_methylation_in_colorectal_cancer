# Specify the projects dir containing all the functions of scRNA-seq analysis pipeline
source("~/Documents/Projects/my_projects/m6A_RNA_methylation_in_colorectal_cancer/library_and_packages/library_handler.R", .GlobalEnv)

# Functions to perfrom analysis tasks
funcset <- list.files("~/Documents/Projects/my_projects/nanopore_RNA_DNA_analysis_functions", pattern = "*.R$", full.names = TRUE)
sapply(funcset, source, .GlobalEnv)

# Directory to store the results
res_dir = "~/Documents/Projects/m6A_RNA_modification_project/RNAseq_analysis_output/"

store_folder = "01_gene_expression_tables"

# Create the output directory if it doesn't exist
if (!dir.exists(file.path(res_dir, store_folder))) {
  dir.create(file.path(res_dir, store_folder), recursive = TRUE, mode = "0777")
}

store_dir = file.path(res_dir, store_folder)

# Load conditions from config.yaml
config_file_path <- "~/Documents/Projects/my_projects/m6A_RNA_methylation_in_colorectal_cancer/library_and_packages/configuration.yaml"

# Data directroy
data_dir <- "~/Documents/Projects/m6A_RNA_modification_project/transcript_count_table/"

# Gene Transcript symbol table
gene_transcript_symbols <- "~/Documents/Projects/m6A_RNA_modification_project/modification_analysis_output/06_fasta_and_gene_transcript_table/gene_transcript_symbol.RData"

# Get the file paths
file_paths <- list.files(data_dir, ".tsv", full.names = TRUE)

# Give the sample names as input
sample_names = c("dld_ctr", "hct_ctr", "hct_ki", "dld_ki")

temp_paths <- "temporary_data" 

prepare_nanocount_data(file_paths, sample_names)

summarized_data <- generate_gene_count_table(temp_paths, config_file_path, gene_transcript_symbols)

save(summarized_data, file = file.path(store_dir, "gene_expression_table.RData"))

unlink(temp_paths, recursive = TRUE)
