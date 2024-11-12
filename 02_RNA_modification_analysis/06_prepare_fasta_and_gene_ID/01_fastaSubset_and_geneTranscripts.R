# Specify the projects dir containing all the functions of scRNA-seq analysis pipeline
source("~/Documents/Projects/my_projects/m6A_RNA_methylation_in_colorectal_cancer/library_and_packages/library_handler.R", .GlobalEnv)

# Functions to perfrom analysis tasks
funcset <- list.files("~/Documents/Projects/my_projects/nanopore_RNA_DNA_analysis_functions", pattern = "*.R$", full.names = TRUE)
sapply(funcset, source, .GlobalEnv)

# Directory to store the results
res_dir = "~/Documents/Projects/m6A_RNA_modification_project/modification_analysis_output/"

store_folder = "06_fasta_and_gene_transcript_table"

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
fasta_path = "~/Documents/Projects/m6A_RNA_modification_project/ref/Homo_sapiens.GRCh38.cdna.all.fa"

# Get fasta subset and gene transcript table
fastaSub_geneTranscript <- fastaSubset_and_geneTransript_table(file_paths, config_file_path, fasta_path)

fasta_subset <- fastaSub_geneTranscript[[1]]

gene_transcript_table <- fastaSub_geneTranscript[[2]]

gene_transcript_table[] <- lapply(gene_transcript_table, function(x) {
  x <- ifelse(x == "", NA, x)  # Replace "" with NA
  type.convert(as.character(x), as.is = TRUE)  # Convert back to integer, factor, etc., if possible
})

# Write to fasta file containing only the reqfeq ID with sites that have a significant modification rate
writeXStringSet(fasta_subset, filepath = file.path(store_dir, "fasta_subset.fasta"))

# Save the gene transcript table
save(gene_transcript_table, file = file.path(store_dir, "gene_transcript_symbol.RData"))
