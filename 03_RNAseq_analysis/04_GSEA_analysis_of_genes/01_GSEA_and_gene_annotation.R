# Specify the projects dir containing all the functions of scRNA-seq analysis pipeline
source("~/Documents/Projects/my_projects/m6A_RNA_methylation_in_colorectal_cancer/library_and_packages/library_handler.R", .GlobalEnv)

# Functions to perfrom analysis tasks
funcset <- list.files("~/Documents/Projects/my_projects/nanopore_RNA_DNA_analysis_functions", pattern = "*.R$", full.names = TRUE)
sapply(funcset, source, .GlobalEnv)

# Directory to store the results
res_dir = "~/Documents/Projects/m6A_RNA_modification_project/RNAseq_analysis_output/"

store_folder = "04_GSEA_condition_characteristic_genes/"

# Create the output directory if it doesn't exist
if (!dir.exists(file.path(res_dir, store_folder))) {
  dir.create(file.path(res_dir, store_folder), recursive = TRUE, mode = "0777")
}

store_dir = file.path(res_dir, store_folder)

# Load conditions from config.yaml
config_file_path <- "~/Documents/Projects/my_projects/m6A_RNA_methylation_in_colorectal_cancer/library_and_packages/configuration.yaml"

# Capture command-line arguments
file_paths <- "~/Documents/Projects/m6A_RNA_modification_project/RNAseq_analysis_output/03_condition_characteristic_genes/characteristic_deg_table.tsv"

# Gene Transcript symbol table
gene_transcript_symbols <- "~/Documents/Projects/m6A_RNA_modification_project/modification_analysis_output/06_fasta_and_gene_transcript_table/gene_transcript_symbol.RData"

# Annotate the genes table
gene_set_annotation <- annotate_genes(file_paths, config_file_path, gene_transcript_symbols)

# Get conditions gene set
gene_set <- gene_set_annotation[1]

# DEG table
deg_table <- gene_set_annotation[[2]]
write.table(deg_table, file = file.path(store_dir, "annotated_DE_genes_table.tsv"), row.names = FALSE, sep = "\t")

# ENTREZ ID
hg_entrez_universe <- gene_set_annotation[[3]]

# Get fasta subset and gene transcript table
perform_GSEA_DEG(gene_set, hg_entrez_universe, deg_table, store_dir)

