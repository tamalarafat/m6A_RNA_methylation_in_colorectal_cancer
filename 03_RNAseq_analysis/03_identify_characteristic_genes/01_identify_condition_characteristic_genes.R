# Specify the projects dir containing all the functions of scRNA-seq analysis pipeline
source("~/Documents/Projects/my_projects/m6A_RNA_methylation_in_colorectal_cancer/library_and_packages/library_handler.R", .GlobalEnv)

# Functions to perfrom analysis tasks
funcset <- list.files("~/Documents/Projects/my_projects/nanopore_RNA_DNA_analysis_functions", pattern = "*.R$", full.names = TRUE)
sapply(funcset, source, .GlobalEnv)

# Directory to store the results
res_dir = "~/Documents/Projects/m6A_RNA_modification_project/RNAseq_analysis_output/"

store_folder = "03_condition_characteristic_genes"

# Create the output directory if it doesn't exist
if (!dir.exists(file.path(res_dir, store_folder))) {
  dir.create(file.path(res_dir, store_folder), recursive = TRUE, mode = "0777")
}

store_dir = file.path(res_dir, store_folder)

# Load conditions from config.yaml
config_file_path <- "~/Documents/Projects/my_projects/m6A_RNA_methylation_in_colorectal_cancer/library_and_packages/configuration.yaml"

# Capture command-line arguments
file_paths <- "~/Documents/Projects/m6A_RNA_modification_project/RNAseq_analysis_output/02_differential_gene_expression_output/DESeq_results.RData"

characteristic_degs <- identify_characteristic_degs(file_paths, config_file_path)

characteristic_genes <- characteristic_degs[[3]]
# write the table in tsv format
write.table(characteristic_genes, file = file.path(store_dir, "characteristic_deg_table.tsv"), row.names = FALSE, sep = "\t")


write.table(characteristic_degs[[1]], file = paste0(store_dir, "/", names(characteristic_degs)[1], "_condition_characteristic_genes.tsv"), row.names = FALSE, sep = "\t")

write.table(characteristic_degs[[2]], file = paste0(store_dir, "/", names(characteristic_degs)[2], "_condition_characteristic_genes.tsv"), row.names = FALSE, sep = "\t")

