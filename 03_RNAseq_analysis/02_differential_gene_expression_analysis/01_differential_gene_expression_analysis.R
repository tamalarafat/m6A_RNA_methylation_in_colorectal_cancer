# Specify the projects dir containing all the functions of scRNA-seq analysis pipeline
source("~/Documents/Projects/my_projects/m6A_RNA_methylation_in_colorectal_cancer/library_and_packages/library_handler.R", .GlobalEnv)

# Functions to perfrom analysis tasks
funcset <- list.files("~/Documents/Projects/my_projects/nanopore_RNA_DNA_analysis_functions", pattern = "*.R$", full.names = TRUE)
sapply(funcset, source, .GlobalEnv)

# Directory to store the results
res_dir = "~/Documents/Projects/m6A_RNA_modification_project/RNAseq_analysis_output/"

store_folder = "02_differential_gene_expression_output"

# Create the output directory if it doesn't exist
if (!dir.exists(file.path(res_dir, store_folder))) {
  dir.create(file.path(res_dir, store_folder), recursive = TRUE, mode = "0777")
}

store_dir = file.path(res_dir, store_folder)

# Load conditions from config.yaml
config_file_path <- "~/Documents/Projects/my_projects/m6A_RNA_methylation_in_colorectal_cancer/library_and_packages/configuration.yaml"

# Capture command-line arguments
file_paths <- "~/Documents/Projects/m6A_RNA_modification_project/RNAseq_analysis_output/01_gene_expression_tables/gene_expression_table.RData"

dds_res <- perform_DEA_deseq(file_paths, config_file_path)

# get and store deseq2 object
deseq_object <- dds_res$dds
save(deseq_object, file = file.path(store_dir, "DESeq_dataset.RData"))

# get and store deseq2 result
deseq_res <- dds_res$res
save(deseq_res, file = file.path(store_dir, "DESeq_results.RData"))

# get and store deg table
deg_table <- extract_DEG_table(deseq_res)
# write the table in tsv format
write.table(deg_table, file = file.path(store_dir, "full_deg_table.tsv"), row.names = FALSE, sep = "\t")
