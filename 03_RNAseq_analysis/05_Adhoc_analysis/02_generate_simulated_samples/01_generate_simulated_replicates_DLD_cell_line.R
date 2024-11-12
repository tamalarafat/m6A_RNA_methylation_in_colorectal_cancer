# Specify the projects dir containing all the functions of scRNA-seq analysis pipeline
source("~/Documents/Projects/my_projects/m6A_RNA_methylation_in_colorectal_cancer/library_and_packages/library_handler.R", .GlobalEnv)

# Functions to perfrom analysis tasks
funcset <- list.files("~/Documents/Projects/my_projects/nanopore_RNA_DNA_analysis_functions", pattern = "*.R$", full.names = TRUE)
sapply(funcset, source, .GlobalEnv)

# Directory to store the results
res_dir = "~/Documents/Projects/m6A_RNA_modification_project/RNAseq_analysis_output/"

store_folder = "06_simulated_gene_expression_tables/"

# Create the output directory if it doesn't exist
if (!dir.exists(file.path(res_dir, store_folder))) {
  dir.create(file.path(res_dir, store_folder), recursive = TRUE, mode = "0777")
}

store_dir = file.path(res_dir, store_folder)

# Capture command-line arguments
file_paths <- "~/Documents/Projects/m6A_RNA_modification_project/RNAseq_analysis_output/01_gene_expression_tables/gene_expression_table.RData"

# Extract the mpi cell line and generate replicates; here two replicates. Storing the output for future references
temp_dld_line = generate_simulated_replicates(input_file = file_paths, 
                                              cell_line = "dld", 
                                              proportion_to_permutate = 0.05, 
                                              variable_std = 0.05, 
                                              num_replicates = 2, 
                                              write_output = FALSE)


temp_count <- as.data.frame(temp_dld_line$counts)

# Save the output
write.table(temp_count, file = file.path(store_dir, "simulated_gene_expression_dld_cell_line.tsv"), row.names = FALSE, sep = "\t")

# Save the output
save(temp_dld_line, file = file.path(store_dir, "simulated_replicates_dld_cell_line.RData"))
