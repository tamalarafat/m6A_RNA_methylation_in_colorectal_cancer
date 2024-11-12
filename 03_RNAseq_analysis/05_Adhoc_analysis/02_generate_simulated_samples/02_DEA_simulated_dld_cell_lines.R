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


# Load conditions from config.yaml
config_file_path <- "~/Documents/Projects/my_projects/m6A_RNA_methylation_in_colorectal_cancer/library_and_packages/configuration.yaml"

# file path
file_paths <- "~/Documents/Projects/m6A_RNA_modification_project/RNAseq_analysis_output/06_simulated_gene_expression_tables/simulated_replicates_dld_cell_line.RData"

# Perform DESeq analysis
dds_res <- extract_simulated_deg_deseq(input_file = file_paths, 
                  reference_group = "control", 
                  control_for_cell_line = FALSE, 
                  filter_gene_exp_threshold = 10, 
                  write_output = FALSE)

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

# Extract up-regulated and down-regulated genes corresponding to the control condition
temp_characteristic_genes <- extract_simulated_characteristic_degs(input_file = deseq_res, 
                                   filter_by_pval = TRUE, 
                                   filter_by_fc = TRUE, 
                                   return_sig_by_pval = TRUE, 
                                   return_sig_by_padj = FALSE, 
                                   log2fc_threshold = 0.5, 
                                   alpha_threshold = 0.05, 
                                   return_result_per_condition = TRUE, 
                                   return_geneIDs = TRUE, 
                                   write_output = FALSE)

characteristic_genes <- temp_characteristic_genes[[3]]
# write the table in tsv format
write.table(characteristic_genes, file = file.path(store_dir, "characteristic_deg_table.tsv"), row.names = FALSE, sep = "\t")


write.table(temp_characteristic_genes[[1]], file = paste0(store_dir, "/", names(temp_characteristic_genes)[1], "_condition_characteristic_genes.tsv"), row.names = FALSE, sep = "\t")

write.table(temp_characteristic_genes[[2]], file = paste0(store_dir, "/", names(temp_characteristic_genes)[2], "_condition_characteristic_genes.tsv"), row.names = FALSE, sep = "\t")

