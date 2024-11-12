# Specify the projects dir containing all the functions of scRNA-seq analysis pipeline
source("~/Documents/Projects/my_projects/m6A_RNA_methylation_in_colorectal_cancer/library_and_packages/library_handler.R", .GlobalEnv)

# Functions to perfrom analysis tasks
funcset <- list.files("~/Documents/Projects/my_projects/nanopore_RNA_DNA_analysis_functions", pattern = "*.R$", full.names = TRUE)
sapply(funcset, source, .GlobalEnv)

# Directory to store the results
res_dir = "~/Documents/Projects/m6A_RNA_modification_project/RNAseq_analysis_output/"

store_folder = "05_cell_line_counts/"

# Create the output directory if it doesn't exist
if (!dir.exists(file.path(res_dir, store_folder))) {
  dir.create(file.path(res_dir, store_folder), recursive = TRUE, mode = "0777")
}

store_dir = file.path(res_dir, store_folder)

# Capture command-line arguments
file_paths <- "~/Documents/Projects/m6A_RNA_modification_project/RNAseq_analysis_output/01_gene_expression_tables/gene_expression_table.RData"

# Extract the counts for MPI cell line
hct_line = extract_cell_line_count(input_file = file_paths, cell_line_name = "hct", write_output = FALSE)

# Save the output
write.table(hct_line, file = file.path(store_dir, "gene_expression_hct_cell_line.tsv"), row.names = FALSE, sep = "\t")

# Get the count matrix
temp_counts = as.data.frame(hct_line$counts)

# Check for mettl3 expression
"ENSG00000165819" %in% rownames(temp_counts)

temp_counts["ENSG00000165819", ]

log2(temp_counts["ENSG00000165819", "hct_ctr"]/temp_counts["ENSG00000165819", "hct_ki"]) # 0.219678 log2

# Filter the data
# remove rows containing only 0's + remove rows if the gene expression is less than 10
keep = ((!rowSums(temp_counts != 0) == 0) & rowSums(temp_counts >= 10) == ncol(temp_counts))

# Subset the data
temp_counts <- temp_counts[keep, ]

temp_counts$log2fc = log2(temp_counts$hct_ctr/temp_counts$hct_ki)

temp_counts = temp_counts[abs(temp_counts$log2fc) > 0.25, ]

"ENSG00000165819" %in% rownames(temp_counts)

