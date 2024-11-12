# Specify the projects dir containing all the functions of scRNA-seq analysis pipeline
projects_dir <- "~/Documents/Projects/Git_repositories/ontextractor/"

# Functions
funcset <- list.files(paste0(projects_dir, c("DataManipulator", "LibraryHandler", "TestStat", "analyses")), pattern = "*.R$", full.names = TRUE)
sapply(funcset, source, .GlobalEnv)

library("sleuth")
library(rhdf5)

# Data input directory
input_dir = "~/Documents/Projects/Git_repositories/effect_of_METTL3_on_m6A_in_cancer/Analyses/Analysis_of_all_samples/01_whole_data_analysis/02_analysis_on_counts/analysis_output/02_count_tables/"

# Directory to store the results
res_dir = "~/Documents/Projects/Git_repositories/effect_of_METTL3_on_m6A_in_cancer/Analyses/Analysis_of_all_samples/01_whole_data_analysis/02_analysis_on_counts/analysis_output/"

# List and load all the input files in the input directory
input_files = list.files(path = input_dir, pattern = ".tsv", full.names = TRUE)

# Replace occurrences of "//" with "/"
input_files <- gsub("//", "/", input_files)

# Assign names to the input files
names(input_files) <- str_replace(string = basename(input_files), pattern = "\\..*", replacement = "")

# Generate a sample table to track the coldata
coldata <- data.frame(sample = factor(names(input_files)), 
                      condition = factor(sub(pattern = ".*_", replacement = "", names(input_files))),
                      path = input_files)

print(coldata)

aa = read_kallisto(input_files[1])
bb = read_kallisto(input_files[2])
temp = aa$abundance

so <- sleuth_prep(coldata)

# Load a count table to get tge transcript ids
temp = read.table(input_files[1], header = TRUE)

h5write(temp, "abundance.h5","df")
h5ls("abundance.h5")

h5f = H5Fopen("abundance.h5")

h5f$df

s2c <- coldata

rownames(s2c) <- NULL

basename(dirname(filePath))




out_dir = getwd()

for (i in (c(1:length(input_files)))){
  if (!dir.exists(str_c(out_dir, "/", names(input_files)[i]))){
    dir.create(str_c(out_dir, "/", names(input_files)[i]), showWarnings = TRUE, recursive = FALSE, mode = "0777")
  }
  
  temp_out = str_c(out_dir, "/", names(input_files)[i], "/")
  
  # Load a count table to get tge transcript ids
  temp = read.table(input_files[i], header = TRUE)
  
  h5write(temp, file = str_c(temp_out, "abundance.h5"),"df")
}

temp_input_new = str_c(out_dir, "/", names(input_files))
names(temp_input_new) <- basename(temp_input_new)

temp_input_new

# Generate a sample table to track the coldata
ss <- data.frame(sample = names(temp_input_new), 
                 condition = sub(pattern = ".*_", replacement = "", names(temp_input_new)),
                 path = temp_input_new)

rownames(ss) <- NULL

so <- sleuth_prepNULLso <- sleuth_prep(ss, extra_bootstrap_summary = TRUE)

library(rhdf5)



file.path("..", "results", sample_id, "kallisto") 

so <- sleuth_prepNULLso <- sleuth_prep(s2c, extra_bootstrap_summary = TRUE, importer = read_tsv)

# List and load all the input files in the input directory
input_files = list.files(path = input_file, pattern = ".tsv", full.names = TRUE)

# Assign names to the input files
names(input_files) <- str_replace(string = basename(input_files), pattern = "\\..*", replacement = "")


if (!dir.exists(str_c(res_dir, "07_simulated_analysis_results"))){
  dir.create(str_c(res_dir, "07_simulated_analysis_results"), showWarnings = TRUE, recursive = FALSE, mode = "0777")
}

# Directory to store the results
res_dir = str_c(res_dir, "07_simulated_analysis_results", "/")

# Extract the mpi cell line and generate replicates; here two replicates. Storing the output for future references
temp_mpi_line = generate_simulated_replicates(input_file = temp_input, 
                                              cell_line = "mpi", 
                                              proportion_to_permutate = 0.05, 
                                              variable_std = 0.05, 
                                              num_replicates = 2, 
                                              write_output = TRUE, 
                                              store_dir = res_dir, 
                                              store_folder = "01_gene_count_tables")

# Input for DESeq analysis; directory of the gene count table
temp_input = str_c(res_dir, "01_gene_count_tables", "/", "mpi_ctr_mpi_ki_simulated_2_replicates_gene_counts.RData")

# Generate DESeq2 object and differentially expressed genes list
extract_deg_deseq(input_file = temp_input, 
                  reference_group = "ctr", 
                  control_for_cell_line = FALSE, 
                  filter_gene_exp_threshold = 10, 
                  write_output = TRUE, 
                  write_deseq_object = TRUE, 
                  store_dir = res_dir, 
                  store_folder = "02_deg_files")



# Directory of the gene count table
temp_input = str_c(res_dir, "02_deg_files", "/", "DEG_DESeq_rep_1_mpi_ctr_rep_1_mpi_ki_rep_2_mpi_ctr_rep_2_mpi_ki_rep_3_mpi_ctr_rep_3_mpi_ki.RData")

# Extract up-regulated and down-regulated genes corresponding to the control condition
extract_sample_characteristic_degs(input_file = temp_input, 
                                   filter_by_pval = TRUE, 
                                   filter_by_fc = TRUE, 
                                   return_sig_by_pval = TRUE, 
                                   return_sig_by_padj = FALSE, 
                                   log2fc_threshold = 0.5, 
                                   alpha_threshold = 0.05, 
                                   return_result_per_condition = TRUE, 
                                   return_geneIDs = TRUE, 
                                   write_output = TRUE, 
                                   store_dir = res_dir, 
                                   store_folder = "03_characteristic_degs")

