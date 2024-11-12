# Specify the projects dir containing all the functions of scRNA-seq analysis pipeline
projects_dir <- "~/Documents/Projects/Git_repositories/ontextractor/"

# Functions
funcset <- list.files(paste0(projects_dir, c("DataManipulator", "LibraryHandler", "TestStat", "analyses")), pattern = "*.R$", full.names = TRUE)
sapply(funcset, source, .GlobalEnv)

# Data input directory
input_dir = "~/Documents/Projects/Git_repositories/effect_of_METTL3_on_m6A_in_cancer/Analyses/Analysis_of_all_samples/01_whole_data_analysis/02_analysis_on_counts/analysis_output/05_deg_files/"

# Directory to store the results
res_dir = "~/Documents/Projects/Git_repositories/effect_of_METTL3_on_m6A_in_cancer/Analyses/Analysis_of_all_samples/01_whole_data_analysis/02_analysis_on_counts/analysis_output/"

# Input for the deseq results
temp_res = str_c(input_dir, "DEG_DESeq_hct_ctr_hct_ki_mpi_ctr_mpi_ki.RData")

# load deseq results
res = loadRData(temp_res)

# create dtaframe with all the genes and their information
temp_df = data.frame(ID = rownames(res), log2FC = res$log2FoldChange, pVal = res$pvalue, adjpVal = res$padj)

# save the data frame
save(temp_df, file = str_c(input_dir, "all_genes_FC.RData"))

# input for the deseq object
temp_dds = str_c(input_dir, "DESeq_object_hct_ctr_hct_ki_mpi_ctr_mpi_ki.RData")

# the DESeq object
dds = loadRData(temp_dds)

# DEGs for the control group - significant ones
ctrldeg = read.table(str_c(res_dir, "06_characteristic_degs/DEG_up_in_ctr_group.csv"))

# DEGs for KI group - significant ones
kideg = read.table(str_c(res_dir, "06_characteristic_degs/DEG_up_in_ki_group.csv"))

## DRACH

# Load all DRACH gens

# load control group DRACH geens

# Load KI group DRACH genes

df = data.frame("peak" = c(9336, 20), "nopeak" = c(18401, 16), row.names = c("chip", "markers")) # 20 out of 23
