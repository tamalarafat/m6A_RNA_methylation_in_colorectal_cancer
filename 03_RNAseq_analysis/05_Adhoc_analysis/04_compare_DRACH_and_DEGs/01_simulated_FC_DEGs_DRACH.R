# Specify the projects dir containing all the functions of scRNA-seq analysis pipeline
projects_dir <- "~/Documents/Projects/Git_repositories/ontextractor/"

# Functions
funcset <- list.files(paste0(projects_dir, c("DataManipulator", "LibraryHandler", "TestStat", "analyses")), pattern = "*.R$", full.names = TRUE)
sapply(funcset, source, .GlobalEnv)

# Clean the query cache of biomaRt
biomaRt::biomartCacheClear()

# Data input directory
input_dir = "~/Documents/Projects/Git_repositories/effect_of_METTL3_on_m6A_in_cancer/Analyses/Analysis_of_all_samples/01_whole_data_analysis/analysis_output/"

# Directory to store the results
res_dir = "~/Documents/Projects/Git_repositories/effect_of_METTL3_on_m6A_in_cancer/Analyses/Analysis_of_all_samples/01_whole_data_analysis/analysis_output/"

# Background gene set tested for the modification analysis - this is the filtered set of transcripts used for t-test
temp_backgroud = loadRData(str_c("~/Documents/Projects/Git_repositories/effect_of_METTL3_on_m6A_in_cancer/Analyses/Analysis_of_all_samples/01_whole_data_analysis/analysis_output/", "05_final_modification_tables/", "filtered_mod_table_curated_joint_sites.RData"))

# get the gene annotation file
df_gene = read.table("~/Documents/Projects/Git_repositories/effect_of_METTL3_on_m6A_in_cancer/Analyses/Analysis_of_all_samples/01_whole_data_analysis/02_analysis_on_counts/analysis_output/03_gene_count_tables/transcript_gene_table.csv")

background_gene = df_gene[match(temp_backgroud$ref_seq_id, df_gene$ensembl_transcript_id_version), ]

# Number of transcripts
length(unique(temp_backgroud$ref_seq_id))

# Number of genes
length(unique(background_gene$ensembl_gene_id))

temp_backgenes = unique(background_gene$ensembl_gene_id)



# DRACH table Load the table with differential modification rate - filtered on all sites for p-value 0.05
temp_drach = loadRData(str_c(input_dir, "12_Annotated_tables/geneID_annotated_DRACH_table_curated_significant_joint_sites.RData"))

# Total transcripts
length(unique(temp_drach$ref_seq_id))

# Contol
temp_ctrl = temp_drach[temp_drach$diff_mod_rate >= 0, ]

# Number of transcripts
length(unique(temp_ctrl$ref_seq_id))

ctrl_gene = df_gene[match(temp_ctrl$ref_seq_id, df_gene$ensembl_transcript_id_version), ]

# Number of genes
length(unique(ctrl_gene$ensembl_gene_id))

# ctrl genes 
drach_ctrl_genes = unique(ctrl_gene$ensembl_gene_id)


# DEG file
deg_dir = "~/Documents/Projects/Git_repositories/effect_of_METTL3_on_m6A_in_cancer/Analyses/Analysis_of_all_samples/01_whole_data_analysis/02_analysis_on_counts/analysis_output/07_simulated_analysis_results/03_characteristic_degs"

# load the degs identified in the control group
crtl_degs = read.table(str_c(deg_dir, "/", "DEG_up_in_ctr_group.csv"))

crtl_degs = crtl_degs[crtl_degs$log2FC >= 0.25, ]

temp_ctrl_degs = crtl_degs$ID

"ENSG00000165819" %in% temp_ctrl_degs

dd = crtl_degs[crtl_degs$padj < 0.05, ]

temp_ctrl_degs = dd$ID

"ENSG00000165819" %in% temp_ctrl_degs

# Common between control genes and DRACH genes
temp_ctrl_common = intersect(drach_ctrl_genes, temp_ctrl_degs)

ctrl_df = data.frame("up" = c(length(drach_ctrl_genes) - length(temp_ctrl_common), length(temp_ctrl_common)), 
                     "nonUp" = c((length(temp_backgenes) - length(drach_ctrl_genes)) - (length(temp_ctrl_degs) - length(temp_ctrl_common)), (length(temp_ctrl_degs) - length(temp_ctrl_common))), 
                     row.names = c("DRACH", "DEGs")) # 182 out of 783

fisher.test(ctrl_df) # p-value = 1.743e-06

fileGenerator(temp_ctrl_common, fileName = "DRACH_DEG_CTRL_Common_Simulated_MPI_padj_05.txt")

# KI
temp_ki = temp_drach[temp_drach$diff_mod_rate < 0, ]

# Number of transcripts
length(unique(temp_ki$ref_seq_id))

ki_gene = df_gene[match(temp_ki$ref_seq_id, df_gene$ensembl_transcript_id_version), ]

# Number of genes
length(unique(ki_gene$ensembl_gene_id))

# ctrl genes 
drach_ki_genes = unique(ki_gene$ensembl_gene_id)

# load the degs identified in the control group
ki_degs = read.table(str_c(deg_dir, "/", "DEG_up_in_ki_group.csv"))

ki_degs = ki_degs[abs(ki_degs$log2FC) >= 0.25, ]

temp_ki_degs = ki_degs$ID

# Common between control genes and DRACH genes
temp_ki_common = intersect(drach_ki_genes, temp_ki_degs) # No overlap was found

ki_df = data.frame("up" = c(length(drach_ki_genes) - length(temp_ki_common), length(temp_ki_common)), 
                     "nonUp" = c((length(temp_backgenes) - length(drach_ki_genes)) - (length(temp_ki_degs) - length(temp_ki_common)), (length(temp_ki_degs) - length(temp_ki_common))), 
                     row.names = c("DRACH", "DEGs"))

fisher.test(ki_df) # p-value = 0.08291

fileGenerator(temp_ki_common, fileName = "DRACH_DEG_KI_Common_Simulated_MPI.txt")
