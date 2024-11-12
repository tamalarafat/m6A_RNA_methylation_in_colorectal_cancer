# Specify the projects dir containing all the functions of scRNA-seq analysis pipeline
source("~/Documents/Projects/my_projects/m6A_RNA_methylation_in_colorectal_cancer/library_and_packages/library_handler.R", .GlobalEnv)

# Functions to perfrom analysis tasks
funcset <- list.files("~/Documents/Projects/my_projects/nanopore_RNA_DNA_analysis_functions", pattern = "*.R$", full.names = TRUE)
sapply(funcset, source, .GlobalEnv)

# Data input directory
input_dir = "~/Documents/Projects/Git_repositories/effect_of_METTL3_on_m6A_in_cancer/Analyses/Analysis_of_all_samples/01_whole_data_analysis/analysis_output/04_filtered_merged_tables_joint_sites/"

# Directory to store the results
res_dir = getwd()

# Load the data mod table
temp_df = loadRData(str_c(input_dir, "filtered_conditions_across_cell_lines_joint_sites.RData"))

# Containing only the columns with modification rate data
df = temp_df[ , grep(pattern = "percent_modified", colnames(temp_df))]

# Change the column names to conditions
colnames(df) = str_replace(string = colnames(df), pattern = "percent_modified_", replacement = "")

df = df[, c(1,3,2,4)]

# Transform the data into long format for plotting
df_long = reshape2::melt(df, variable.name = "condition", value.name = "mod_rate")

# gene transcript table:
temp_gt = loadRData("~/Documents/Projects/Git_repositories/effect_of_METTL3_on_m6A_in_cancer/Analyses/Analysis_of_all_samples/01_whole_data_analysis/02_analysis_on_counts/03_DET_analysis/gene_transcript_count.RData")

temp_trans = unique(temp_df$ref_seq_id)

temp_gene = temp_gt[match(temp_trans, temp_gt$TXNAME), ]

# Boxplot of differential modification rate for all the shared sites across the samples of two conditions
p <- ggplot(df_long, aes(x = condition, y = mod_rate)) +
  geom_boxplot(coef = 6) +
  #geom_dotplot(binaxis='y', stackdir='center', dotsize = 0.4, binwidth = 1/40) + 
  #geom_jitter(shape = 16, position = position_jitter(0.1)) +
  theme_classic() +
  labs(
    title = str_c(nrow(df), "sites;", length(unique(temp_df$ref_seq_id)), "transcripts;", length(unique(temp_gene$GENEID)), "genes", sep = " "),
    x = "Samples",
    y = "Modification Rate"
  ) +
  scale_x_discrete(labels = c("hct_ctrl" = "HCT control", "dld_ctrl" = "DLD control",
                            "hct_ki" = "HCT ki", "dld_ki" = "DLD ki")) + 
  theme(
    line = element_blank(),
    panel.border = element_blank(),
    axis.title = element_text(size = 18),
    axis.text = element_text(size = 18),
    title = element_text(size = 18))

ggsave(filename = "01_Figure_2.png", plot = p, width = 12, height = 6, dpi = 300)

