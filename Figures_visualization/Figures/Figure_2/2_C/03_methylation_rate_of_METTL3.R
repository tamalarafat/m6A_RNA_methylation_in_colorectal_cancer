# Specify the projects dir containing all the functions of scRNA-seq analysis pipeline
source("~/Documents/Projects/my_projects/m6A_RNA_methylation_in_colorectal_cancer/library_and_packages/library_handler.R", .GlobalEnv)

# Functions to perfrom analysis tasks
funcset <- list.files("~/Documents/Projects/my_projects/nanopore_RNA_DNA_analysis_functions", pattern = "*.R$", full.names = TRUE)
sapply(funcset, source, .GlobalEnv)

# Data input directory
input_dir = "~/Documents/Projects/Git_repositories/effect_of_METTL3_on_m6A_in_cancer/Analyses/Analysis_of_all_samples/01_whole_data_analysis/analysis_output/03_filtered_merged_tables/"

# Directory to store the results
res_dir = getwd()

# Load the data mod table
temp_df = loadRData(str_c(input_dir, "filtered_conditions_across_cell_lines_all_possible_sites.RData"))

# Containing only the columns with modification rate data
df = temp_df[ , grep(pattern = "percent_modified", colnames(temp_df))]

# Change the column names to conditions
colnames(df) = str_replace(string = colnames(df), pattern = "percent_modified_", replacement = "")

df = df[, c(1,3,2,4)]

# gene transcript table:
temp_gt = loadRData("~/Documents/Projects/Git_repositories/effect_of_METTL3_on_m6A_in_cancer/Analyses/Analysis_of_all_samples/01_whole_data_analysis/02_analysis_on_counts/03_DET_analysis/gene_transcript_count.RData")

temp_gt$GENEID = str_replace(string = temp_gt$GENEID, pattern = "\\..*", replacement = "")

temp_trans = unique(temp_df$ref_seq_id)

temp_gene = temp_gt[grep(pattern = "ENSG00000165819", x = temp_gt$GENEID), ]

temp_df = temp_df[grep(pattern = str_c(temp_gene$TXNAME, collapse = "|"), x = temp_df$ref_seq_id), ]

# Containing only the columns with modification rate data
df = temp_df[ , grep(pattern = str_c(c("ref_seq_id", "start_pos", "percent_modified", "N_valid_cov"), collapse = "|"), colnames(temp_df))]

save(df, file = "METTL3_methylation.RData")

df = loadRData("METTL3_methylation.RData")

df = df[rowSums(is.na(df[, grep(pattern = "percent_modified", x = colnames(df))])) == 0, ]

length(unique(df$ref_seq_id))

# 3228 sites detected - 10 transcripts
df = df[(rowSums(df[, grep(pattern = "N_valid_cov", x = colnames(df))] > 30) != 0) & (rowSums(df[, grep(pattern = "percent_modified", x = colnames(df))] >= 10) != 0), ]

length(unique(df$ref_seq_id))

#311 after filtering - 4

# df = df[!rowSums(is.na(df)) == 4, ] # 4478 sites

dfv = df[, grep(pattern = "N_valid_cov", x = colnames(df))]

dfg = df[(rowSums(df[, grep(pattern = "N_valid_cov", x = colnames(df))] > 30) != 0) & (rowSums(df[, grep(pattern = "percent_modified", x = colnames(df))] >= 10) != 0), ]

mdf = dfg[, grep(pattern = str_c(c("ref_seq_id", "start_pos", "percent_modified"), collapse = "|"), colnames(dfg))]

mdf = mdf[rowSums(is.na(mdf[, grep(pattern = "percent_modified", x = colnames(mdf))])) == 0, ]

df = mdf

df = df[c("34516534", "4840348", "17426996", "33771052"), ]

rownames(df) <- NULL

df = df[,c(1, 3, 5, 4, 6)]

# Change the colNULL# Change the column names to conditions
colnames(df) = str_replace(string = colnames(df), pattern = "percent_modified_", replacement = "")

# Transform the data into long format for plotting
df_long = reshape2::melt(df, variable.name = "condition", value.name = "mod_rate")

# Change line types by groups (supp)
p <- ggplot(df_long, aes(x = condition, y = mod_rate, group = ref_seq_id)) +
  geom_line(aes(color = ref_seq_id), linewidth = 1.5) +
  geom_point(aes(color = ref_seq_id), size = 1.5) +
  labs(
    x = "Samples",
    y = "Modification rate",
    color = "METTL3 transcripts"
  ) +
  scale_x_discrete(labels = c("hct_ctrl" = "HCT control", "dld_ctrl" = "DLD control",
                              "hct_ki" = "HCT ki", "dld_ki" = "DLD ki")) +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border = element_blank(),
    axis.line = element_line(size = 0.5, colour = "#00A08A"),
    panel.background = element_blank(), # Background of the entire plot
    axis.title = element_text(size = 22, color = "black"),
    axis.ticks.length = unit(.20, "cm"), 
    axis.text = element_text(size = 24, colour = "black"),
    legend.position="top",
    legend.title = element_text(size = 22, color = "black"),
    legend.key.size = unit(2, "line"), 
    legend.text = element_text(size = 22)) + 
  guides(fill = guide_legend(title = "Transcripts"), color = guide_legend(override.aes = list(size = 4)))

ggsave(filename = "03_Figure_2.png", plot = p, width = 18, height = 10, dpi = 300)
