# Specify the projects dir containing all the functions of scRNA-seq analysis pipeline
source("~/Documents/Projects/my_projects/m6A_RNA_methylation_in_colorectal_cancer/library_and_packages/library_handler.R", .GlobalEnv)

# Functions to perfrom analysis tasks
funcset <- list.files("~/Documents/Projects/my_projects/nanopore_RNA_DNA_analysis_functions", pattern = "*.R$", full.names = TRUE)
sapply(funcset, source, .GlobalEnv)

# Data input directory
input_dir = "~/Documents/Projects/Git_repositories/effect_of_METTL3_on_m6A_in_cancer/Analyses/Analysis_of_all_samples/01_whole_data_analysis/analysis_output/12_Annotated_tables/"

# Directory to store the results
res_dir = "~/Documents/Projects/Git_repositories/effect_of_METTL3_on_m6A_in_cancer/Analyses/Analysis_of_all_samples/01_whole_data_analysis/analysis_output/"

# Load the table with differential modification rate - filtered on all sites for p-value 0.05
temp_df = loadRData("../Figure_3/all_significant_sites.RData")

# control group
temp_control = temp_df[temp_df$diff_mod_rate > 0 & temp_df$DRACH == "YES", ]

# Your vector of sequences
DRACH <-  temp_control$mers

# Convert the vector into a matrix where each column is a position in the sequences
seq_matrix <- do.call(rbind, strsplit(DRACH, ""))

# Calculate the proportion of each base at each position
proportions <- apply(seq_matrix, 2, function(x) {
  table(factor(x, levels = c("A", "C", "G", "U"))) / length(x)
})

# Convert the matrix of proportions to a data frame
proportion_df <- as.data.frame(t(proportions))

proportion_df$position = as.factor(rownames(proportion_df))

# Transform the data into long format for plotting
df_long = reshape2::melt(proportion_df, variable.name = "position", value.name = "proportion")

# Rename columns to represent positions
colnames(df_long)[2] <- "bases"

# Display the result
print(proportion_df)

p <- ggplot(data = df_long, aes(x = position, y = proportion)) + 
  geom_bar(stat = "identity", position = "fill", aes(fill = bases)) +
  xlab("Position") + ylab("Proportion of bases") + 
  scale_x_discrete(labels = c("1" = "D", "2" = "R", "3" = "A", "4" = "C", "5" = "H")) +
  scale_fill_manual(values = c("A" = grp_col[2], "U" = grp_col[1], "G" = grp_col[3], "C" = grp_col[4])) + 
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border = element_blank(),
    #axis.line = element_line(color = grp_col[2], size = 2),
    panel.background = element_blank(), # Background of the entire plot
    axis.title = element_text(size = 22, color = "black"),
    axis.ticks.length = unit(0, "cm"), 
    axis.text = element_text(size = 22, colour = "black"),
    legend.position = "top",
    legend.title = element_blank(),
    legend.key.size = unit(1, "line"), 
    legend.text = element_text(size = 22)) + 
  guides(color = guide_legend(override.aes = list(size = 4)))

ggsave(filename = "02_Figure_4_control_DRACH_proportion.png", plot = p, width = 12, height = 6, dpi = 300)
