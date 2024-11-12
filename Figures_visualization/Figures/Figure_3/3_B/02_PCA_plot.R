# Specify the projects dir containing all the functions of scRNA-seq analysis pipeline
source("~/Documents/Projects/my_projects/m6A_RNA_methylation_in_colorectal_cancer/library_and_packages/library_handler.R", .GlobalEnv)

# Functions to perfrom analysis tasks
funcset <- list.files("~/Documents/Projects/my_projects/nanopore_RNA_DNA_analysis_functions", pattern = "*.R$", full.names = TRUE)
sapply(funcset, source, .GlobalEnv)

# Data input directory
input_dir = "~/Documents/Projects/Git_repositories/effect_of_METTL3_on_m6A_in_cancer/Analyses/Analysis_of_all_samples/01_whole_data_analysis/analysis_output/08_genes_and_their_function/"

# Directory to store the results
res_dir = "~/Documents/Projects/Git_repositories/effect_of_METTL3_on_m6A_in_cancer/Analyses/Analysis_of_all_samples/01_whole_data_analysis/analysis_output/"

# Load the table with differential modification rate - filtered on all sites for p-value 0.05
temp_df = loadRData(str_c("~/Documents/Projects/Git_repositories/effect_of_METTL3_on_m6A_in_cancer/Analyses/Analysis_of_all_samples/01_whole_data_analysis/analysis_output/07_final_modification_tables_with_added_motif/", "motif_added_mod_table_curated_significant_joint_sites.RData"))

# Containing only the columns with modification rate data
df = temp_df[ , grep(pattern = "percent_modified", colnames(temp_df))]

# Change the column names to conditions
colnames(df) = str_replace(string = colnames(df), pattern = "percent_modified_", replacement = "")

colnames(df)

# df = df[rowSums(df != 0) == 4, ]

res <- prcomp(t(df), center = TRUE)

sample_colors <- c("hct_ctrl" = "#00A08A", "dld_ctrl" =  "#046C9A","hct_ki" = "#FD6467", "dld_ki" = "#F98400")

plot(
  res$x[,1], res$x[,2],
  col = sample_colors,
  pch = 19,
  main = "Methylation rate across samples",
  xlab = paste("PC1 (", round(100 * summary(res)$importance[2,1]), "%)", sep = ""),
  ylab = paste("PC2 (", round(100 * summary(res)$importance[2,2]), "%)", sep = ""), 
  cex = text(res$x[,1],res$x[,2],labels=c("hct_ctrl" = "HCT control", "dld_ctrl" = "DLD control",
                                          "hct_ki" = "HCT ki", "dld_ki" = "DLD ki"), offset = c(1, 1), pos = 4, cex = 1), cex.axis = 1.5, cex.lab = 1.5, cex.main = 1.5)


# Perform PCA
res <- prcomp(t(df), center = TRUE, scale. = TRUE)

# Define sample colors
sample_colors <- c("hct_ctrl" = "#00A08A", "dld_ctrl" =  "#046C9A","hct_ki" = "#FD6467", "dld_ki" = "#F98400")

# Determine the range for the axes, expanded by a certain factor
x_range <- range(res$x[,1]) * 1.5  # Expand by 10%
y_range <- range(res$x[,2]) * 1.5  # Expand by 10%


# 1. Open jpeg file
jpeg("02_Figure_3.jpg", width = "2400", height = "2400", bg = "white", res = 300)

# 2. Create the plot
# Plot with expanded axes and larger points
plot(
  res$x[,1], res$x[,2],
  col = sample_colors,
  pch = 19,  # Shape of the points
  main = "Methylation rate across samples",
  xlab = paste("PC1 (", round(100 * summary(res)$importance[2,1]), "%)", sep = ""),
  ylab = paste("PC2 (", round(100 * summary(res)$importance[2,2]), "%)", sep = ""), 
  xlim = x_range,  # Set expanded x-axis limits
  ylim = y_range,  # Set expanded y-axis limits
  cex = 2,         # Increase the size of the dots
  cex.axis = 1.5,  # Increase the size of the axis labels
  cex.lab = 1.5,   # Increase the size of the axis titles
  cex.main = 1.5   # Increase the size of the main title
)

# Add text labels
text(
  res$x[,1], res$x[,2],
  labels = c("hct_ctrl" = "HCT control", "dld_ctrl" = "DLD control",
             "hct_ki" = "HCT ki", "dld_ki" = "DLD ki"),
  pos = 3,    # Position of the text (4 means to the right of the point)
  offset = 1, # Distance of the text from the point
  cex = 1   # Size of the text labels
)

# 3. Close the file
dev.off()
