# Specify the projects dir containing all the functions of scRNA-seq analysis pipeline
source("~/Documents/Projects/my_projects/m6A_RNA_methylation_in_colorectal_cancer/library_and_packages/library_handler.R", .GlobalEnv)

# Functions to perfrom analysis tasks
funcset <- list.files("~/Documents/Projects/my_projects/nanopore_RNA_DNA_analysis_functions", pattern = "*.R$", full.names = TRUE)
sapply(funcset, source, .GlobalEnv)

# Data input directory
input_dir = "~/Documents/Projects/Git_repositories/effect_of_METTL3_on_m6A_in_cancer/Analyses/Analysis_of_all_samples/01_whole_data_analysis/analysis_output/11_DRACH_motif/"

# Directory to store the results
res_dir = "~/Documents/Projects/Git_repositories/effect_of_METTL3_on_m6A_in_cancer/Analyses/Analysis_of_all_samples/01_whole_data_analysis/analysis_output/"

# Load the table with differential modification rate - filtered on all sites for p-value 0.05
temp_df = loadRData(str_c(input_dir, "DRACH_added_mod_table_curated_significant_joint_sites.RData"))

# Load the table containing gene information
genes_table = loadRData(str_c(res_dir, "/", "08_genes_and_their_function/curated_sig_joint_sites_gene_transcript_ids.RData"))

temp_df$ensembl_id = genes_table$ensembl_gene_id[match(temp_df$ref_seq_id, genes_table$ensembl_transcript_id_version)]

save(temp_df, file = "all_significant_sites.RData")

#volcano plot
df <- temp_df %>%
  mutate(point_size = -log10(pVal_ztest),
         diff_mod = ifelse(diff_mod_rate > 0, "control", "ki"))

# Limit the number of labeled points to the top 30 most significant
top_n <- 10

df <- df[order(df$pVal_ztest, decreasing = FALSE),]

top_points <- head(df, top_n)

# Custom color palette for conditions
grp_col <- c("control" = "#00A08A", "ki" = "#F98400")  # Replace with your colors

p <- ggplot(df, aes(x = diff_mod_rate, y = point_size, color = diff_mod)) + 
  geom_point(aes(size = point_size), alpha = 0.6) +
  
  # Repelling text labels for better separation
  geom_text_repel(data = top_points, aes(label = mers), nudge_x = -15, size = 5, force = 5) + 
  
  # Adding title with dynamic content
  ggtitle(paste(nrow(df), "sites;", length(unique(df$ref_seq_id)), 
                "transcripts;", length(unique(df$ensembl_id)), "genes", sep = " ")) +
  
  # Labels for X and Y axes
  xlab("Differential modification rate") +
  ylab("p-Value (-log10)") +
  theme_classic() + 
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border = element_blank(),
    panel.background = element_blank(), # Background of the entire plot
    axis.line = element_blank(),
    axis.ticks.length = unit(0, "cm"), 
    axis.title = element_text(size = 18),
    axis.text.x = element_text(size = 18),
    axis.text.y = element_blank(),
    plot.title = element_text(size = 18),
    legend.position = "top",   # Move legend to the top
    legend.title = element_text(size = 18, color = "black"),
    legend.key.size = unit(2, "line"),
    legend.text = element_text(size = 18, color = "black")
  ) +
  
  # Custom manual color scale using grp_col
  scale_color_manual(values = grp_col, labels = c("control" = "Control", "ki" = "Knock-down")) +
  
  # Customizing the legend
  guides(
    color = guide_legend(title = "Conditions", override.aes = list(size = 4)),  # Title for conditions
    size = "none"  # Disable size legend
  )

ggsave(filename = "03_Figure_3.png", plot = p, width = 8, height = 8, dpi = 300)


