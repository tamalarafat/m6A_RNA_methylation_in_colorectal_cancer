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
temp_df = loadRData("all_significant_sites.RData")

# control group
temp_control = temp_df[temp_df$diff_mod_rate > 0, ]

# how many k-mers were detected
temp_control = temp_control[!is.na(temp_control$mers), ]

# Total 182 motifs
temp_motif = as.data.frame(sort(table(temp_control$mers), decreasing = TRUE))

temp_top = temp_control[temp_control$mers %in% "GGACU", ]

p <- ggseqlogo(temp_top$mers, method = "probability") +
  labs(title = paste(nrow(temp_top), "sites;", length(unique(temp_top$ref_seq_id)), 
                     "transcripts;", length(unique(temp_top$ensembl_id)), "genes", sep = " ")) + 
  theme_classic() + 
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border = element_blank(),
    panel.background = element_blank(), # Background of the entire plot
    axis.line = element_blank(),
    axis.ticks.length = unit(0, "cm"), 
    axis.title = element_blank(),
    axis.text.x = element_blank(),
    axis.text.y = element_blank(),
    plot.title = element_text(size = 18),
    legend.position = "top",   # Move legend to the top
    legend.title = element_text(size = 18, color = "black"),
    legend.key.size = unit(2, "line"),
    legend.text = element_text(size = 18, color = "black")
  )
ggsave(filename = "06_Figure_3_control_motif_top.png", plot = p, width = 12, height = 6, dpi = 300)

# control group
temp_ki = temp_df[temp_df$diff_mod_rate < 0, ]

# how many k-mers were detected
temp_ki = temp_ki[!is.na(temp_ki$mers), ]

# Total 182 motifs
temp_motif = as.data.frame(sort(table(temp_ki$mers), decreasing = TRUE))

temp_top = temp_ki[temp_ki$mers %in% "CCAUG", ]

p <- ggseqlogo(temp_top$mers, method = "probability") +
  labs(title = paste(nrow(temp_top), "sites;", length(unique(temp_top$ref_seq_id)), 
                     "transcripts;", length(unique(temp_top$ensembl_id)), "genes", sep = " ")) + 
  theme_classic() + 
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border = element_blank(),
    panel.background = element_blank(), # Background of the entire plot
    axis.line = element_blank(),
    axis.ticks.length = unit(0, "cm"), 
    axis.title = element_blank(),
    axis.text.x = element_blank(),
    axis.text.y = element_blank(),
    plot.title = element_text(size = 18),
    legend.position = "top",   # Move legend to the top
    legend.title = element_text(size = 18, color = "black"),
    legend.key.size = unit(2, "line"),
    legend.text = element_text(size = 18, color = "black")
  )
ggsave(filename = "06_Figure_3_ki_motif_top.png", plot = p, width = 12, height = 6, dpi = 300)
