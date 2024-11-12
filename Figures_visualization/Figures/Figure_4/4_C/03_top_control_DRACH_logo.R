# Specify the projects dir containing all the functions of scRNA-seq analysis pipeline
source("~/Documents/Projects/my_projects/m6A_RNA_methylation_in_colorectal_cancer/library_and_packages/library_handler.R", .GlobalEnv)

# Functions to perfrom analysis tasks
funcset <- list.files("~/Documents/Projects/my_projects/nanopore_RNA_DNA_analysis_functions", pattern = "*.R$", full.names = TRUE)
sapply(funcset, source, .GlobalEnv)

# Load the table with differential modification rate - filtered on all sites for p-value 0.05
temp_df = loadRData("../Figure_3/all_significant_sites.RData")

# control group
temp_control = temp_df[temp_df$diff_mod_rate > 0 & temp_df$DRACH == "YES", ]

# how many k-mers were detected
temp_control = temp_control[!is.na(temp_control$mers), ]

# Total 182 motifs
temp_motif = as.data.frame(sort(table(temp_control$mers), decreasing = TRUE))

temp_top = temp_control[temp_control$mers %in% temp_motif$Var1[1], ]

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
ggsave(filename = "03_Figure_4_control_DRACH_top.png", plot = p, width = 12, height = 6, dpi = 300)


temp_top = temp_control[temp_control$mers %in% temp_motif$Var1[2], ]

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
ggsave(filename = "03_Figure_4_control_DRACH_top2.png", plot = p, width = 12, height = 6, dpi = 300)


temp_top = temp_control[temp_control$mers %in% temp_motif$Var1[3], ]

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
ggsave(filename = "03_Figure_4_control_DRACH_top3.png", plot = p, width = 12, height = 6, dpi = 300)



temp_top = temp_control[temp_control$mers %in% temp_motif$Var1[4], ]

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
ggsave(filename = "03_Figure_4_control_DRACH_top4.png", plot = p, width = 12, height = 6, dpi = 300)


