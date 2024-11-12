# Specify the projects dir containing all the functions of scRNA-seq analysis pipeline
source("~/Documents/Projects/my_projects/m6A_RNA_methylation_in_colorectal_cancer/library_and_packages/library_handler.R", .GlobalEnv)

# Functions to perfrom analysis tasks
funcset <- list.files("~/Documents/Projects/my_projects/nanopore_RNA_DNA_analysis_functions", pattern = "*.R$", full.names = TRUE)
sapply(funcset, source, .GlobalEnv)

# Load the table containing gene information
genes_table = read.table("~/Documents/Projects/Git_repositories/effect_of_METTL3_on_m6A_in_cancer/Analyses/Analysis_of_all_samples/01_whole_data_analysis/02_analysis_on_counts/analysis_output/03_gene_count_tables/transcript_gene_table.csv")

# DEGs
res = loadRData("~/Documents/Projects/Git_repositories/effect_of_METTL3_on_m6A_in_cancer/Analyses/Analysis_of_all_samples/01_whole_data_analysis/02_analysis_on_counts/analysis_output/05_deg_files/DEG_DESeq_hct_ctr_hct_ki_mpi_ctr_mpi_ki.RData")

# Generate table with all DEG information
temp_df = data.frame(ID = rownames(res), log2FC = res$log2FoldChange, pVal = res$pvalue, adjpVal = res$padj)

temp_df$hgnc_symbol = genes_table$hgnc_symbol[match(temp_df$ID, genes_table$ensembl_gene_id)]

temp_df$entrezgene_id = genes_table$entrezgene_id[match(temp_df$ID, genes_table$ensembl_gene_id)]

# DEGs in control and DRACH
ctrl_deg = read.table("~/Documents/Projects/Git_repositories/effect_of_METTL3_on_m6A_in_cancer/Analyses/Analysis_of_all_samples/01_whole_data_analysis/02_analysis_on_counts/analysis_output/06_characteristic_degs/DEG_up_in_ctr_group.csv")

ctrl_deg$hgnc_symbol = genes_table$hgnc_symbol[match(ctrl_deg$ID, genes_table$ensembl_gene_id)]

ctrl_deg$entrezgene_id = genes_table$entrezgene_id[match(ctrl_deg$ID, genes_table$ensembl_gene_id)]

# Load the table with differential modification rate - filtered on all sites for p-value 0.05
temp_mod = loadRData("../Figure_3/all_significant_sites.RData")

# control group
temp_control = temp_mod[temp_mod$diff_mod_rate > 0 & temp_mod$DRACH == "YES", ]

# how many k-mers were detected
temp_control = temp_control[!is.na(temp_control$mers), ]

temp_DRACH_genes = intersect(ctrl_deg$ID, temp_control$ensembl_id)

#volcano plot
df <- temp_df %>%
  mutate(point_size = -log10(pVal),
         group = ifelse(log2FC > 0, "control", "ki"))

# order the df
df <- df[order(df$pVal, decreasing = FALSE),]

top_points <- df[which(df$ID %in% temp_DRACH_genes), ]

# Custom color palette for conditions
grp_col <- c("control" = "#00A08A", "ki" = "#F98400")  # Replace with your colors

p <- ggplot(df, aes(x = log2FC, y = point_size, color = group)) + 
  geom_point(aes(size = point_size), alpha = 0.6) +
  
  geom_point(data = top_points, aes(color = "#FD6467", size = 0.5), alpha = 1) + 
  
  # Repelling text labels for better separation
  geom_text_repel(data = top_points, aes(label = hgnc_symbol), 
                  size = 4,
                  box.padding = unit(0.5, "lines"),
                  point.padding = unit(0.1, "lines"),
                  nudge_y = 1.8,
                  force = 30) +
  # Labels for X and Y axes
  xlab("Fold change (log2 scale)") +
  ylab("p-Value (-log10)") +
  theme_classic() + 
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border = element_blank(),
    panel.background = element_blank(), # Background of the entire plot
    axis.line = element_blank(),
    axis.ticks.length = unit(0, "cm"), 
    axis.title = element_text(size = 18, face = "bold"),
    axis.text.x = element_text(size = 18, face = "bold"),
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

p

ggsave(filename = "05_Figure_5.png", plot = p, width = 8, height = 8, dpi = 300)



