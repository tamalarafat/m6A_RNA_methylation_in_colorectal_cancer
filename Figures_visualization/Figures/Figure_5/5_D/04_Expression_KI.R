# Specify the projects dir containing all the functions of scRNA-seq analysis pipeline
source("~/Documents/Projects/my_projects/m6A_RNA_methylation_in_colorectal_cancer/library_and_packages/library_handler.R", .GlobalEnv)

# Functions to perfrom analysis tasks
funcset <- list.files("~/Documents/Projects/my_projects/nanopore_RNA_DNA_analysis_functions", pattern = "*.R$", full.names = TRUE)
sapply(funcset, source, .GlobalEnv)

# Load the table with differential modification rate - filtered on all sites for p-value 0.05
dds = loadRData("~/Documents/Projects/Git_repositories/effect_of_METTL3_on_m6A_in_cancer/Analyses/Analysis_of_all_samples/01_whole_data_analysis/02_analysis_on_counts/analysis_output/05_deg_files/DESeq_object_hct_ctr_hct_ki_mpi_ctr_mpi_ki.RData")

# DEGs
DEGs = read.table("~/Documents/Projects/Git_repositories/effect_of_METTL3_on_m6A_in_cancer/Analyses/Analysis_of_all_samples/01_whole_data_analysis/02_analysis_on_counts/analysis_output/06_characteristic_degs/DEG_up_in_ki_group.csv")

DEGs = DEGs[order(DEGs$log2FC, decreasing = FALSE), ]

# Load the table containing gene information
genes_table = read.table("~/Documents/Projects/Git_repositories/effect_of_METTL3_on_m6A_in_cancer/Analyses/Analysis_of_all_samples/01_whole_data_analysis/02_analysis_on_counts/analysis_output/03_gene_count_tables/transcript_gene_table.csv")

DEGs$hgnc_symbol = genes_table$hgnc_symbol[match(DEGs$ID, genes_table$ensembl_gene_id)]

DEGs$entrezgene_id = genes_table$entrezgene_id[match(DEGs$ID, genes_table$ensembl_gene_id)]

topGenes = DEGs$ID[1:5]

geneCounts <- plotCounts(dds, gene = topGenes[1], intgroup=c("cell", "condition"), returnData = TRUE)

# base
ggplot(geneCounts, aes(x = condition, y = count, color = cell, group = cell)) +
  scale_y_log10() + geom_point(size = 3) + geom_line() + 
  # Adding title with dynamic content
  ggtitle(str_c(DEGs$hgnc_symbol[1], "; log2FC:", 
                round(DEGs$log2FC[1], digits = 2), "; p-value: ", 
                round(DEGs$pvalue[1], digits = 2), "; adj.pval: ", 
                round(DEGs$padj[1], digits = 2))) +
  
  # Labels for X and Y axes
  labs(color = "Cell line") + 
  xlab("Conditions") +
  ylab("Expression count") +
  scale_x_discrete(labels = c("ctr" = "Control", "ki" = "KI")) +
  scale_color_manual(values = c("hct" = grp_col[2], "mpi" = grp_col[1]), labels = c("hct" = "HCT", "mpi" = "DLD")) + 
  theme_classic() + 
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border = element_blank(),
    axis.line = element_line(color = grp_col[2], size = 0.5),
    panel.background = element_blank(), # Background of the entire plot
    axis.title = element_text(size = 22, color = "black"),
    axis.ticks.length = unit(0.1, "cm"), 
    axis.text = element_text(size = 22, colour = "black"),
    plot.title = element_text(size = 18),
    legend.position = "top",   # Move legend to the top
    legend.title = element_text(size = 18, color = "black"),
    legend.key.size = unit(2, "line"),
    legend.text = element_text(size = 18, color = "black")
  )


# 1
geneCounts <- plotCounts(dds, gene = topGenes[1], intgroup=c("cell", "condition"), returnData = TRUE)

p <- ggplot(geneCounts, aes(x = condition, y = count, color = cell, group = cell)) +
  scale_y_log10() + geom_point(size = 3) + geom_line() + 
  # Adding title with dynamic content
  ggtitle(str_c(DEGs$hgnc_symbol[1], "; log2FC:", 
                round(DEGs$log2FC[1], digits = 2), "; p-value: ", 
                round(DEGs$pvalue[1], digits = 2), "; adj.pval: ", 
                round(DEGs$padj[1], digits = 2))) +
  
  # Labels for X and Y axes
  labs(color = "Cell line") + 
  xlab("Conditions") +
  ylab("Expression count") +
  scale_x_discrete(labels = c("ctr" = "Control", "ki" = "KI")) +
  scale_color_manual(values = c("hct" = grp_col[2], "mpi" = grp_col[1]), labels = c("hct" = "HCT", "mpi" = "DLD")) + 
  theme_classic() + 
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border = element_blank(),
    axis.line = element_line(color = grp_col[2], size = 0.5),
    panel.background = element_blank(), # Background of the entire plot
    axis.title = element_blank(),
    axis.ticks.length = unit(0.05, "cm"), 
    axis.text.x = element_blank(),
    axis.text.y = element_text(size = 18, color = "black"),
    plot.title = element_text(size = 18),
    legend.position = "none",   # Move legend to the top
    legend.title = element_text(size = 18, color = "black"),
    legend.key.size = unit(2, "line"),
    legend.text = element_text(size = 18, color = "black")
  )

ggsave(filename = "04_Figure_5_1.png", plot = p, width = 6, height = 6, dpi = 300)

# 2
geneCounts <- plotCounts(dds, gene = topGenes[2], intgroup=c("cell", "condition"), returnData = TRUE)

p <- ggplot(geneCounts, aes(x = condition, y = count, color = cell, group = cell)) +
  scale_y_log10() + geom_point(size = 3) + geom_line() + 
  # Adding title with dynamic content
  ggtitle(str_c(DEGs$hgnc_symbol[2], "; log2FC:", 
                round(DEGs$log2FC[2], digits = 2), "; p-value: ", 
                round(DEGs$pvalue[2], digits = 2), "; adj.pval: ", 
                round(DEGs$padj[2], digits = 2))) +
  
  # Labels for X and Y axes
  labs(color = "Cell line") + 
  xlab("Conditions") +
  ylab("Expression count") +
  scale_x_discrete(labels = c("ctr" = "Control", "ki" = "KI")) +
  scale_color_manual(values = c("hct" = grp_col[2], "mpi" = grp_col[1]), labels = c("hct" = "HCT", "mpi" = "DLD")) + 
  theme_classic() + 
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border = element_blank(),
    axis.line = element_line(color = grp_col[2], size = 0.5),
    panel.background = element_blank(), # Background of the entire plot
    axis.title = element_blank(),
    axis.ticks.length = unit(0.05, "cm"), 
    axis.text.x = element_blank(),
    axis.text.y = element_text(size = 18, color = "black"),
    plot.title = element_text(size = 18),
    legend.position = "none",   # Move legend to the top
    legend.title = element_text(size = 18, color = "black"),
    legend.key.size = unit(2, "line"),
    legend.text = element_text(size = 18, color = "black")
  )

ggsave(filename = "04_Figure_5_2.png", plot = p, width = 6, height = 6, dpi = 300)

# 3
geneCounts <- plotCounts(dds, gene = topGenes[3], intgroup=c("cell", "condition"), returnData = TRUE)

p <- ggplot(geneCounts, aes(x = condition, y = count, color = cell, group = cell)) +
  scale_y_log10() + geom_point(size = 3) + geom_line() + 
  # Adding title with dynamic content
  ggtitle(str_c(DEGs$hgnc_symbol[3], "; log2FC:", 
                round(DEGs$log2FC[3], digits = 2), "; p-value: ", 
                round(DEGs$pvalue[3], digits = 2), "; adj.pval: ", 
                round(DEGs$padj[3], digits = 2))) +
  
  # Labels for X and Y axes
  labs(color = "Cell line") + 
  xlab("Conditions") +
  ylab("Expression count") +
  scale_x_discrete(labels = c("ctr" = "Control", "ki" = "KI")) +
  scale_color_manual(values = c("hct" = grp_col[2], "mpi" = grp_col[1]), labels = c("hct" = "HCT", "mpi" = "DLD")) + 
  theme_classic() + 
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border = element_blank(),
    axis.line = element_line(color = grp_col[2], size = 0.5),
    panel.background = element_blank(), # Background of the entire plot
    axis.title = element_blank(),
    axis.ticks.length = unit(0.05, "cm"), 
    axis.text.x = element_blank(),
    axis.text.y = element_text(size = 18, color = "black"),
    title = element_text(size = 14),
    legend.position = "none",   # Move legend to the top
    legend.title = element_text(size = 18, color = "black"),
    legend.key.size = unit(2, "line"),
    legend.text = element_text(size = 18, color = "black")
  )

ggsave(filename = "04_Figure_5_3.png", plot = p, width = 6, height = 6, dpi = 300)

