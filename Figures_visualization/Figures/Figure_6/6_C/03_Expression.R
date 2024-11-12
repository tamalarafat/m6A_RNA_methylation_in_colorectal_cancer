# Specify the projects dir containing all the functions of scRNA-seq analysis pipeline
source("~/Documents/Projects/my_projects/m6A_RNA_methylation_in_colorectal_cancer/library_and_packages/library_handler.R", .GlobalEnv)

# Functions to perfrom analysis tasks
funcset <- list.files("~/Documents/Projects/my_projects/nanopore_RNA_DNA_analysis_functions", pattern = "*.R$", full.names = TRUE)
sapply(funcset, source, .GlobalEnv)

# Load the table with differential modification rate - filtered on all sites for p-value 0.05
dds = loadRData("~/Documents/Projects/Git_repositories/effect_of_METTL3_on_m6A_in_cancer/Analyses/Analysis_of_all_samples/01_whole_data_analysis/02_analysis_on_counts/analysis_output/07_simulated_analysis_results/02_deg_files/DESeq_object_rep_1_mpi_ctr_rep_1_mpi_ki_rep_2_mpi_ctr_rep_2_mpi_ki_rep_3_mpi_ctr_rep_3_mpi_ki.RData")

# Load the table containing gene information
genes_table = read.table("~/Documents/Projects/Git_repositories/effect_of_METTL3_on_m6A_in_cancer/Analyses/Analysis_of_all_samples/01_whole_data_analysis/02_analysis_on_counts/analysis_output/03_gene_count_tables/transcript_gene_table.csv")

# KI
KI = read.table("~/Documents/Projects/Git_repositories/effect_of_METTL3_on_m6A_in_cancer/Analyses/Analysis_of_all_samples/01_whole_data_analysis/02_analysis_on_counts/analysis_output/07_simulated_analysis_results/03_characteristic_degs/DEG_up_in_ki_group.csv")

KI = KI[abs(KI$log2FC) >= 0.25, ]

KI = KI[abs(KI$padj) < 0.05, ]

# DEGs in control and DRACH
DEGs = read.table("~/Documents/Projects/Git_repositories/effect_of_METTL3_on_m6A_in_cancer/Analyses/Analysis_of_all_samples/01_whole_data_analysis/02_analysis_on_counts/analysis_output/07_simulated_analysis_results/03_characteristic_degs/DEG_up_in_ctr_group.csv")

DEGs$hgnc_symbol = genes_table$hgnc_symbol[match(DEGs$ID, genes_table$ensembl_gene_id)]

DEGs$entrezgene_id = genes_table$entrezgene_id[match(DEGs$ID, genes_table$ensembl_gene_id)]

DEGs = DEGs[DEGs$log2FC >= 0.25, ]

DEGs = DEGs[DEGs$padj < 0.05, ]

# Load the table with differential modification rate - filtered on all sites for p-value 0.05
temp_mod = loadRData("../Figure_3/all_significant_sites.RData")

# control group
temp_control = temp_mod[temp_mod$diff_mod_rate > 0 & temp_mod$DRACH == "YES", ]

# how many k-mers were detected
temp_control = temp_control[!is.na(temp_control$mers), ]

temp_DRACH_genes = intersect(DEGs$ID, temp_control$ensembl_id)

topGenes = temp_DRACH_genes

# METTL3
geneCounts <- plotCounts(dds, gene = "ENSG00000165819", returnData = TRUE)

# METTL3 index
temp_ident = which(DEGs$ID %in% "ENSG00000165819")

p <- ggplot(geneCounts, aes(x = condition, y = count, color = condition, group = condition)) +
  scale_y_log10() + geom_point(size = 3) + geom_line() + 
  # Adding title with dynamic content
  ggtitle(str_c(DEGs$hgnc_symbol[temp_ident], "; log2FC:", 
                round(DEGs$log2FC[temp_ident], digits = 2), "; p-value: ", 
                round(DEGs$pvalue[temp_ident], digits = 2), "; adj.pval: ", 
                round(DEGs$padj[temp_ident], digits = 2))) +
  
  # Labels for X and Y axes
  labs(color = "DLD cell line") + 
  xlab("Conditions") +
  ylab("Expression count") +
  scale_x_discrete(labels = c("ctr" = "Control", "ki" = "KI")) +
  scale_color_manual(values = c("ctr" = grp_col[2], "ki" = grp_col[1]), labels = c("ctr" = "Control", "ki" = "Knock-down")) + 
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

ggsave(filename = "03_Figure_6_METTL3.png", plot = p, width = 6, height = 6, dpi = 300)

# 1
temp_ident = 1

geneCounts <- plotCounts(dds, gene = topGenes[temp_ident], intgroup=c("cell", "condition"), returnData = TRUE)

p <- ggplot(geneCounts, aes(x = condition, y = count, color = condition, group = condition)) +
  scale_y_log10() + geom_point(size = 3) + geom_line() + 
  # Adding title with dynamic content
  ggtitle(str_c(DEGs$hgnc_symbol[temp_ident], "; log2FC:", 
                round(DEGs$log2FC[temp_ident], digits = 2), "; p-value: ", 
                round(DEGs$pvalue[temp_ident], digits = 2), "; adj.pval: ", 
                round(DEGs$padj[temp_ident], digits = 2))) +
  
  # Labels for X and Y axes
  labs(color = "DLD cell line") + 
  xlab("Conditions") +
  ylab("Expression count") +
  scale_x_discrete(labels = c("ctr" = "Control", "ki" = "KI")) +
  scale_color_manual(values = c("ctr" = grp_col[2], "ki" = grp_col[1]), labels = c("ctr" = "Control", "ki" = "Knock-down")) + 
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

ggsave(filename = "03_Figure_6_1.png", plot = p, width = 6, height = 6, dpi = 300)

# 2
temp_ident = 2

geneCounts <- plotCounts(dds, gene = topGenes[temp_ident], intgroup=c("cell", "condition"), returnData = TRUE)

p <- ggplot(geneCounts, aes(x = condition, y = count, color = condition, group = condition)) +
  scale_y_log10() + geom_point(size = 3) + geom_line() + 
  # Adding title with dynamic content
  ggtitle(str_c(DEGs$hgnc_symbol[temp_ident], "; log2FC:", 
                round(DEGs$log2FC[temp_ident], digits = 2), "; p-value: ", 
                round(DEGs$pvalue[temp_ident], digits = 2), "; adj.pval: ", 
                round(DEGs$padj[temp_ident], digits = 2))) +
  
  # Labels for X and Y axes
  labs(color = "DLD cell line") + 
  xlab("Conditions") +
  ylab("Expression count") +
  scale_x_discrete(labels = c("ctr" = "Control", "ki" = "KI")) +
  scale_color_manual(values = c("ctr" = grp_col[2], "ki" = grp_col[1]), labels = c("ctr" = "Control", "ki" = "Knock-down")) + 
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

ggsave(filename = "03_Figure_6_2.png", plot = p, width = 6, height = 6, dpi = 300)

# 3
temp_ident = 3

geneCounts <- plotCounts(dds, gene = topGenes[temp_ident], intgroup=c("cell", "condition"), returnData = TRUE)

p <- ggplot(geneCounts, aes(x = condition, y = count, color = condition, group = condition)) +
  scale_y_log10() + geom_point(size = 3) + geom_line() + 
  # Adding title with dynamic content
  ggtitle(str_c(DEGs$hgnc_symbol[temp_ident], "; log2FC:", 
                round(DEGs$log2FC[temp_ident], digits = 2), "; p-value: ", 
                round(DEGs$pvalue[temp_ident], digits = 2), "; adj.pval: ", 
                round(DEGs$padj[temp_ident], digits = 2))) +
  
  # Labels for X and Y axes
  labs(color = "DLD cell line") + 
  xlab("Conditions") +
  ylab("Expression count") +
  scale_x_discrete(labels = c("ctr" = "Control", "ki" = "KI")) +
  scale_color_manual(values = c("ctr" = grp_col[2], "ki" = grp_col[1]), labels = c("ctr" = "Control", "ki" = "Knock-down")) + 
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

ggsave(filename = "03_Figure_6_3.png", plot = p, width = 6, height = 6, dpi = 300)




