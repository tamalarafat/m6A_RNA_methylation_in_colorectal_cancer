# Specify the projects dir containing all the functions of scRNA-seq analysis pipeline
source("~/Documents/Projects/my_projects/m6A_RNA_methylation_in_colorectal_cancer/library_and_packages/library_handler.R", .GlobalEnv)

# Functions to perfrom analysis tasks
funcset <- list.files("~/Documents/Projects/my_projects/nanopore_RNA_DNA_analysis_functions", pattern = "*.R$", full.names = TRUE)
sapply(funcset, source, .GlobalEnv)

# Load the table with differential modification rate - filtered on all sites for p-value 0.05
dds = loadRData("~/Documents/Projects/Git_repositories/effect_of_METTL3_on_m6A_in_cancer/Analyses/Analysis_of_all_samples/01_whole_data_analysis/02_analysis_on_counts/analysis_output/05_deg_files/DESeq_object_hct_ctr_hct_ki_mpi_ctr_mpi_ki.RData")

# DEGs
res = loadRData("~/Documents/Projects/Git_repositories/effect_of_METTL3_on_m6A_in_cancer/Analyses/Analysis_of_all_samples/01_whole_data_analysis/02_analysis_on_counts/analysis_output/05_deg_files/DEG_DESeq_hct_ctr_hct_ki_mpi_ctr_mpi_ki.RData")

DEGs = data.frame(ID = rownames(res), log2FC = res$log2FoldChange, pVal = res$pvalue, adjpVal = res$padj)

# Load the table containing gene information
genes_table = read.table("~/Documents/Projects/Git_repositories/effect_of_METTL3_on_m6A_in_cancer/Analyses/Analysis_of_all_samples/01_whole_data_analysis/02_analysis_on_counts/analysis_output/03_gene_count_tables/transcript_gene_table.csv")

DEGs$hgnc_symbol = genes_table$hgnc_symbol[match(DEGs$ID, genes_table$ensembl_gene_id)]

DEGs$entrezgene_id = genes_table$entrezgene_id[match(DEGs$ID, genes_table$ensembl_gene_id)]

# METTL3
geneCounts <- plotCounts(dds, gene = "ENSG00000165819", intgroup=c("cell", "condition"), returnData = TRUE)

# METTL3 index
temp_ind = which(DEGs$ID %in% "ENSG00000165819")

p <- ggplot(geneCounts, aes(x = condition, y = count, color = cell, group = cell)) +
  scale_y_log10() + geom_point(size = 3) + geom_line() + 
  # Adding title with dynamic content
  ggtitle(str_c(DEGs$hgnc_symbol[temp_ind], "; log2FC:", 
                round(DEGs$log2FC[temp_ind], digits = 2), "; p-value: ", 
                round(DEGs$pVal[temp_ind], digits = 2), "; adj.pval: ", 
                round(DEGs$adjpVal[temp_ind], digits = 2))) +
  
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
    axis.title = element_text(size = 18, color = "black"),
    axis.ticks.length = unit(0.05, "cm"), 
    axis.text.x = element_text(size = 18, color = "black"),
    axis.text.y = element_text(size = 18, color = "black"),
    title = element_text(size = 14),
    legend.position = "top",   # Move legend to the top
    legend.title = element_text(size = 18, color = "black"),
    legend.key.size = unit(2, "line"),
    legend.text = element_text(size = 18, color = "black")
  )

ggsave(filename = "03_Figure_5_METTL3.png", plot = p, width = 6, height = 6, dpi = 300)







