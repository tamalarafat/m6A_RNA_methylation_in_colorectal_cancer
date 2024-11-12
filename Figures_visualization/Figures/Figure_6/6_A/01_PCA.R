# Specify the projects dir containing all the functions of scRNA-seq analysis pipeline
source("~/Documents/Projects/my_projects/m6A_RNA_methylation_in_colorectal_cancer/library_and_packages/library_handler.R", .GlobalEnv)

# Functions to perfrom analysis tasks
funcset <- list.files("~/Documents/Projects/my_projects/nanopore_RNA_DNA_analysis_functions", pattern = "*.R$", full.names = TRUE)
sapply(funcset, source, .GlobalEnv)

# Load the table with differential modification rate - filtered on all sites for p-value 0.05
dds = loadRData("~/Documents/Projects/Git_repositories/effect_of_METTL3_on_m6A_in_cancer/Analyses/Analysis_of_all_samples/01_whole_data_analysis/02_analysis_on_counts/analysis_output/07_simulated_analysis_results/02_deg_files/DESeq_object_rep_1_mpi_ctr_rep_1_mpi_ki_rep_2_mpi_ctr_rep_2_mpi_ki_rep_3_mpi_ctr_rep_3_mpi_ki.RData")

vsd <- vst(dds, blind = FALSE)

head(assay(vsd), 3)

p <- plotPCA(vsd, intgroup = c("condition",  "cell"), ntop = 2000)

p <- p + 
  geom_point(size = 5) + 
  labs(color = "DLD cell line") + 
  xlab("PC1 (91% variance)") +
  ylab("PC2 (5% variance)") +
  scale_color_manual(values = c("ctr:rep" = grp_col[2], "ki:rep" = grp_col[1]), labels = c("Control", "Knock-down")) + 
  theme_classic() + 
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border = element_blank(),
    axis.line = element_line(color = grp_col[2], size = 0.5),
    panel.background = element_blank(), # Background of the entire plot
    axis.title = element_text(size = 22, color = "black"),
    axis.ticks.length = unit(0, "cm"), 
    axis.text = element_text(size = 22, colour = "black"),
    legend.position = "top",
    legend.title = element_text(size = 16, colour = "black", face = "bold"),
    legend.key.size = unit(1, "line"), 
    legend.text = element_text(size = 22)) + 
  guides(color = guide_legend(override.aes = list(size = 4)))

p

ggsave(filename = "01_Figure_6.png", plot = p, width = 10, height = 12, dpi = 300)

