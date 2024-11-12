# Specify the projects dir containing all the functions of scRNA-seq analysis pipeline
source("~/Documents/Projects/my_projects/m6A_RNA_methylation_in_colorectal_cancer/library_and_packages/library_handler.R", .GlobalEnv)

# Functions to perfrom analysis tasks
funcset <- list.files("~/Documents/Projects/my_projects/nanopore_RNA_DNA_analysis_functions", pattern = "*.R$", full.names = TRUE)
sapply(funcset, source, .GlobalEnv)

# Load the table with differential modification rate - filtered on all sites for p-value 0.05
dds = loadRData("~/Documents/Projects/Git_repositories/effect_of_METTL3_on_m6A_in_cancer/Analyses/Analysis_of_all_samples/01_whole_data_analysis/02_analysis_on_counts/analysis_output/05_deg_files/DESeq_object_hct_ctr_hct_ki_mpi_ctr_mpi_ki.RData")

vsd <- vst(dds, blind = FALSE)

head(assay(vsd), 3)

p <- plotPCA(vsd, intgroup = c("condition",  "cell"), ntop = 2000)

p <- p + 
  geom_point(size = 5) + 
  labs(color = "Samples") + 
  xlab("PC1 (97% variance)") +
  ylab("PC2 (3% variance)") +
  scale_color_manual(values = c("ctr:hct" = grp_col[2], "ctr:mpi" = grp_col[1], "ki:hct" = grp_col[3], "ki:mpi" = grp_col[4]), labels = c("Control HCT", "Control DLD", "KI HCT", "KI DLD")) + 
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

ggsave(filename = "01_Figure_5.png", plot = p, width = 10, height = 12, dpi = 300)


# Another way  
pcaData <- plotPCA(vsd, intgroup = c("condition",  "cell"), returnData = TRUE)
pcaData

percentVar <- round(100 * attr(pcaData, "percentVar"))

p <- ggplot(pcaData, aes(x = PC1, y = PC2, color = condition, shape = cell)) +
  geom_point(size = 5) +
  xlab(paste0("PC1 (", percentVar[1], "% variance)")) +
  ylab(paste0("PC2 (", percentVar[2], "% variance)")) + 
  coord_fixed() +
  scale_color_manual(values = c("ctr:hct" = grp_col[2], "ctr:mpi" = grp_col[1], "ki:hct" = grp_col[3], "ki:mpi" = grp_col[4]), labels = c("Control HCT", "Control DLD", "KI HCT", "KI DLD")) + 
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

