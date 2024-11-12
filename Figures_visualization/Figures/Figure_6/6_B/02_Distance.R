# Specify the projects dir containing all the functions of scRNA-seq analysis pipeline
source("~/Documents/Projects/my_projects/m6A_RNA_methylation_in_colorectal_cancer/library_and_packages/library_handler.R", .GlobalEnv)

# Functions to perfrom analysis tasks
funcset <- list.files("~/Documents/Projects/my_projects/nanopore_RNA_DNA_analysis_functions", pattern = "*.R$", full.names = TRUE)
sapply(funcset, source, .GlobalEnv)

# Load the table with differential modification rate - filtered on all sites for p-value 0.05
dds = loadRData("~/Documents/Projects/Git_repositories/effect_of_METTL3_on_m6A_in_cancer/Analyses/Analysis_of_all_samples/01_whole_data_analysis/02_analysis_on_counts/analysis_output/07_simulated_analysis_results/02_deg_files/DESeq_object_rep_1_mpi_ctr_rep_1_mpi_ki_rep_2_mpi_ctr_rep_2_mpi_ki_rep_3_mpi_ctr_rep_3_mpi_ki.RData")

vsd <- vst(dds, blind = FALSE)

head(assay(vsd), 3)

sampleDists <- dist(t(assay(vsd)))
sampleDists

sampleDistMatrix <- as.matrix( sampleDists )
rownames(sampleDistMatrix) <- paste(vsd$names)
colnames(sampleDistMatrix) <- NULL
colors <- colorRampPalette( rev(brewer.pal(9, "Blues")) )(255)

p <- pheatmap(sampleDistMatrix,
              clustering_distance_rows = sampleDists,
              clustering_distance_cols = sampleDists,
              col = colors, 
              legend_breaks = c(0, 20, 40, 60, 80, 100), 
              cluster_cols = FALSE, 
              cluster_rows = FALSE,
              show_rownames = TRUE, 
              show_colnames = TRUE
)

ggsave(filename = "02_Figure_6.png", plot = p, width = 6, height = 6, dpi = 300)

