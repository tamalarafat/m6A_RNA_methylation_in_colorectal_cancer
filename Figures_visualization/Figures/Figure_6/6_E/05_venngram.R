# Specify the projects dir containing all the functions of scRNA-seq analysis pipeline
source("~/Documents/Projects/my_projects/m6A_RNA_methylation_in_colorectal_cancer/library_and_packages/library_handler.R", .GlobalEnv)

# Functions to perfrom analysis tasks
funcset <- list.files("~/Documents/Projects/my_projects/nanopore_RNA_DNA_analysis_functions", pattern = "*.R$", full.names = TRUE)
sapply(funcset, source, .GlobalEnv)

# DEGs in control and DRACH
ctrl_deg = read.table("~/Documents/Projects/Git_repositories/effect_of_METTL3_on_m6A_in_cancer/Analyses/Analysis_of_all_samples/01_whole_data_analysis/02_analysis_on_counts/analysis_output/07_simulated_analysis_results/03_characteristic_degs/DEG_up_in_ctr_group.csv")

ctrl_deg = ctrl_deg[ctrl_deg$log2FC >= 0.25, ]

ctrl_deg = ctrl_deg[ctrl_deg$padj < 0.05, ]

# Load the table with differential modification rate - filtered on all sites for p-value 0.05
temp_mod = loadRData("../Figure_3/all_significant_sites.RData")

# control group
temp_control = temp_mod[temp_mod$diff_mod_rate > 0 & temp_mod$DRACH == "YES", ]

# how many k-mers were detected
temp_control = temp_control[!is.na(temp_control$mers), ]

x = list(drach_genes = temp_control$ensembl_id,
         ctrl_genes = ctrl_deg$ID)

p <- ggvenn::ggvenn(
  x, 
  fill_color = c(grp_col[1], grp_col[2], grp_col[6]),
  fill_alpha = 0.7,
  stroke_size = 0.5, 
  set_name_size = 0,
  text_size = 10
)

ggsave(plot = p, filename = "05_Figure_6.png", width = 10, height = 10, dpi = 300, bg = "white")
