# Specify the projects dir containing all the functions of scRNA-seq analysis pipeline
source("~/Documents/Projects/my_projects/m6A_RNA_methylation_in_colorectal_cancer/library_and_packages/library_handler.R", .GlobalEnv)

# Functions to perfrom analysis tasks
funcset <- list.files("~/Documents/Projects/my_projects/nanopore_RNA_DNA_analysis_functions", pattern = "*.R$", full.names = TRUE)
sapply(funcset, source, .GlobalEnv)

# Load the table containing gene information
genes_table = read.table("~/Documents/Projects/Git_repositories/effect_of_METTL3_on_m6A_in_cancer/Analyses/Analysis_of_all_samples/01_whole_data_analysis/02_analysis_on_counts/analysis_output/03_gene_count_tables/transcript_gene_table.csv")

# DEGs in control and DRACH
ctrl_deg = read.table("~/Documents/Projects/Git_repositories/effect_of_METTL3_on_m6A_in_cancer/Analyses/Analysis_of_all_samples/01_whole_data_analysis/02_analysis_on_counts/analysis_output/07_simulated_analysis_results/03_characteristic_degs/DEG_up_in_ctr_group.csv")

ctrl_deg$hgnc_symbol = genes_table$hgnc_symbol[match(ctrl_deg$ID, genes_table$ensembl_gene_id)]

ctrl_deg$entrezgene_id = genes_table$entrezgene_id[match(ctrl_deg$ID, genes_table$ensembl_gene_id)]

ctrl_deg = ctrl_deg[ctrl_deg$log2FC >= 0.25, ]

ctrl_deg = ctrl_deg[ctrl_deg$padj < 0.05, ]

# Load the table with differential modification rate - filtered on all sites for p-value 0.05
temp_mod = loadRData("../Figure_3/all_significant_sites.RData")

# control group
temp_control = temp_mod[temp_mod$diff_mod_rate > 0 & temp_mod$DRACH == "YES", ]

# how many k-mers were detected
temp_control = temp_control[!is.na(temp_control$mers), ]

temp_DRACH_genes = intersect(ctrl_deg$ID, temp_control$ensembl_id)

# Subset the control data
temp_ent_id = ctrl_deg$entrezgene_id[ctrl_deg$ID %in% temp_DRACH_genes]

# Directory to store the results
res_dir = "~/Documents/Projects/Git_repositories/effect_of_METTL3_on_m6A_in_cancer/Analyses/Analysis_of_all_samples/01_whole_data_analysis/analysis_output/"

# Load all entrez ids of the hg
hg_entrz = read.table(str_c(res_dir, "08_genes_and_their_function/hg_complete_set_entrez_id.txt"))

hg_entrz = as.character(str_sort(hg_entrz[, 1], numeric = TRUE))

hg_entrz = hg_entrz[!is.na(hg_entrz)]

# performing the over representation analysis (ORA) for Gene Ontology class - Biological processes
GeneSet_ORA_BP <- enrichGO(gene = temp_ent_id,
                           universe = hg_entrz,
                           OrgDb = org.Hs.eg.db,
                           keyType = "ENTREZID", 
                           ont = "BP",
                           pAdjustMethod = "BH",
                           pvalueCutoff = 0.05,
                           readable = TRUE,
                           pool = FALSE)

p <- dotplot(GeneSet_ORA_BP, showCategory = 10) + 
  labs(colour = "p.adj") + 
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        axis.line = element_line(color = grp_col[2]),
        axis.title.x = element_text(size = 22, face = "bold", color = "black"), 
        axis.ticks.length = unit(.30, "cm"), 
        axis.text.x = element_text(size = 22, color = "black", face = "bold"),
        axis.text.y = element_text(size = 22, color = "black", face = "bold"),
        legend.key = element_rect(size = 22),
        legend.text = element_text(size = 22),
        legend.title = element_text(size = 22, face = "bold"),
        legend.spacing = unit(2.0, 'cm'),
        legend.key.size = unit(3,"line"),
        legend.position = "none")

ggsave(filename = "06_Figure_6.png", plot = p, width = 14, height = 14, dpi = 300)

GeneSet_ORA_BP <- pairwise_termsim(GeneSet_ORA_BP, method = "JC")

p2 <- emapplot(GeneSet_ORA_BP, color = "qvalue", cex_label_category = 2.5, cex_category = 2)  + theme_bw() + 
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        axis.title = element_blank(), 
        text = element_text(size = 24, face = "bold"),
        axis.line = element_blank(),
        axis.ticks.length = unit(0, "cm"), 
        axis.text = element_blank(),
        legend.key = element_rect(size = 22),
        legend.spacing = unit(2.0, 'cm'),
        legend.key.size = unit(3,"line"),
        legend.position = "none")
ggsave(filename = "06_Figure_6_network.png", plot = p2, width = 22, height = 22, dpi = 300)
