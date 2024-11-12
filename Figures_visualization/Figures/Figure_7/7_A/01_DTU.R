# Specify the projects dir containing all the functions of scRNA-seq analysis pipeline
source("~/Documents/Projects/my_projects/m6A_RNA_methylation_in_colorectal_cancer/library_and_packages/library_handler.R", .GlobalEnv)

# Functions to perfrom analysis tasks
funcset <- list.files("~/Documents/Projects/my_projects/nanopore_RNA_DNA_analysis_functions", pattern = "*.R$", full.names = TRUE)
sapply(funcset, source, .GlobalEnv)

# Load the DRIMSeq object
dds <- loadRData("~/Documents/Projects/Git_repositories/effect_of_METTL3_on_m6A_in_cancer/Analyses/Analysis_of_all_samples/01_whole_data_analysis/02_analysis_on_counts/03_DET_analysis/DTU.RData")

res <- DRIMSeq::results(dds)

head(res)

res.txp <- DRIMSeq::results(dds, level="feature")
head(res.txp)

no.na <- function(x) ifelse(is.na(x), 1, x)
res$pvalue <- no.na(res$pvalue)
res.txp$pvalue <- no.na(res.txp$pvalue)

idx <- which(res$pvalue < 0.05)

temp_df = res[idx,]

temp_df = temp_df[order(temp_df$lr, decreasing = TRUE), ]

# Load the table containing gene information
genes_table = read.table("~/Documents/Projects/Git_repositories/effect_of_METTL3_on_m6A_in_cancer/Analyses/Analysis_of_all_samples/01_whole_data_analysis/02_analysis_on_counts/analysis_output/03_gene_count_tables/transcript_gene_table.csv")

# Temp indexing - METTL3
temp_gene = "ENSG00000165819.12"

temp_id = which(genes_table$ensembl_gene_id %in% str_replace(temp_gene, pattern = "\\..*", replacement = ""))[1]

p <- plotProportions(dds, temp_gene, "condition") + 
  # Adding title with dynamic content
  ggtitle(genes_table$hgnc_symbol[temp_id]) +
  # Labels for X and Y axes
  xlab("Transcript id") +
  ylab("Proportion usage") 

ggsave(filename = "01_DTU_1.png", plot = p, width = 6, height = 6, dpi = 300)


# Temp indexing - 1
temp_gene = temp_df$gene_id[1]

temp_id = which(genes_table$ensembl_gene_id %in% str_replace(temp_gene, pattern = "\\..*", replacement = ""))[1]

p <- plotProportions(dds, temp_gene, "condition") + 
  # Adding title with dynamic content
  ggtitle(genes_table$hgnc_symbol[temp_id]) +
  # Labels for X and Y axes
  xlab("Transcript id") +
  ylab("Proportion usage") 

ggsave(filename = "01_DTU_2.png", plot = p, width = 6, height = 6, dpi = 300)


# Temp indexing - 2
temp_gene = temp_df$gene_id[2]

temp_id = which(genes_table$ensembl_gene_id %in% str_replace(temp_gene, pattern = "\\..*", replacement = ""))[1]

p <- plotProportions(dds, temp_gene, "condition") + 
  # Adding title with dynamic content
  ggtitle(genes_table$hgnc_symbol[temp_id]) +
  # Labels for X and Y axes
  xlab("Transcript id") +
  ylab("Proportion usage") 

ggsave(filename = "01_DTU_3.png", plot = p, width = 6, height = 6, dpi = 300)

# Temp indexing - 4
temp_gene = temp_df$gene_id[3]

temp_id = which(genes_table$ensembl_gene_id %in% str_replace(temp_gene, pattern = "\\..*", replacement = ""))[1]

p <- plotProportions(dds, temp_gene, "condition") + 
  # Adding title with dynamic content
  ggtitle(genes_table$hgnc_symbol[temp_id]) +
  # Labels for X and Y axes
  xlab("Transcript id") +
  ylab("Proportion usage") 

ggsave(filename = "01_DTU_4.png", plot = p, width = 6, height = 6, dpi = 300)
