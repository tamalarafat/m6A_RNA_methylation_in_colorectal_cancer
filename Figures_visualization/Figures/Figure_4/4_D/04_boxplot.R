# Specify the projects dir containing all the functions of scRNA-seq analysis pipeline
source("~/Documents/Projects/my_projects/m6A_RNA_methylation_in_colorectal_cancer/library_and_packages/library_handler.R", .GlobalEnv)

# Functions to perfrom analysis tasks
funcset <- list.files("~/Documents/Projects/my_projects/nanopore_RNA_DNA_analysis_functions", pattern = "*.R$", full.names = TRUE)
sapply(funcset, source, .GlobalEnv)

# Load the table with differential modification rate - filtered on all sites for p-value 0.05
temp_df = loadRData("../Figure_3/all_significant_sites.RData")

# control group
df = temp_df[temp_df$DRACH == "YES", ]

# Containing only the columns with modification rate data
df = df[ , grep(pattern = "avg_mod_rate", colnames(df))]

# Change the column names to conditions
colnames(df) = str_replace(string = colnames(df), pattern = "avg_mod_rate_", replacement = "")

# Transform the data into long format for plotting
df_long = reshape2::melt(df, variable.name = "condition", value.name = "mod_rate")

# gene transcript table:
temp_gt = loadRData("~/Documents/Projects/Git_repositories/effect_of_METTL3_on_m6A_in_cancer/Analyses/Analysis_of_all_samples/01_whole_data_analysis/02_analysis_on_counts/03_DET_analysis/gene_transcript_count.RData")

temp_trans = unique(temp_df$ref_seq_id)

temp_gene = temp_gt[match(temp_trans, temp_gt$TXNAME), ]

# Variance test (checks for equality)
res.ftest <- var.test(mod_rate ~ condition, data = df_long)
res.ftest

# Compute t-test
res <- t.test(df$ctrl, df$ki, var.equal = FALSE, alternative = "two.sided")
res$p.value

# The p-value of the test is 2.2e-16, which is less than the significance level alpha = 0.05. 
# We can conclude that control group’s average modification rate is significantly different from KI condition's average modification rate with a p-value = 2.2e-16

# Boxplot of differential modification rate for all the shared sites across the samples of two conditions
p <- ggplot(df_long, aes(x = condition, y = mod_rate)) +
  geom_boxplot(coef = 6) +
  #geom_dotplot(binaxis='y', stackdir='center', dotsize = 0.4, binwidth = 1/40) + 
  #geom_jitter(shape = 16, position = position_jitter(0.1)) +
  theme_classic() +
  labs(
    title = str_c(nrow(df), "DRACH sites;", length(unique(temp_df$ref_seq_id)), "transcripts;", length(unique(temp_gene$GENEID)), "genes; \nsignificant difference between conditions: p-val (t-test): 7.2e-219", sep = " "),
    x = "Conditions",
    y = "Modification rate"
  ) +
  scale_x_discrete(labels = c("ctrl" = "Control", "ki" = "Knock-down")) + 
  theme(
    line = element_blank(),
    panel.border = element_blank(),
    title = element_text(size = 10),
    axis.title = element_text(size = 18),
    axis.text = element_text(size = 18))

ggsave(filename = "04_Figure_4.png", plot = p, width = 6, height = 6, dpi = 300)

# Variance test (checks for equality)
res.ftest <- var.test(mod_rate ~ condition, data = df_long)
res.ftest

# Compute t-test
res <- t.test(df$ctrl, df$ki, var.equal = FALSE, alternative = "two.sided")
res$p.value

# The p-value of the test is 2.2e-16, which is less than the significance level alpha = 0.05. 
# We can conclude that control group’s average modification rate is significantly different from KI condition's average modification rate with a p-value = 2.2e-16
