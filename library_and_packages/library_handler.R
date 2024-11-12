# Set CRAN mirror for package installation
options(repos = c(CRAN = "https://cloud.r-project.org"))

# Check for BiocManager and install if necessary
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

# Function to install CRAN packages if not already installed
install_if_missing <- function(pkg) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    install.packages(pkg, dependencies = FALSE)  # Avoid installing dependencies unless necessary
  }
}

# Function to install Bioconductor packages without updating existing packages
install_bioc_if_missing <- function(pkg) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    BiocManager::install(pkg, ask = FALSE, update = FALSE)
  }
}

# Install required CRAN packages
cran_packages <- c("ggplot2", "ggrepel", "ggseqlogo", "reshape2", "magrittr", "stringr", "readr", "yaml", "dplyr", "seqgendiff")
lapply(cran_packages, install_if_missing)

# Install required Bioconductor packages
bioc_packages <- c("enrichplot", "org.Hs.eg.db", "clusterProfiler", "tximport", "DESeq2", "Biostrings", "biomaRt", "seqinr")
lapply(bioc_packages, install_bioc_if_missing)

# Load necessary libraries
suppressMessages(library(ggplot2))
suppressMessages(library(ggrepel))
suppressMessages(library(ggseqlogo))
suppressMessages(library(reshape2))
suppressMessages(library(magrittr))
suppressMessages(library(stringr))
suppressMessages(library(readr))
suppressMessages(library(yaml))
suppressMessages(library(dplyr))
suppressMessages(library(Biostrings))
suppressMessages(library(biomaRt))
suppressMessages(library(seqinr))
suppressMessages(library(enrichplot))
suppressMessages(library(org.Hs.eg.db))
suppressMessages(library(clusterProfiler))
suppressMessages(library(tximport))
suppressMessages(library(DESeq2))
suppressMessages(library(seqgendiff))

# Define the colors
grp_col <- c("#FD6467", "#00A08A", "#F98400", "#046C9A", "#075149FF", "#FFA319FF", "#00BF81", "#767676FF", "#FD8CC1FF")
