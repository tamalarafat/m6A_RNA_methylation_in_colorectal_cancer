#!/bin/bash

# Directory containing the aligned basecalled reads
data_dir="/ds/storage/common/ytamal_work/Project_14052024_ONT_RNAseq/PAS98818_FLO_PRO004RA_SQK_RNA004_HCT_ctr/analysis_of_the_sample/Basecalled_aligned_data/basecalled_reads.bam"

# Output directory of the indexed bam file
output_dir="/ds/storage/common/ytamal_work/Project_14052024_ONT_RNAseq/PAS98818_FLO_PRO004RA_SQK_RNA004_HCT_ctr/analysis_of_the_sample/analysis_output"

# Create output directory if it doesn't exist
mkdir -p "$output_dir"

# Sort (by coordinate) the bam file with samtools
sorted_bam="$output_dir/basecalled_reads_sorted.bam"

samtools sort -o "$sorted_bam" "$data_dir"

# Index the sorted bam file with samtools
samtools index "$sorted_bam"
