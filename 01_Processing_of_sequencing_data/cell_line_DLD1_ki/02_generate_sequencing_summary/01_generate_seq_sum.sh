#!/bin/bash

# Directory containing the aligned basecalled reads
data_dir="/ds/storage/common/ytamal_work/Project_14052024_ONT_RNAseq/PAS97975_FLO_PRO004RA_SQK_RNA004_MPI_Ki/analysis_of_the_sample/Basecalled_aligned_data/basecalled_reads.bam"

# Output directory of the indexed bam file
output_dir="/ds/storage/common/ytamal_work/Project_14052024_ONT_RNAseq/PAS97975_FLO_PRO004RA_SQK_RNA004_MPI_Ki/analysis_of_the_sample/analysis_output"

# Create output directory if it doesn't exist
mkdir -p "$output_dir"

# Convert all pod5 files to fast5 files
time $ONTrna_dorado summary "$data_dir" > "$output_dir/DLD_ki_summary.tsv"

