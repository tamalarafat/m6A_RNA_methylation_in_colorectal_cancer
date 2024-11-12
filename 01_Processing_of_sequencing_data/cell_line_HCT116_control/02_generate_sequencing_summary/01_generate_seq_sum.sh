#!/bin/bash

data_dir="/ds/storage/common/ytamal_work/Project_14052024_ONT_RNAseq/PAS98818_FLO_PRO004RA_SQK_RNA004_HCT_ctr/analysis_of_the_sample/Basecalled_aligned_data/basecalled_reads.bam"

output_dir="/ds/storage/common/ytamal_work/Project_14052024_ONT_RNAseq/PAS98818_FLO_PRO004RA_SQK_RNA004_HCT_ctr/analysis_of_the_sample/analysis_output"

# Create output directory if it doesn't exist
mkdir -p "$output_dir"

# Convert all pod5 files to fast5 files
time $ONTrna_dorado summary "$data_dir" > "$output_dir/HCT_ctr_summary.tsv"

