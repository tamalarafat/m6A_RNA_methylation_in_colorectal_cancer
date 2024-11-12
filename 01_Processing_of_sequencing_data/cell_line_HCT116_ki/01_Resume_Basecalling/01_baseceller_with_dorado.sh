#!/bin/bash

# Directory containing the pod5 files
data_dir="/ds/storage/common/ytamal_work/Project_14052024_ONT_RNAseq/PAS98400_FLO_PRO004RA_SQK_RNA004_HCT_Ki/analysis_of_the_sample/pod5_pass/"

# Directory - the output will be saved
output_dir="/ds/storage/common/ytamal_work/Project_14052024_ONT_RNAseq/PAS98400_FLO_PRO004RA_SQK_RNA004_HCT_Ki/analysis_of_the_sample/Basecalled_aligned_data/"

# Directory with the file to resume the basecalling from
resume_file="/ds/storage/common/ytamal_work/Project_14052024_ONT_RNAseq/PAS98400_FLO_PRO004RA_SQK_RNA004_HCT_Ki/analysis_of_the_sample/Basecalled_aligned_data/basecalled_reads.bam"

REF_DIR="/ds/storage/common/ytamal_work/Project_14052024_ONT_RNAseq/References/Transcripts/Homo_sapiens.GRCh38.cdna.all.fa"

# Create output directory if it doesn't exist
mkdir -p "$output_dir"

# Resume the basecalling from previous run
time $ONTrna_dorado basecaller rna004_130bps_hac@v5.0.0 "$data_dir" --modified-bases m6A --reference "$REF_DIR" --resume-from "$resume_file" > "$output_dir/basecalled_reads_resumed.bam"
