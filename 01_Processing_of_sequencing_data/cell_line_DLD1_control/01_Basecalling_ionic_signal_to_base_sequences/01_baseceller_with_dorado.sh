#!/bin/bash

data_dir="/ds/storage/common/ytamal_work/Project_14052024_ONT_RNAseq/PAS98171_FLO_PRO004RA_SQK_RNA004_MPI_ctr/analysis_of_the_sample/pod5_pass/"

output_dir="/ds/storage/common/ytamal_work/Project_14052024_ONT_RNAseq/PAS98171_FLO_PRO004RA_SQK_RNA004_MPI_ctr/analysis_of_the_sample/Basecalled_aligned_data/"

REF_DIR="/ds/storage/common/ytamal_work/Project_14052024_ONT_RNAseq/References/Transcripts/Homo_sapiens.GRCh38.cdna.all.fa"

# Create output directory if it doesn't exist
mkdir -p "$output_dir"

time $ONTrna_dorado basecaller rna004_130bps_hac@v5.0.0 "$data_dir" --modified-bases m6A --reference "$REF_DIR" > "$output_dir/basecalled_reads.bam"


