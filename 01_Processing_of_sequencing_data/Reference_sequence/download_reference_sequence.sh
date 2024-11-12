#!/bin/bash

# Download and extract the transcriptome reference sequence
wget -P http://ftp.ensembl.org/pub/release-104/fasta/homo_sapiens/cdna/Homo_sapiens.GRCh38.cdna.all.fa.gz
gunzip -c Homo_sapiens.GRCh38.cdna.all.fa.gz > Homo_sapiens.GRCh38.cdna.all.fa


