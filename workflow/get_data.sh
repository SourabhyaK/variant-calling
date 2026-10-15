#!/bin/bash
# Downloads the real E. coli REL606 reference genome and the real trimmed
# FASTQ subset used in the Data Carpentry "Wrangling Genomics" lesson.
# Run this once, from the project root, before `snakemake` or inside the container.
set -e

mkdir -p data/ref_genome data/trimmed_fastq

echo "Downloading E. coli REL606 reference genome (real NCBI assembly GCA_000017985.1)..."
curl -L -o data/ref_genome/ecoli_rel606.fasta.gz \
  https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/017/985/GCA_000017985.1_ASM1798v1/GCA_000017985.1_ASM1798v1_genomic.fna.gz
gunzip -f data/ref_genome/ecoli_rel606.fasta.gz

echo "Downloading real trimmed FASTQ subset (3 samples, Ara-3 strain, 3 timepoints)..."
curl -L -o sub.tar.gz https://ndownloader.figshare.com/files/14418248
tar xvf sub.tar.gz
mv sub/* data/trimmed_fastq/
rm -rf sub sub.tar.gz

echo "Done. Reference genome: $(grep '^>' data/ref_genome/ecoli_rel606.fasta)"
echo "FASTQ files ready in data/trimmed_fastq/"
