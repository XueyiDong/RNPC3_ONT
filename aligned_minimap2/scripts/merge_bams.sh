#!/bin/bash
#SBATCH --time=24:00:00
#SBATCH --cpus-per-task=8
#SBATCH --mem=200G
#SBATCH --job-name=merge_bams
#SBATCH --mail-user=dong.x@wehi.edu.au
#SBATCH --mail-type=BEGIN,END,FAIL

cd /stornext/General/data/user_managed/grpu_mritchie_1/XueyiDong/RNPC3/aligned_minimap2
module load samtools
samtools merge --threads 7 -o merged_NT.bam barcode01_fail.bam barcode01_pass.bam barcode02_fail.bam barcode02_pass.bam barcode03_fail.bam barcode03_pass.bam barcode04_fail.bam barcode04_pass.bam
samtools merge --threads 7 -o merged_RNPC3.bam barcode05_fail.bam barcode05_pass.bam barcode06_fail.bam barcode06_pass.bam barcode07_fail.bam barcode07_pass.bam barcode08_fail.bam barcode08_pass.bam