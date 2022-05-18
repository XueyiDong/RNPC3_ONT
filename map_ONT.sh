#!/bin/bash
#SBATCH --time=24:00:00
#SBATCH --cpus-per-task=8
#SBATCH --mem=20G
#SBATCH --job-name=map_RNPC3
#SBATCH --mail-user=dong.x@wehi.edu.au
#SBATCH --mail-type=BEGIN,END,FAIL

cd /stornext/General/data/user_managed/grpu_mritchie_1/XueyiDong/RNPC3/aligned_minimap2
# unixhome=/wehisan/home/allstaff/d/dong.x
fq=/stornext/Projects/promethion/promethion_access/lab_ritchie/A549_RNPC3/long_term/
genome=/stornext/General/data/user_managed/grpu_mritchie_1/XueyiDong/annotation/Human/CHM13v2/GCF_009914755.1_T2T-CHM13v2.0_genomic.fna
bed=/stornext/General/data/user_managed/grpu_mritchie_1/XueyiDong/annotation/Human/CHM13v2/anno.bed
# using minimap2 2.24
export PATH=$PATH:/wehisan/home/allstaff/d/dong.x/Programs/minimap2

module load samtools/1.7

minimap2 -ax splice -k14 --junc-bed  $bed $genome $fq/$1.fq.gz | samtools view -b | samtools sort > $1.bam
samtools index $1.bam $1.bai

# for sample in barcode{01..08} 
# do minimap2 -ax splice -k14 --junc-bed  $bed $genome $fq/$sample\_pass.fq.gz | samtools view -b | samtools sort > $sample.pass.bam
# samtools index $sample.pass.bam $sample.pass.bai
# minimap2 -ax splice -k14 --junc-bed  $bed $genome $fq/$sample\_fail.fq.gz | samtools view -b | samtools sort > $sample.fail.bam
# samtools index $sample.fail.bam $sample.fail.bai
# done