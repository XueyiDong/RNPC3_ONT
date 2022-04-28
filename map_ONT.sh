#!/bin/bash
#SBATCH --time=48:00:00
#SBATCH --cpus-per-task=8
#SBATCH --mem=150G
#SBATCH --job-name=map_RNPC3
#SBATCH --output=map_RNPC3_fail.out
#SBATCH --mail-user=dong.x@wehi.edu.au
#SBATCH --mail-type=BEGIN,END,FAIL

cd /stornext/General/data/user_managed/grpu_mritchie_1/XueyiDong/long_read_benchmark/ONT/bam
# unixhome=/wehisan/home/allstaff/d/dong.x
fq=/stornext/Projects/promethion/promethion_access/lab_ritchie/A549_RNPC3/long_term/
genome=/stornext/General/data/user_managed/grpu_mritchie_1/XueyiDong/annotation/Human/CHM13v2/GCF_009914755.1_T2T-CHM13v2.0_genomic.fna
bed=/stornext/General/data/user_managed/grpu_mritchie_1/XueyiDong/annotation/Human/CHM13v2/anno.bed
# export PATH=$PATH:/wehisan/home/allstaff/d/dong.x/Programs/minimap2

module load samtools/1.7
module load minimap2/2.4

for sample in barcode{01..08} 
# do minimap2 -ax splice -uf -k14 --junc-bed  $bed $genome $fq/$sample\_pass.fq.gz | samtools view -b | samtools sort > $fq/aligned_minimap2/$sample.pass.bam
# samtools index $fq/aligned_minimap2/$sample.pass.bam $fq/aligned_minimap2/$sample.pass.bai
do minimap2 -ax splice -uf -k14 --junc-bed  $bed $genome $fq/$sample\_fail.fq.gz | samtools view -b | samtools sort > $fq/aligned_minimap2/$sample.fail.bam
samtools index $fq/aligned_minimap2/$sample.fail.bam $fq/aligned_minimap2/$sample.fail.bai
done