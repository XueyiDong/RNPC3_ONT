#!/bin/bash
#SBATCH --time=24:00:00
#SBATCH --cpus-per-task=8
#SBATCH --mem=100G
#SBATCH --mail-user=dong.x@wehi.edu.au
#SBATCH --mail-type=END,FAIL

cd /stornext/General/data/user_managed/grpu_mritchie_1/XueyiDong/RNPC3/tx_aligned_bambu_minimap2
unixhome=/wehisan/home/allstaff/d/dong.x
FQ=/stornext/Projects/promethion/promethion_access/lab_ritchie/A549_RNPC3/long_term/
ref=../salmon/bambu_transcripts.fa
export PATH=$PATH:/wehisan/home/allstaff/d/dong.x/Programs/minimap2

module load samtools/1.7

# for sample in barcode{01..08}
# do minimap2 -ax map-ont -t 8 --sam-hit-only $ref $FQ/$sample.fq.gz | samtools view -b | samtools sort > $sample.bam
# samtools index $sample.bam $sample.bai
# done

minimap2 -ax map-ont -t 8 --sam-hit-only $ref $FQ/$1.fq.gz | samtools view -b | samtools sort > $1.bam
samtools index $1.bam $.bai