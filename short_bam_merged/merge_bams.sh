#!/bin/bash
#SBATCH --time=24:00:00
#SBATCH --cpus-per-task=8
#SBATCH --mem=200G
#SBATCH --job-name=merge_bams
#SBATCH --mail-user=dong.x@wehi.edu.au
#SBATCH --mail-type=BEGIN,END,FAIL

DIR=/stornext/General/data/user_managed/grpu_mritchie_1/Mei/intron/bam/
cd /stornext/General/data/user_managed/grpu_mritchie_1/XueyiDong/RNPC3/short_bam_merged
module load samtools
samtools merge --threads 7 -o merged_NT.bam --threads 7 $DIR/NT_siRNA_1_HFWKTDRXY_CGTTATTCTA-AACCTTATGG_L001_R1.fastq.gz.subread.bam $DIR/NT_siRNA_2_HFWKTDRXY_AGATCCATTA-TGGTAGAGAT_L001_R1.fastq.gz.subread.bam $DIR/NT_siRNA_3_HFWKTDRXY_GTCCTGGATA-TTCGCCACCG_L001_R1.fastq.gz.subread.bam $DIR/NT_siRNA_4_HFWKTDRXY_CAGTGGCACT-CCTATTGTTA_L001_R1.fastq.gz.subread.bam
samtools merge --threads 7 -o merged_RNPC3.bam  --threads 7 $DIR/RNPC3_siRNA_18_1_HFWKTDRXY_AGTGTTGCAC-CGTGTACCAG_L001_R1.fastq.gz.subread.bam $DIR/RNPC3_siRNA_18_2_HFWKTDRXY_GACACCATGT-TACACGTTGA_L001_R1.fastq.gz.subread.bam $DIR/RNPC3_siRNA_18_3_HFWKTDRXY_CCTGTCTGTC-TCACAACAGT_L001_R1.fastq.gz.subread.bam $DIR/RNPC3_siRNA_18_4_HFWKTDRXY_TGATGTAAGA-AAGGACGCAC_L001_R1.fastq.gz.subread.bam