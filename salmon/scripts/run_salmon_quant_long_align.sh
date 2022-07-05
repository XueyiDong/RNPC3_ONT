#!/bin/bash
#SBATCH --time=48:00:00
#SBATCH --cpus-per-task=16
#SBATCH --mem=50G
#SBATCH --mail-user=dong.x@wehi.edu.au
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --output=salmon_quant_long_bambu.log

# this is the final version used

cd /stornext/General/data/user_managed/grpu_mritchie_1/XueyiDong/RNPC3/salmon
module load samtools

export PATH=/home/users/allstaff/dong.x/Programs/salmon/salmon-1.9.0_linux_x86_64/bin:$PATH

FQ=/stornext/Projects/promethion/promethion_access/lab_ritchie/A549_RNPC3/long_term/
mkdir -p ../tx_aligned_bambu_minimap2/primary

# for SAMPLE in barcode{01..08}
# do
# 	mkdir -p ./long_bambu/$SAMPLE
# 	samtools view -F 256 -b -@ 8 ../tx_aligned_bambu_minimap2/$SAMPLE\_pass.bam > ../tx_aligned_bambu_minimap2/primary/$SAMPLE\_pass.bam
# 	samtools view -F 256 -b -@ 8 ../tx_aligned_bambu_minimap2/$SAMPLE\_fail.bam > ../tx_aligned_bambu_minimap2/primary/$SAMPLE\_fail.bam
# 	salmon quant -t bambu_transcripts.fa -l A -a ../tx_aligned_bambu_minimap2/primary/$SAMPLE\_pass.bam ../tx_aligned_bambu_minimap2/primary/$SAMPLE\_fail.bam -o ./long_bambu/$SAMPLE -p 16 --numBootstraps 100 -v
# done

mkdir -p ./long_bambu/$1
samtools view -F 256 -b -@ 8 ../tx_aligned_bambu_minimap2/$1\_pass.bam > ../tx_aligned_bambu_minimap2/primary/$1\_pass.bam
samtools view -F 256 -b -@ 8 ../tx_aligned_bambu_minimap2/$1\_fail.bam > ../tx_aligned_bambu_minimap2/primary/$1\_fail.bam
salmon quant -t bambu_transcripts.fa -l A -a ../tx_aligned_bambu_minimap2/primary/$1\_pass.bam ../tx_aligned_bambu_minimap2/primary/$1\_fail.bam -o ./long_bambu/$1 -p 16 --numBootstraps 100 -v