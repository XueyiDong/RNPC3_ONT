#!/bin/bash
#SBATCH --time=2:00:00
#SBATCH --cpus-per-task=12
#SBATCH --mem=50G
#SBATCH --mail-user=dong.x@wehi.edu.au
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --output=build_index.log

# module load anaconda3
# source activate SQANTI3.env
export PATH=/home/users/allstaff/dong.x/Programs/salmon/salmon-1.9.0_linux_x86_64/bin:$PATH
cd /stornext/General/data/user_managed/grpu_mritchie_1/XueyiDong/RNPC3/salmon
# gffread -w bambu_transcripts.fa -g /stornext/General/data/user_managed/grpu_mritchie_1/XueyiDong/annotation/Human/CHM13v2/GCF_009914755.1_T2T-CHM13v2.0_genomic.fna ../bambu/out/extended_annotations.gtf
# grep "^>" /stornext/General/data/user_managed/grpu_mritchie_1/XueyiDong/annotation/Human/CHM13v2/GCF_009914755.1_T2T-CHM13v2.0_genomic.fna | cut -d " " -f 1 > decoys.txt
# sed -i.bak -e 's/>//g' decoys.txt
# cat bambu_transcripts.fa /stornext/General/data/user_managed/grpu_mritchie_1/XueyiDong/annotation/Human/CHM13v2/GCF_009914755.1_T2T-CHM13v2.0_genomic.fna > gentrome.fa
salmon index -t gentrome.fa -d decoys.txt -p 12 -i salmon_index