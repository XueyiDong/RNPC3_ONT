#!/bin/bash
#SBATCH --time=48:00:00
#SBATCH --cpus-per-task=16
#SBATCH --mem=50G
#SBATCH --mail-user=dong.x@wehi.edu.au
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --output=salmon_quant_short_bambu.log

cd /stornext/General/data/user_managed/grpu_mritchie_1/XueyiDong/RNPC3/salmon

export PATH=/home/users/allstaff/dong.x/Programs/salmon/salmon-1.9.0_linux_x86_64/bin:$PATH

FQ=/stornext/General/data/user_managed/grpu_mritchie_1/JoanHeath/AGRF_CAGRF21056775_HFWKTDRXY

for SAMPLE in NT_siRNA_1_HFWKTDRXY_CGTTATTCTA-AACCTTATGG NT_siRNA_2_HFWKTDRXY_AGATCCATTA-TGGTAGAGAT NT_siRNA_3_HFWKTDRXY_GTCCTGGATA-TTCGCCACCG NT_siRNA_4_HFWKTDRXY_CAGTGGCACT-CCTATTGTTA RNPC3_siRNA_18_1_HFWKTDRXY_AGTGTTGCAC-CGTGTACCAG RNPC3_siRNA_18_2_HFWKTDRXY_GACACCATGT-TACACGTTGA RNPC3_siRNA_18_3_HFWKTDRXY_CCTGTCTGTC-TCACAACAGT RNPC3_siRNA_18_4_HFWKTDRXY_TGATGTAAGA-AAGGACGCAC
do
	mkdir -p ./short_bambu/$SAMPLE
	salmon quant -i salmon_index -l A -1 $FQ/$SAMPLE\_L001_R1.fastq.gz -2 $FQ/$SAMPLE\_L001_R2.fastq.gz --validateMappings -o ./short_bambu/$SAMPLE -p 16 --numBootstraps 100
done