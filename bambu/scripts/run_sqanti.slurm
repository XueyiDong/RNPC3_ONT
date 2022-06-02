#!/bin/bash
#SBATCH --time=48:00:00
#SBATCH --cpus-per-task=16
#SBATCH --mem=100G
#SBATCH --job-name=SQANTI3_RNPC3
#SBATCH --output=sqanti_qc_RNPC3.out
#SBATCH --mail-user=dong.x@wehi.edu.au
#SBATCH --mail-type=BEGIN,END,FAIL

# Author: Mei Du, Xueyi Dong
# Using SQANTI3 version 5.0.0 (https://github.com/ConesaLab/SQANTI3)

module load anaconda3

OUT_DIR="/stornext/General/data/user_managed/grpu_mritchie_1/XueyiDong/RNPC3/bambu/sqanti_out"
BAMUB_OUT="/stornext/General/data/user_managed/grpu_mritchie_1/XueyiDong/RNPC3/bambu/out"
ANNO_DIR="/stornext/General/data/user_managed/grpu_mritchie_1/XueyiDong/annotation/Human/CHM13v2/"

cd ${OUT_DIR}
source activate SQANTI3.env

export PYTHONPATH=$PYTHONPATH:~/Programs/cDNA_Cupcake/sequence/
export PYTHONPATH=$PYTHONPATH:~/Programs/cDNA_Cupcake/

python ~/Programs/SQANTI3-5.0/sqanti3_qc.py \
$BAMUB_OUT/extended_annotations.gtf \
$ANNO_DIR/GCF_009914755.1_T2T-CHM13v2.0_genomic.gtf \
$ANNO_DIR/GCF_009914755.1_T2T-CHM13v2.0_genomic.fna \
--cpus 16 -c ./STAR_mapping > sqanti3_qc.out  
#--short_reads ../scripts/short_read.fofn > sqanti3_qc.out

