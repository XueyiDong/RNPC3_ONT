#!/bin/bash
#SBATCH --time=48:00:00
#SBATCH --cpus-per-task=16
#SBATCH --mem=600G
#SBATCH --job-name=ONT_bambu
#SBATCH --output=ONT_bambu.out
#SBATCH --mail-user=dong.x@wehi.edu.au
#SBATCH --mail-type=BEGIN,END,FAIL

# Author: Mei Du, Xueyi Dong
# Using bambu version 1.0.3 (https://github.com/GoekeLab/bambu)

cd /stornext/General/data/user_managed/grpu_mritchie_1/XueyiDong/RNPC3/bambu/outputs/

module load R/4.1.1
Rscript /stornext/General/data/user_managed/grpu_mritchie_1/XueyiDong/RNPC3/bambu/scripts/run_bambu.R
