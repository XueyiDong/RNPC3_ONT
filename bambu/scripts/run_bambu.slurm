#!/bin/bash
#SBATCH --time=48:00:00
#SBATCH --cpus-per-task=16
#SBATCH --mem=600G
#SBATCH --job-name=ONT_bambu
#SBATCH --output=ONT_bambu.out
#SBATCH --mail-user=du.m@wehi.edu.au
#SBATCH --mail-type=END,FAIL

# Author: Mei Du, Xueyi Dong
# Using bambu version 1.0.3 (https://github.com/GoekeLab/bambu)


module load R/4.0.5
Rscript run_bambu.R
