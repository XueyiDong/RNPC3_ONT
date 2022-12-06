# bambu pipeline

# BiocManager::install(c("bambu", "Rsamtools", "BSgenome"))
BiocManager::install("bambu", version = "3.12")
library(Rsamtools)
library(bambu)
library(BSgenome)

# genome
# create fai file first with samtools
CHM13v2_seq = FaFile("/stornext/General/data/user_managed/grpu_mritchie_1/XueyiDong/annotation/Human/CHM13v2/GCF_009914755.1_T2T-CHM13v2.0_genomic.fna")

# annotations
CHM13v2_gtf = prepareAnnotations("/stornext/General/data/user_managed/grpu_mritchie_1/XueyiDong/annotation/Human/CHM13v2/GCF_009914755.1_T2T-CHM13v2.0_genomic.gtf")

# BAM files - ONT
bam_loc <- "/stornext/General/data/user_managed/grpu_mritchie_1/XueyiDong/RNPC3/aligned_minimap2"
bam_file <- list.files(bam_loc, pattern = "bam")
# bams <- sapply(file.path(bam_loc, bam_file), function(x){
#   BamFile(x)
# })
# merged_ont <- BamFileList(bams)

ont_iso_analysis = bambu(reads = file.path(bam_loc, bam_file), 
                         annotations = CHM13v2_gtf, genome = CHM13v2_seq, 
                         ncore = 16, rcOutDir = "/vast/scratch/users/du.m/bambu")
saveRDS(ont_iso_analysis, "/vast/scratch/users/du.m/bambu/bambu_out.RDS")
writeBambuOutput(ont_iso_analysis, path = "/vast/scratch/users/du.m/bambu")
sessionInfo()