# bambu pipeline

BiocManager::install(c("bambu", "Rsamtools", "BSGenome"))
library(Rsamtools)
library(bambu)
library(BSgenome)

# genome
# create fai file first with samtools
CHM13v2_seq = FaFile("/stornext/General/data/user_managed/grpu_mritchie_1/XueyiDong/annotation/Human/CHM13v2/GCF_009914755.1_T2T-CHM13v2.0_genomic.fna")

# annotations
CHM13v2_gtf = prepareAnnotations("/stornext/General/data/user_managed/grpu_mritchie_1/XueyiDong/annotation/Human/CHM13v2/GCF_009914755.1_T2T-CHM13v2.0_genomic.gtf")

# BAM files - ONT
bam_loc <- "/stornext/Projects/promethion/promethion_access/lab_ritchie/A549_RNPC3/long_term/aligned_minimap2"
bam_file <- list.files(bam_loc, pattern = "bam")
bams <- sapply(file.path(bam_loc, bam_file), function(x){
  BamFile(x, index=gsub("bam", "bai", x))
})
merged_ont <- BamFileList(bams)

ont_iso_analysis = bambu(reads = merged_ont, annotations = CHM13v2_gtf, genome = CHM13v2_seq, ncore = 8)
saveRDS(ont_iso_analysis, "/stornext/General/data/user_managed/grpu_mritchie_1/XueyiDong/RNPC3/bambu/outputs/bambu_out.RDS")
writeBambuOutput(ont_iso_analysis, path = "/stornext/General/data/user_managed/grpu_mritchie_1/XueyiDong/RNPC3/bambu/outputs", prefix = "ONT_")