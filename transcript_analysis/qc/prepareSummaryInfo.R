library(data.table)
# Files <- c("../../data/20220411_0251_2-A3-D3_PAH98828_35e1272e/sequencing_summary_PAH98828_cc5dea17.txt.gz",
#            "../../data/20220411_0423_1-A5-D5_PAI06629_97116803/sequencing_summary_PAI06629_a42e68a2.txt.gz",
#            "../../data/20220411_0423_2-E7-H7_PAI05373_1e0ae7af/sequencing_summary_PAI05373_77f46716.txt.gz")
setwd("/stornext/General/data/user_managed/grpu_mritchie_1/XueyiDong/RNPC3/transcript_analysis/qc")
Files <- c("sequencing_summary_PAH98828_cc5dea17.txt",
           "sequencing_summary_PAI06629_a42e68a2.txt",
           "sequencing_summary_PAI05373_77f46716.txt")
for (f in Files){
  Table <- fread(f, header = TRUE, sep = "\t")
  Read_Id <- as.character(Table$read_id)
  Read_length <- as.numeric(Table$sequence_length_template)
  Qscore <- as.numeric(Table$mean_qscore_template)
  Pass_filtering <- Table$passes_filtering
  Barcode <- as.character(Table$barcode_arrangement)
  Table <- cbind(Read_Id, Read_length, Qscore, Pass_filtering, Barcode)
  saveRDS(Table, paste0("rds/summaryInfo", substr(f, 20, 27), ".RDS"))
}
rds <- list.files("rds")
summaryInfo <- data.frame()
for (i in 1:length(rds)){
  tmp <- readRDS(file.path("rds", rds[i]))
  summaryInfo <- rbind(summaryInfo, tmp)
  rm(tmp)
}
saveRDS(summaryInfo, "rds/summaryInfoRDS")