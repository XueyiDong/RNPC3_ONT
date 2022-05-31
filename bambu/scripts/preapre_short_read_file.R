dir = "/stornext/General/data/user_managed/grpu_mritchie_1/JoanHeath/AGRF_CAGRF21056775_HFWKTDRXY"
files <- list.files(dir, "fastq.gz$")
fofn <- data.frame(
  r1 = file.path(dir, files[seq(1, 15, 2)]),
  r2 = file.path(dir, files[seq(2, 16, 2)])
)
write.table(fofn, "short_read.fofn", sep = " ", 
            col.names = FALSE, row.names = FALSE,
            quote = FALSE)
