# target: add gene row to gtf format
# gff containing gene records transformed by gffread
# command: gffread -E --keep-genes -F extended_annotations.gtf -o extended_annotations.gff

gff <- "../bambu/out/extended_annotations.gff"
gtf <- "../bambu/out/extended_annotations.gtf"
# gff: with gene row
gff <- read.table(gff)
# gtf: without gene row
gtf <- read.table(gtf, sep = "\t")

info <- limma::strsplit2(gtf$V9, ";")
# make geneinfo df
idx = which(!duplicated(info[,1]))
gene <- gff[gff$V3=="gene", ]
gene$V9 <- paste0(info[idx, 1], ";")
# calc nrow per gene
idx2 = c(idx[-1], nrow(info)+1)
gene.row = idx2 - idx
out <- rbind(gtf, gene)
# calculate order
ord.other = rep(1:length(idx), gene.row) + 1:nrow(info)
ord.gene = 1:length(idx) + idx -1
out <- out[order(c(ord.other, ord.gene)),]
write.table(out, "../bambu/out/extended_annotations.genes.gtf",
            sep = "\t",
            quote = FALSE,
            row.names = FALSE, col.names = FALSE)
