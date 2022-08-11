library(ggplot2)
library(edgeR)
library(MetBrewer)
library(scales)
library(RColorBrewer)
library(tidyverse)
library(cowplot)
library(viridis)

# DIR="/stornext/General/data/user_managed/grpu_mritchie_1/XueyiDong/long_read_benchmark"
# load DGE lists
# s <- catchSalmon(file.path(DIR, "ONT/salmon_bs", list.files(file.path(DIR, "/ONT/salmon_bs"))))
# dge <- DGEList(counts=s$counts/s$annotation$Overdispersion, genes=s$annotation)
# s.short <- catchSalmon(file.path(DIR, "illumina/salmon_bs", list.files(file.path(DIR, "illumina/salmon_bs"))))
# dge.short <- DGEList(counts = s.short$counts/s.short$annotation$Overdispersion, genes = s.short$annotation)

# dge <- readRDS("dge.rds")
# dge.short <- readRDS("dge.short.rds")

# organize long-read counts ----
counts <- read.delim("../../bambu/out/counts_transcript.txt")
rownames(counts) <- counts$TXNAME
dge <- DGEList(counts = counts[, c(-1, -2)])
dge$counts[, seq(2, 16, 2)] <-dge$counts[, seq(1, 15, 2)] + dge$counts[, seq(2, 16, 2)] 
dge$samples$lib.size[seq(2, 16, 2)] <- dge$samples$lib.size[seq(1, 15, 2)] + dge$samples$lib.size[seq(2, 16, 2)]
dge <- dge[, -seq(1, 15, 2)]
dge$samples$group <- rep(c("control", "targeted"), c(4, 4))
colnames(dge) <- gsub("_pass", "", colnames(dge))

# organize short-read counts ----
counts.short <- catchSalmon(path = paste0("../../salmon/short_bambu/", list.files("../../salmon/short_bambu/")))
# inflate by overdisp
dge.short <- DGEList(counts = counts.short$counts / counts.short$annotation$Overdispersion, genes = counts.short$annotation)
dge.short$samples$group <- rep(c("control", "targeted"), c(4, 4))
colnames(dge.short) <- paste(rep(c("NT", "RNPC3"), each = 4), rep(1:4, 2), sep = "_")

# read num plot----
# organize read num stat 
read.stat <- data.frame(
  sample = rep(paste(rep(c("NT", "RNPC3"), each = 4), rep(1:4, 2), sep = "-"), 2),
  raw_reads = c(21460386+7649032, 32984104+12383754, 27061720+10829457, 28627983+11681329, 28141261+9703764, 26939112+10246879, 29952450+11023531, 24600719+11215100,
                58963456, 40652346, 74909438, 33936801, 72561084, 75934404, 70456874, 53348132),
  # mapped_reads = c(18948958+4369497, 27968074+6417123, 22410004+5199474, 24569023+5743348, 24859748+5647440, 24050590+5700168, 26030588+5913487, 20093397+4699682,
  #                  27495741, 18453726, 34247708, 15494205, 32423391, 34238487, 32054224, 22274216),
  read_counts = c(dge$samples$lib.size, colSums(counts.short$counts)),
  dataset = rep(c("ONT", "Illumina"), c(8,8))
)

read.stat <- reshape2::melt(read.stat, id.vars = c("sample", "dataset"))
# read num plot
pdf("plots/readNum.pdf", height = 5, width = 8)
ggplot(read.stat, aes(x=variable, y=value, fill=sample, label = value))+
  geom_bar(stat="identity") +
  facet_grid(cols=vars(dataset)) +
  geom_text(aes(label = label_number_si()(value)), size = 4, colour = "white", position = position_stack(vjust = 0.5)) +
  ylab("Number") +
  xlab("Category") +
  theme_bw() +
  theme(text = element_text(size = 20), axis.text.x = element_text(angle = 30, hjust = 1)) +
  scale_y_continuous(labels = unit_format(unit = "M", scale = 1e-6))+
  scale_fill_manual(values = met.brewer("Tiepolo", 8))
dev.off()

# quantification correlation matrix heatmap by biotype-------
## use raw count -----
m <- match(rownames(dge.short), rownames(dge))
table(is.na(m))
dge.all <- DGEList(counts = cbind(dge$counts[m,], dge.short$counts))
dge.all$samples$group <- rep(c("control_long", "targeted_long", "control_short", "targeted_short"), each = 4)
filt <- filterByExpr(dge.all)
dge.all <- dge.all[filt,]
dge.all$genes <- dge.short$genes[match(rownames(dge.all), rownames(dge.short)),]
colnames(dge.all$genes)[3] <- "Overdispersion.short"
cormat <- cor(dge.all$counts)
library(pheatmap)
anno <- data.frame(
  group = rep(rep(c("NT", "RNPC3"), c(4, 4)), 2),
  dataset = rep(c("ONT", "Illumina"), c(8, 8))
)
anno_colours = list(
  group = c(NT = "#8b3a2b", RNPC3 = "#235070"),
  dataset = c(ONT = "#438DAC", Illumina = "#FCB344")
)
rownames(anno) <- rownames(cormat)
pheatmap(cormat,
         cluster_cols = FALSE,
         cluster_rows = FALSE,
         show_colnames = FALSE,
         show_rownames = FALSE,
         annotation_col = anno,
         annotation_row = anno,
         annotation_colors = anno_colours,
         scale = "none",
         display_numbers = TRUE
)


# long vs short quantification------
# long CPM vs short TPM
# filter
# cpm.long <- dge[filterByExpr(dge),] %>% calcNormFactors %>% cpm
cpm.long <- cpm(dge.all[, 1:8])
tpm3 <- function(counts,len) {
  x <- counts/len
  return(t(t(x)*1e6/colSums(x)))
}
# tpm.short <- dge.short[filterByExpr(dge.short),] %>% calcNormFactors %>% {tpm3(.[["counts"]], .[["genes"]][,"Length"])}
tpm.short <- tpm3(dge.all$counts[,9:16], dge.all$genes$Length)
# m <- match(rownames(dge.short), rownames(dge))
quant <- data.frame(
  TPM_short = c(log(rowMeans(tpm.short[, 1:4]) + 0.5),
                log(rowMeans(tpm.short[, 5:8]) + 0.5)),
  CPM_long = c(log(rowMeans(cpm.long[, 1:4]) + 0.5),
               log(rowMeans(cpm.long[, 5:8]) + 0.5)),
  group = rep(c("NT", "RNPC3"), rep(nrow(tpm.short), 2))
)
# quant <- na.omit(quant)

# correlation
cor(quant$TPM_short, quant$CPM_long)
cor(quant$TPM_short[quant$group=="NT"], quant$CPM_long[quant$group=="NT"])
cor(quant$TPM_short[quant$group=="RNPC3"], quant$CPM_long[quant$group=="RNPC3"])
pdf("plots/longVsShortQuant.pdf", height = 5, width = 8)
ggplot(quant, aes(x = CPM_long, y = TPM_short))+
  stat_binhex(bins=100) +
  scale_fill_viridis(trans = "log10", option = "A")+
  annotate(geom="text", x=3, y=12,
           label=paste0("Pearson's r=", round(cor(quant$TPM_short, quant$CPM_long), 3)),
           size=7)+
  labs(x = expression("log"[2]*"CPM ONT read counts"),
       y = expression("log"[2]*"TPM Illumina read counts")
  ) +
  theme_bw() +
  theme(text=element_text(size = 20))
dev.off()

# DISREGARD BELOW PARTS! ----

# biotype ONT----
# gene biotype info
dge.human <- dge[grep("^ENST", rownames(dge)), ]
txid <- strsplit2(rownames(dge.human$counts), "\\|")[,1]
txid <- strsplit2(txid, "\\.")[,1]

library(AnnotationHub)
ah <- AnnotationHub()
EnsDb.Hsapiens.v104 <- query(ah, c("EnsDb", "Homo Sapiens", 104))[[1]]
biotype<- mapIds(
  x = EnsDb.Hsapiens.v104,
  # NOTE: Need to remove gene version number prior to lookup.
  keys = txid,
  keytype = "TXID",
  column = "TXBIOTYPE")
dge.human$genes$biotype <- biotype
# deal with biotype
# http://asia.ensembl.org/info/genome/genebuild/biotypes.html
dge.human$genes$biotype[grepl("pseudogene$", dge.human$genes$biotype)] <- "pseudogene"
dge.human$genes$biotype[grepl("^TR", dge.human$genes$biotype)] <- "IG_or_TR_gene"
dge.human$genes$biotype[grepl("^IG", dge.human$genes$biotype)] <- "IG_or_TR_gene"
dge.human$genes$biotype[dge.human$genes$biotype %in% c("miRNA", "misc_RNA", 
                                                       "piRNA", "rRNA", "siRNA",
                                                       "snRNA", "snoRNA", "scaRNA",
                                                       "tRNA", "vault_RNA", "scRNA",
                                                       "sRNA", "Mt_rRNA", "Mt_tRNA",
                                                       "ribozyme"
)] <- "ncRNA"
saveRDS(dge.human$genes, "txInfo.long.RDS")

# for each sample
biotype_sum <- sapply(1:6, function(x){
  typesum = aggregate(dge.human$counts[,x], by=list(dge.human$genes$biotype), FUN=sum, simplify=TRUE)
  return(typesum)
}, simplify=FALSE)
biotype_sum <- do.call("rbind", biotype_sum)
biotype_sum$sample <- rep(c("H1975-1", "H1975-2", "H1975-3", "HCC827-1", "HCC827-2", "HCC827-5"), rep(10, 6))
colnames(biotype_sum) <- c("biotype", "total_count", "sample")
# order the bars
ord = aggregate(biotype_sum$total_count, by = list(biotype_sum$biotype), FUN = sum, simplify = TRUE)
ord = ord[order(ord$x), ]

pdf("plots/biotype.pdf", height = 5, width = 8)
ggplot(biotype_sum, aes(x=sample, y=total_count, fill=factor(biotype, levels=ord$Group.1))) +
  geom_bar(stat="identity", position = "fill") +
  theme_bw() +
  theme(text = element_text(size = 20), axis.text.x = element_text(angle = 90, vjust = 0.5)) +
  scale_fill_brewer(palette = "Set3") +
  labs(fill = "Transcript biotype", x = "Sample", y = "Proportion of count")
dev.off()

# calculate proportion for ONT----
biotype_sum$proportion <- sapply(1:nrow(biotype_sum), function(x){
  sample.sum = sum(biotype_sum$total_count[biotype_sum$sample == biotype_sum$sample[x]])
  return(biotype_sum$total_count[x] / sample.sum)
}, simplify = TRUE)

#biotype Illumina---- 
dge.short.human <- dge[grep("^ENST", rownames(dge.short)), ]
txid <- strsplit2(rownames(dge.short.human$counts), "\\|")[,1]
txid <- strsplit2(txid, "\\.")[,1]

biotype<- mapIds(
  x = EnsDb.Hsapiens.v104,
  # NOTE: Need to remove gene version number prior to lookup.
  keys = txid,
  keytype = "TXID",
  column = "TXBIOTYPE")
dge.short.human$genes$biotype <- biotype
saveRDS(biotype, "biotype.RDS")
# deal with biotype
dge.short.human$genes$biotype[grepl("pseudogene$", dge.short.human$genes$biotype)] <- "pseudogene"
dge.short.human$genes$biotype[grepl("^TR", dge.short.human$genes$biotype)] <- "IG_or_TR_gene"
dge.short.human$genes$biotype[grepl("^IG", dge.short.human$genes$biotype)] <- "IG_or_TR_gene"
dge.short.human$genes$biotype[dge.short.human$genes$biotype %in% c("miRNA", "misc_RNA", 
                                                                   "piRNA", "rRNA", "siRNA",
                                                                   "snRNA", "snoRNA", "scaRNA",
                                                                   "tRNA", "vault_RNA", "scRNA",
                                                                   "sRNA", "Mt_rRNA", "Mt_tRNA",
                                                                   "ribozyme"
)] <- "ncRNA"
saveRDS(dge.short.human$genes, "txInfo.short.RDS")

# for each sample
biotype_sum.short <- sapply(1:6, function(x){
  typesum = aggregate(dge.short.human$counts[,x], by=list(dge.short.human$genes$biotype), FUN=sum, simplify=TRUE)
  return(typesum)
}, simplify=FALSE)
biotype_sum.short <- do.call("rbind", biotype_sum.short)

biotype_sum.short$sample <- rep(c("H1975-1", "H1975-2", "H1975-3", "HCC827-1", "HCC827-2", "HCC827-5"), rep(10, 6))
# biotype_sum.short$sample <- rep(paste0("barcode0", 1:6), rep(10, 6))
colnames(biotype_sum.short) <- c("biotype", "total_count", "sample")
# order the bars
ord = aggregate(biotype_sum$total_count, by = list(biotype_sum$biotype), FUN = sum, simplify = TRUE)
ord = ord[order(ord$x), ]

pdf("plots/biotype_short.pdf", height = 5, width = 8)
ggplot(biotype_sum.short, aes(x=sample, y=total_count, fill=factor(biotype, levels=ord$Group.1))) +
  geom_bar(stat="identity", position = "fill") +
  theme_bw() +
  theme(text = element_text(size = 20), axis.text.x = element_text(angle = 90, vjust = 0.5)) +
  scale_fill_brewer(palette = "Set3") +
  labs(fill = "Transcript biotype", x = "Sample", y = "Proportion of count")
dev.off()

# calculate proportion for Illumina----
biotype_sum.short$proportion <- sapply(1:nrow(biotype_sum.short), function(x){
  sample.sum = sum(biotype_sum.short$total_count[biotype_sum.short$sample == biotype_sum.short$sample[x]])
  return(biotype_sum.short$total_count[x] / sample.sum)
}, simplify = TRUE)

#biotype long and short---- 
biotype_sum.all <- rbind(biotype_sum, biotype_sum.short)
biotype_sum.all$dataset <- rep(c("ONT", "Illumina"), c(nrow(biotype_sum), nrow(biotype_sum.short)))

pdf("plots/biotype_all.pdf", height = 4, width = 8)
ggplot(biotype_sum.all, aes(x=sample, y=total_count, fill=factor(biotype, levels=ord$Group.1))) +
  geom_bar(stat="identity", position = "fill") +
  facet_grid(cols=vars(dataset)) +
  theme_bw() +
  theme(text = element_text(size = 20), 
        axis.text.x = element_text(angle = 45, hjust = 1),
  ) +
  scale_fill_brewer(palette = "Set3") +
  labs(fill = "Transcript biotype", x = "Sample", y = "Proportion of count")
dev.off()




# quantification correlation matrix heatmap by biotype-------
## use long cpm vs short tpm----
### all----
cpm <- cpm(dge.all[, 1:6])
tpm <- tpm3(dge.all$counts[, 7:12], dge.short$genes$Length[match(rownames(dge.all), rownames(dge.short))])
quant.all <- cbind(cpm, tpm)
cormat2 <- cor(quant.all)
pdf("plots/corHeatmap.pdf", height = 8, width = 9)
pheatmap(cormat2,
         color = colorRampPalette(brewer.pal(n = 7, name = "PuBuGn"))(100),
         cluster_cols = FALSE,
         cluster_rows = FALSE,
         show_colnames = FALSE,
         show_rownames = FALSE,
         annotation_col = anno,
         annotation_row = anno,
         annotation_colors = anno_colours,
         scale = "none",
         display_numbers = TRUE,
         number_color = "black",
         fontsize = 16,
         fontsize_number = 12
)
dev.off()
### coding----
# filter by biotype
tx.sel <- names(biotype)[biotype %in% c("protein_coding", "lncRNA")]
dge.pure.sel <- dge.pure[substr(rownames(dge.pure), 1, 15) %in% tx.sel, ]
dge.short.pure.sel <- dge.short.pure[substr(rownames(dge.short.pure), 1, 15) %in% tx.sel, ]
cpm.sel <- cpm(dge.pure.sel)
tpm.sel <- tpm3(dge.short.pure.sel$counts, dge.short.pure.sel$genes$Length)
quant.all.sel <- cbind(cpm.sel[match(rownames(tpm.sel), rownames(cpm.sel)), ], tpm.sel) %>% na.omit
cormat3 <- cor(quant.all.sel)
pdf("plots/corHeatmapCoding.pdf", height = 7, width = 8)
pheatmap(cormat3,
         color = colorRampPalette(brewer.pal(n = 7, name = "PuBuGn"))(100),
         cluster_cols = FALSE,
         cluster_rows = FALSE,
         show_colnames = FALSE,
         show_rownames = FALSE,
         annotation_col = anno,
         annotation_row = anno,
         annotation_colors = anno_colours,
         scale = "none",
         display_numbers = TRUE,
         number_color = "black",
         fontsize = 16,
         fontsize_number = 12
)
dev.off()

### other----
tx.sel <- names(biotype)[!(biotype %in% c("protein_coding", "lncRNA"))]
dge.pure.sel <- dge.pure[substr(rownames(dge.pure), 1, 15) %in% tx.sel, ]
dge.short.pure.sel <- dge.short.pure[substr(rownames(dge.short.pure), 1, 15) %in% tx.sel, ]
cpm.sel <- cpm(dge.pure.sel)
tpm.sel <- tpm3(dge.short.pure.sel$counts, dge.short.pure.sel$genes$Length)
quant.all.sel <- cbind(cpm.sel[match(rownames(tpm.sel), rownames(cpm.sel)), ], tpm.sel) %>% na.omit
cormat3 <- cor(quant.all.sel)
pdf("plots/corHeatmapOther.pdf", height = 8, width = 9)
pheatmap(cormat3,
         color = colorRampPalette(brewer.pal(n = 7, name = "PuBuGn"))(100),
         cluster_cols = FALSE,
         cluster_rows = FALSE,
         show_colnames = FALSE,
         show_rownames = FALSE,
         annotation_col = anno,
         annotation_row = anno,
         annotation_colors = anno_colours,
         scale = "none",
         display_numbers = TRUE,
         number_color = "black",
         fontsize = 16,
         fontsize_number = 12
)
dev.off()

### sequins----
cormat4 <- cor(quant.all[grepl("^R", rownames(quant.all)),])
pdf("plots/corHeatmapSequin.pdf", height = 8, width = 9)
pheatmap(cormat4,
         color = colorRampPalette(brewer.pal(n = 7, name = "PuBuGn"))(100),
         cluster_cols = FALSE,
         cluster_rows = FALSE,
         show_colnames = FALSE,
         show_rownames = FALSE,
         annotation_col = anno,
         annotation_row = anno,
         annotation_colors = anno_colours,
         scale = "none",
         display_numbers = TRUE,
         number_color = "black",
         fontsize = 16,
         fontsize_number = 12
)
dev.off()

## by number of tx-------
### 1----
cormat.tx1 <- cor(quant.all[dge.all$genes$nTranscript==1,])
pdf("plots/corHeatmaptx1.pdf", height = 8, width = 9)
pheatmap(cormat.tx1,
         color = colorRampPalette(brewer.pal(n = 7, name = "PuBuGn"))(100),
         cluster_cols = FALSE,
         cluster_rows = FALSE,
         show_colnames = FALSE,
         show_rownames = FALSE,
         annotation_col = anno,
         annotation_row = anno,
         annotation_colors = anno_colours,
         scale = "none",
         display_numbers = TRUE,
         number_color = "black",
         fontsize = 16,
         fontsize_number = 12
)
dev.off()

### 2 - 5----
cormat.tx2_5 <- cor(quant.all[dge.all$genes$nTranscript>=2 & dge.all$genes$nTranscript <= 5,])
pdf("plots/corHeatmaptx2_5.pdf", height = 8, width = 9)
pheatmap(cormat.tx2_5,
         color = colorRampPalette(brewer.pal(n = 7, name = "PuBuGn"))(100),
         cluster_cols = FALSE,
         cluster_rows = FALSE,
         show_colnames = FALSE,
         show_rownames = FALSE,
         annotation_col = anno,
         annotation_row = anno,
         annotation_colors = anno_colours,
         scale = "none",
         display_numbers = TRUE,
         number_color = "black",
         fontsize = 16,
         fontsize_number = 12
)
dev.off()

### > 5----
cormat.tx5 <- cor(quant.all[dge.all$genes$nTranscript > 5,])
pdf("plots/corHeatmaptx5.pdf", height = 8, width = 9)
pheatmap(cormat.tx5,
         color = colorRampPalette(brewer.pal(n = 7, name = "PuBuGn"))(100),
         cluster_cols = FALSE,
         cluster_rows = FALSE,
         show_colnames = FALSE,
         show_rownames = FALSE,
         annotation_col = anno,
         annotation_row = anno,
         annotation_colors = anno_colours,
         scale = "none",
         display_numbers = TRUE,
         number_color = "black",
         fontsize = 16,
         fontsize_number = 12
)
dev.off()





# length bias plot--------------------
dge.pure$genes$totalCount <- rowSums(dge.pure$counts)
lm.long <- lm(dge.pure$genes$totalCount~dge.pure$genes$Length)
summary(lm.long)
# Multiple R-squared:  0.0007356,	Adjusted R-squared:  0.0007191 
pdf("plots/longLenBias.pdf", height = 5, width = 8)
ggplot(dge.pure$genes, aes(x = Length, y = totalCount))+
  stat_binhex() +
  geom_smooth(formula = y~x, method="lm") +
  scale_fill_viridis(trans = "log10", option = "A")+
  scale_x_continuous(trans = "log10") +
  scale_y_continuous(trans = "log10") +
  annotate(geom="text", x=max(dge.pure$genes$Length) * 0.1, y=max(dge.pure$genes$totalCount) * 0.8,
           label=paste0("Pearson's r=", round(cor(dge.pure$genes$Length, dge.pure$genes$totalCount), 3)),
           size=7)+
  labs(x = "Transcript length",
       y = "Total read count"
  ) +
  theme_bw() +
  theme(text=element_text(size = 20))
dev.off()

p <- lapply(1:6, function(x){
  dat <- data.frame(
    Length = dge.pure$genes$Length,
    Count = dge.pure$counts[,x]
  )
  p <- ggplot(dat, aes(x = Length, y = Count))+
    stat_binhex() +
    geom_smooth(formula = y~x, method="lm") +
    scale_fill_viridis(trans = "log10", option = "A")+
    scale_x_continuous(trans = "log10") +
    scale_y_continuous(trans = "log10") +
    annotate(geom="text", x=max(dat$Length) * 0.05, y=max(dat$Count) * 0.8,
             label=paste0("Pearson's r=", round(cor(dat$Length, dat$Count), 3)),
             size=7)+
    labs(x = "Transcript length",
         y = "Total read count"
    ) +
    theme_bw() +
    theme(text=element_text(size = 20), axis.text.x = element_text(angle = 30, hjust = 1))
  return(p)
})
pdf("plots/longLenBiasSamp.pdf", height = 9, width = 16)
plot_grid(p[[1]], p[[2]], p[[3]], p[[4]], p[[5]], p[[6]],
          labels = paste(rep(c("H1975", "HCC827"), c(3, 3)), c(1, 2, 3, 1, 2, 5), sep = "-"))
dev.off()

dge.short.pure$genes$totalCount <- rowSums(dge.short.pure$counts)
lm.short <- lm(dge.short.pure$genes$totalCount~dge.short.pure$genes$Length)
summary(lm.short)
# Multiple R-squared:  3.275e-06,	Adjusted R-squared:  -2.818e-05
pdf("plots/shortLenBias.pdf", height = 5, width = 8)
ggplot(dge.short.pure$genes, aes(x = Length, y = totalCount))+
  stat_binhex() +
  geom_smooth(formula = y~x, method="lm") +
  scale_fill_viridis(trans = "log10", option = "A")+
  scale_x_continuous(trans = "log10") +
  scale_y_continuous(trans = "log10") +
  annotate(geom="text", x=max(dge.short.pure$genes$Length) * 0.1, y=max(dge.short.pure$genes$totalCount) * 0.8,
           label=paste0("Pearson's r=", round(cor(dge.short.pure$genes$Length, dge.short.pure$genes$totalCount), 3)),
           size=7)+
  labs(x = "Transcript length",
       y = "Total read count"
  ) +
  theme_bw() +
  theme(text=element_text(size = 20))
dev.off()

p <- lapply(1:6, function(x){
  dat <- data.frame(
    Length = dge.short.pure$genes$Length,
    Count = dge.short.pure$counts[,x]
  )
  p <- ggplot(dat, aes(x = Length, y = Count))+
    stat_binhex() +
    geom_smooth(formula = y~x, method="lm") +
    scale_fill_viridis(trans = "log10", option = "A")+
    scale_x_continuous(trans = "log10") +
    scale_y_continuous(trans = "log10") +
    annotate(geom="text", x=max(dat$Length) * 0.05, y=max(dat$Count) * 0.8,
             label=paste0("Pearson's r=", round(cor(dat$Length, dat$Count), 3)),
             size=7)+
    labs(x = "Transcript length",
         y = "Total read count"
    ) +
    theme_bw() +
    theme(text=element_text(size = 20), axis.text.x = element_text(angle = 30, hjust = 1))
  return(p)
})
pdf("plots/ShortLenBiasSamp.pdf", height = 9, width = 16)
plot_grid(p[[1]], p[[2]], p[[3]], p[[4]], p[[5]], p[[6]],
          labels = paste(rep(c("H1975", "HCC827"), c(3, 3)), c(1, 2, 3, 1, 2, 5), sep = "-"))
dev.off()