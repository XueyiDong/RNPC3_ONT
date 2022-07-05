# TO DO LIST ----
# 1. make heatmaps
# 2. add option in function for adding alignment track
# 3. Select Stephen's highlighted genes

## DONE LIST ----
# 1. add short read coverage
# 2. inport Stephen's gene lists

# Start up ----
library(superintronic)
suppressPackageStartupMessages(library(plyranges))
library(ggplot2)
library(pheatmap)

# prepare annotation----
gff <- "../bambu/out/extended_annotations.genes.edited.gtf"
gr <- read_gff(gff, genome_info = "CHM13v2")
# parts <- collect_parts(gr)
# saveRDS(parts, "complete_annotation.rds")
parts <- readRDS("complete_annotation.rds")

# design for long read ----
bam_dir = "../aligned_minimap2"
design <- data.frame(
  Sample = paste0(rep(paste(rep(c("NT", "targeted"), c(4, 4)), rep(1:4, 2), sep = "_"), rep(2, 8)), rep(c("_fail", "_pass"), 8)),
  Bam = file.path(bam_dir, list.files(bam_dir, ".bam$")),
  Replicate = rep(rep(1:4, 2), rep(2, 8)),
  siRNA = rep(c("NT", "RNPC3"), c(8,8))
)
# compute coverage for long read ----
# cvg <- compute_coverage_long(design, source="Bam")
# saveRDS(cvg, "complete_coverage.rds")
cvg <- readRDS("complete_coverage.rds")

# design for short read ----
bam_short_dir = "/stornext/General/data/user_managed/grpu_mritchie_1/Mei/intron/bam"
design.short <- data.frame(
  Sample = paste0(rep(c("NT", "targeted"), each = 4), rep(1:4, 2)),
  Bam = file.path(bam_short_dir, list.files(bam_short_dir, ".bam$")),
  Replicate = rep(1:4, 2),
  siRNA = rep(c("NT", "RNPC3"), each = 4)
)
# compute coverage for short read ----
# cvg.short <- compute_coverage_long(design.short, source="Bam")
# saveRDS(cvg.short, "complete_coverage_short.RDS")
cvg.short <- readRDS("complete_coverage_short.RDS")

# transform intron coord into new version----
# write minor intron information into bed format
tmap <- read.delim("../bambu/out/bambu_comp.extended_annotations.gtf.tmap")
refmap <- read.delim("../bambu/out/bambu_comp.extended_annotations.gtf.refmap")
refmap$qry_gene_id <- limma::strsplit2(refmap$qry_id_list, "|", fixed = TRUE)[,1]
refmap$qry_id <- limma::strsplit2(refmap$qry_id_list, "|", fixed = TRUE)[,2]
mi_gene <- xlsx::read.xlsx("../HomoSapiens_IntronInfo.xlsx", sheetIndex = 1)
# mi <- GRanges(
#   seqnames = Rle(paste0("chr", mi_gene$Chromosome)),
#   ranges = IRanges(as.numeric(limma::strsplit2(mi_gene$Intron.coordinates, "-")[,1]),
#                    end = as.numeric(limma::strsplit2(mi_gene$Intron.coordinates, "-")[,2]),
#                    names = mi_gene$Gene.Name)
# )
# write_bed(mi, "mi.bed")
# #intron coord converted using online tool: https://genome.ucsc.edu/cgi-bin/hgLiftOver
mi <- read_bed("hglft_genome_b5cb_ee5900.bed")
mi$gene_id <- tmap$qry_gene_id[match(mi$name, tmap$ref_gene_id)]

# Visualize coverage using Gviz ----
library(Gviz)

make_cov_track <- function(cvg, gene, dataset){
  if (gene %in% parts$gene_id){
    # cat("Making plot for gene", gene, ".\n")
    features <- parts %>% filter(gene_id == gene)
    cvg_over_features <- cvg %>% 
      select(-Bam) %>% 
      join_parts(features)
    if(length(unique(BiocGenerics::strand(unnest_parts(features)))) > 1) {
      cat("Gene", gene, "strand not unique.\n")
    } else{
      # calculate mean coverage score by group
      data <- cvg_over_features %>% plyranges::group_by(siRNA) %>% disjoin_ranges_directed(score = mean(score))
      dTrack <- lapply(unique(data$siRNA), function(x){
        DataTrack(range = data %>% filter(siRNA == x),
                  options(ucscChromosomeNames=FALSE), 
                  name = paste(x, dataset),
                  data = "score")
      })
      return(dTrack)
    }
  }
}

plot_cov_genes2 <- function(gene, cov_data = c("long", "short"), anno_col = "transcript_id"){
    if (gene %in% parts$gene_id){
      # cat("Making plot for gene", gene, ".\n")
      features <- parts %>% filter(gene_id == gene)
      if(length(unique(BiocGenerics::strand(unnest_parts(features)))) > 1) {
        cat("Gene", gene, "strand not unique.\n")
      } else{
        # axis
        axisTrack <- GenomeAxisTrack()
        # transcripts annotation
        trackList = GeneRegionTrack(gr %>% 
                                      filter(gene_id==gene) %>% 
                                      GenomicFeatures::makeTxDbFromGRanges(),
                                    transcriptAnnotation = "transcript",
                                    options(ucscChromosomeNames=FALSE))
        # add coverage histogram track
        if("long" %in% cov_data){
          dTrack_long <- make_cov_track(cvg, gene, "long")
          trackList <- append(dTrack_long, trackList)
        }
        if("short" %in% cov_data){
          dTrack_short <- make_cov_track(cvg.short, gene, "short")
          trackList <- append(dTrack_short, trackList)
        }
        # add highlight for minor intron region
        ref_gene <- tmap$ref_gene_id[match(gene, tmap$qry_gene_id)][1]
        ht <- HighlightTrack(trackList = trackList,
                             range = mi %>% filter(name == ref_gene))
        plotTracks(list(ht, axisTrack),
                   type="h", groupAnnotation = "group",
                   main = paste0(ref_gene, " (", gene, ")")
                   )
      }
    } else {
      cat("Gene", gene, "not found.\n")
    }
}

# heatmap of DTU genes ----
# use CPM
counts <- read.delim("../bambu/out/counts_transcript.txt")
rownames(counts) <- counts$TXNAME
counts[, seq(4, 18, 2)] <- counts[, seq(4, 18, 2)] + counts[, seq(3, 17, 2)]
counts <- counts[, -seq(3, 17, 2)]
colnames(counts) <- gsub("_pass", "", colnames(counts))
cpm <- cpm(counts[, c(-1, -2)], log=TRUE)
make_gen_heatmap <- function(gene, anno.row=NULL){
  anno <- gr %>% filter(gene_id == gene, type == "transcript") %>% as.data.frame
  if(nrow(anno) == 1){
    cat("Gene ", gene, "only have one transcript.\n")
  } else {
    dat <- cpm[anno$transcript_id, ]
    anno_col <- data.frame(
      siRNA = rep(c("NT", "RNPC3"), each = 4)
    )
    rownames(anno_col) <- paste0("barcode0", 1:8)
    if(is.null(anno.row)){
      pheatmap(dat,
               cluster_rows = FALSE,
               cluster_cols = FALSE,
               scale = "column",
               annotation_col = anno_col,
               main = gene)
    } else {
      anno_row <- as.data.frame(anno[, anno.row]) %>%
        dplyr::mutate_if(is.logical, as.character)
      colnames(anno_row) <- anno.row
      rownames(anno_row) <- anno$transcript_id
      pheatmap(dat,
               cluster_rows = FALSE,
               cluster_cols = FALSE,
               scale = "column",
               annotation_col = anno_col,
               annotation_row = anno_row,
               main = gene)
    }
  }
}
# Current issue: cannot find how the transcripts were ordered in Gviz, even after digging into the plot object.


# add DTE/DTU info and isoform category to annotation ----
tt <- read.csv("results/DE_topTags.csv")
ts <- read.csv("results/DTU_transcript_topSplice.csv")
gr$isDTE <- gr$transcript_id %in% (tt %>% dplyr::filter(FDR < 0.05) %>% dplyr::pull(qry_id))
gr$isDTU <- gr$transcript_id %in% (ts %>% dplyr::filter(FDR < 0.05) %>% dplyr::pull(qry_id))
gr$transcript_id_status <- gr$transcript_id
gr$transcript_id_status[gr$isDTE] <- paste(gr$transcript_id_status[gr$isDTE], "(DTE)")
gr$transcript_id_status[gr$isDTU] <- paste(gr$transcript_id_status[gr$isDTU], "(DTU)")

gr$class_code <- tmap$class_code[match(gr$transcript_id, tmap$qry_id)]
gr$transcript_id_status <- paste(gr$transcript_id_status, gr$class_code)

# make plot for DTE/DTU using gviz ----
# genes <- list()
for(i in c("DTEgenesMI", "DTEgenesNotMI", "DTUgenesMI", "DTUgenesNotMI")){
  cat("Working on", i, ".\n")
  # genes[[i]] <- readRDS(paste0(i, ".RDS"))
  pdf(paste0("plots/cov_", i, "2.pdf"))
  for(j in genes[[i]]){
    plot_cov_genes2(j)
  }
  dev.off()
}

# plot for non DTE/DTU minor intron genes
mi_nosig_genes <- unique(na.omit(tmap$qry_gene_id[match(mi_gene$Gene.Name, tmap$ref_gene_id)]))
mi_nosig_genes <- mi_nosig_genes[!(mi_nosig_genes %in% union(genes$DTEgenesMI, genes$DTUgenesMI))]
mi_nosig_genes <- mi_nosig_genes[sapply(mi_nosig_genes, function(x){
  !!sum(c("m", "n", "j") %in% filter(gr, gene_id == x)$class_code)
}, simplify = TRUE)]
pdf("plots/cov_notSigMI_mnj.pdf")
for(i in mi_nosig_genes){
  plot_cov_genes2(i, anno_col = "transcript_id_status")
}
dev.off()

# make plot for Stephen's gene lists ----
## extract tables from pdf ----
library(tabulizer)
minor_ir <- as.data.frame(extract_tables("../metadata/Short-read splicing analysis.pdf", pages = 1)[[1]])
colnames(minor_ir) <- minor_ir[1 ,]
minor_ir <- minor_ir[-1, ]
colnames(minor_ir)[6:10] <- colnames(minor_ir)[5:9]
colnames(minor_ir)[4] <- "Coord start"
colnames(minor_ir)[5] <- "Coord end"
minor_as <- as.data.frame(extract_tables("../metadata/Short-read splicing analysis.pdf", pages = 2)[[1]])
colnames(minor_as) <- minor_as[1, ]
minor_as <- minor_as[-1, ]
colnames(minor_as)[6:10] <- colnames(minor_as)[5:9]
colnames(minor_as)[4] <- "Coord start" 
colnames(minor_as)[5] <- "Coord end"
ir_finder <- as.data.frame(extract_tables("../metadata/Short-read splicing analysis.pdf", pages = 3)[[1]])
colnames(ir_finder) <- ir_finder[1, ]
ir_finder <- ir_finder[-1, ]
## make plot ----
### minor ir ----
genes <- unique(refmap$qry_gene_id[refmap$ref_gene_id %in% minor_ir$`Gene Name`])
pdf("plots/list/minor_ir_coverage.pdf")
for(i in genes){
  plot_cov_genes2(i)
}
dev.off()
pdf("plots/list/minor_ir_heatmap.pdf")
for(i in genes){
  make_gen_heatmap(i, c("class_code", "isDTU", "isDTE"))
}
dev.off()
### minor as ----
genes <- limma::strsplit2(minor_as$)

# test----
plot_cov_genes2("TMEM80_1", anno_col = "transcript_id_status")
plot_cov_genes2("TMEM80_1")
plot_cov_genes2("UBL5_1")
plot_cov_genes2("SPCS2_1")
plot_cov_genes2("KRTCAP2_1")
pdf("test.pdf")
for(i in c("UBL5_1", "SPCS2_1")){
  plot_cov_genes2(i)
}
dev.off()

pdf("test2.pdf", width = 10)
# par(mfrow=c(1, 2))
plot_cov_genes2("TMEM80_1", cov_data = "long")
make_gen_heatmap("TMEM80_1", c("isDTU", 'isDTE', "class_code"))
dev.off()
