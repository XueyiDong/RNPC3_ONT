# TO DO LIST ----
# 1. add option in function for adding alignment track
# 2. make short read heatmap

## DONE LIST ----
# 1. add short read coverage
# 2. inport Stephen's gene lists
# 3. make heatmaps
# 4. Select Stephen's highlighted genes

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
mi$gene_id <- refmap$qry_gene_id[match(mi$name, refmap$ref_gene_id)]

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
                  data = "score",
                  type = "h")
      })
      return(dTrack)
    }
  }
}

plot_cov_genes2 <- function(gene, cov_data = c("long", "short"), anno_col = "transcript_id", 
                            hightlight_range = mi, alignment_track = FALSE){
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
        # add alignment track
        if(alignment_track){
          # anno_range <- gr %>% filter(gene_id == gene, type == "gene")
          alTrack_long_NT <- AlignmentsTrack("../aligned_minimap2/merged_NT.bam",
                                             isPaired = FALSE,
                                             options(ucscChromosomeNames=FALSE),
                                             type = "pileup",
                                             name = "long NT")
          alTrack_long_RNPC3 <- AlignmentsTrack("../aligned_minimap2/merged_RNPC3.bam",
                                                isPaired = FALSE,
                                                options(ucscChromosomeNames=FALSE),
                                                type = "pileup",
                                                name = "long RNPC3")
          alTrack_short_NT <- AlignmentsTrack("../short_bam_merged/merged_NT.bam",
                                              isPaired = TRUE,
                                              options(ucscChromosomeNames=FALSE),
                                              type = "pileup",
                                              name = "short NT")
          alTrack_short_RNPC3 <- AlignmentsTrack("../short_bam_merged/merged_RNPC3.bam",
                                                 isPaired = TRUE,
                                                 options(ucscChromosomeNames=FALSE),
                                                 type = "pileup",
                                                 name = "short RNPC3")
          trackList<- Reduce(append, c(alTrack_long_NT, alTrack_long_RNPC3, alTrack_short_NT, alTrack_short_RNPC3, trackList))
        }
        # add highlight for minor intron region
        ref_gene <- tmap$ref_gene_id[match(gene, tmap$qry_gene_id)][1]
        ht <- HighlightTrack(trackList = trackList,
                             range = hightlight_range %>% filter(name == ref_gene))
        plotTracks(list(ht, axisTrack),
                   groupAnnotation = "group",
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
make_gen_heatmap <- function(gene, anno.row=NULL, scale = "none", 
                             cluster_rows = TRUE, cluster_cols = FALSE){
  anno <- gr %>% filter(gene_id == gene, type == "transcript") %>% as.data.frame
  if(nrow(anno) == 1){
    cat("Gene", gene, "only have one transcript.\n")
    dat <- cpm[anno$transcript_id, ]
    dat <- data.frame(
      sample = names(dat),
      logCPM = dat,
      group = rep(c("NT", "RNPC3"), each = 4)
    )
    p <- ggplot(dat, aes(x = sample, y = logCPM, fill = group)) +
      geom_bar(stat="identity")+
      ggtitle(gene)
    print(p)
  } else {
    dat <- cpm[anno$transcript_id, ]
    anno_col <- data.frame(
      siRNA = rep(c("NT", "RNPC3"), each = 4)
    )
    rownames(anno_col) <- paste0("barcode0", 1:8)
    if(is.null(anno.row)){
      pheatmap(dat,
               cluster_rows = cluster_rows,
               cluster_cols = cluster_cols,
               scale = scale,
               annotation_col = anno_col,
               main = gene)
    } else {
      anno_row <- as.data.frame(anno[, anno.row]) %>%
        dplyr::mutate_if(is.logical, as.character)
      colnames(anno_row) <- anno.row
      rownames(anno_row) <- anno$transcript_id
      pheatmap(dat,
               cluster_rows = cluster_rows,
               cluster_cols = cluster_cols,
               scale = scale,
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

# WARNING: this should be run after imporint Stephen's lists
gr$ref_gene_id <- refmap$ref_gene_id[match(gr$gene_id, refmap$qry_gene_id)]
gr$minor_ir <- gr$ref_gene_id %in% minor_ir$`Gene Name`
gr$minor_as <- gr$ref_gene_id %in% limma::strsplit2(minor_as$`AS event – Gene ID – Gene Name`, "-", fixed=TRUE)[,3]
gr$ir_finder <- gr$ref_gene_id %in%ir_finder$Gene

# make plot for Stephen's gene lists ----
library(tabulizer)

## minor ir ----
### extract lists ----
minor_ir <- as.data.frame(extract_tables("../metadata/Short-read splicing analysis.pdf", pages = 1)[[1]])
colnames(minor_ir) <- minor_ir[1 ,]
minor_ir <- minor_ir[-1, ]
colnames(minor_ir)[6:10] <- colnames(minor_ir)[5:9]
colnames(minor_ir)[4] <- "Coord start"
colnames(minor_ir)[5] <- "Coord end"

#### prepare minor intron region annotation ----
# # minor_ir_range <- GRanges(
# #   seqnames = minor_ir$Chr,
# #   ranges = IRanges(as.numeric(minor_ir$'Coord start'),
# #                    end = as.numeric(minor_ir$'Coord end'),
# #                    names = minor_ir$'Gene Name')
# # )
# # write_bed(minor_ir_range, "minor_ir_range.bed")
# # #intron coord converted using online tool: https://genome.ucsc.edu/cgi-bin/hgLiftOver
# minor_ir_range <- read_bed("hglft_genome_ceec_3de600.bed")
# minor_ir_range$gene_id <- refmap$qry_gene_id[match(minor_ir_range$name, refmap$ref_gene_id)]

### make plots ----
genes <- unique(refmap$qry_gene_id[refmap$ref_gene_id %in% minor_ir$`Gene Name`])
pdf("plots/list/minor_ir_coverage.pdf")
for(i in genes){
  plot_cov_genes2(i)
}
dev.off()
pdf("plots/list/minor_ir_heatmap.pdf")
for(i in genes){
  make_gen_heatmap(i, c("class_code", "isDTU", "isDTE", "minor_ir", "minor_as", "ir_finder"))
}
dev.off()

## minor as ----
### extract lists ----
minor_as <- as.data.frame(extract_tables("../metadata/Short-read splicing analysis.pdf", pages = 2)[[1]])
colnames(minor_as) <- minor_as[1, ]
minor_as <- minor_as[-1, ]
minor_as <- minor_as[, -9]
colnames(minor_as)[5:9] <- colnames(minor_as)[4:8]
colnames(minor_as)[3] <- "Coord start" 
colnames(minor_as)[4] <- "Coord end"
### make plots ----
genes <- unique(limma::strsplit2(minor_as[,1], "-", fixed=TRUE)[,3])
genes <- unique(refmap$qry_gene_id[refmap$ref_gene_id %in% genes])
pdf("plots/list/minor_as_coverage.pdf")
for(i in genes){
  plot_cov_genes2(i)
}
dev.off()
pdf("plots/list/minor_as_heatmap.pdf")
for(i in genes){
  make_gen_heatmap(i, c("class_code", "isDTU", "isDTE", "minor_ir", "minor_as", "ir_finder"))
}
dev.off()

## ir_finder ----
### extract lists ----
ir_finder <- as.data.frame(extract_tables("../metadata/Short-read splicing analysis.pdf", pages = 3)[[1]])
colnames(ir_finder) <- ir_finder[1, ]
ir_finder <- ir_finder[-1, ]
### prepare minor intron region annotation ----
ir_finder_range <- GRanges(
  seqnames = paste0("chr", ir_finder$Chr),
  ranges = IRanges(as.numeric(ir_finder$Start),
                   end = as.numeric(ir_finder$End),
                   names = ir_finder$Gene)
)
write_bed(ir_finder_range, "ir_finder_range.bed")
#intron coord converted using online tool: https://genome.ucsc.edu/cgi-bin/hgLiftOver
ir_finder_range <- read_bed("hglft_genome_2b537_3f5770.bed")
ir_finder_range$gene_id <- refmap$qry_gene_id[match(ir_finder_range$name, refmap$ref_gene_id)]
### make plots ----
genes <- unique(refmap$qry_gene_id[refmap$ref_gene_id %in% ir_finder$Gene])
pdf("plots/list/ir_finder_coverage2.pdf")
for(i in genes){
  plot_cov_genes2(i, hightlight_range = ir_finder_range)
}
dev.off()
pdf("plots/list/ir_finder_heatmap.pdf")
for(i in genes){
  make_gen_heatmap(i, c("class_code", "isDTU", "isDTE", "minor_ir", "minor_as", "ir_finder"))
}
dev.off()

## Gene examples ----
gene_examples <- xlsx::read.xlsx("../metadata/Gene examples and thoughts.xlsx", sheetIndex = 1)
genes <- refmap$qry_gene_id[match(na.omit(gene_examples$Gene), refmap$ref_gene_id)]
pdf("plots/list/gene_examples_coverage.pdf")
for(i in genes){
  plot_cov_genes2(i)
}
dev.off()
pdf("plots/list/gene_examples_heatmap.pdf")
for(i in genes){
  make_gen_heatmap(i, c("class_code", "isDTU", "isDTE", "minor_ir", "minor_as", "ir_finder"))
}
dev.off()

# make plot for DTE/DTU using gviz ----
genes <- list()
for(i in c("DTEgenesMI", "DTEgenesNotMI", "DTUgenesMI", "DTUgenesNotMI")){
  cat("Working on", i, ".\n")
  genes[[i]] <- readRDS(paste0(i, ".RDS"))
  pdf(paste0("plots/cov_", i, "2.pdf"))
  for(j in genes[[i]]){
    plot_cov_genes2(j)
  }
  dev.off()
  pdf(paste0("plots/heatmap_", i, "2.pdf"))
  for(j in genes[[i]]){
    make_gen_heatmap(j, c("class_code", "isDTU", "isDTE", "minor_ir", "minor_as", "ir_finder"))
  }
  dev.off()
}

# plot for non DTE/DTU minor intron genes
mi_nosig_genes <- unique(na.omit(refmap$qry_gene_id[match(mi_gene$Gene.Name, refmap$ref_gene_id)]))
mi_nosig_genes <- mi_nosig_genes[!(mi_nosig_genes %in% union(genes$DTEgenesMI, genes$DTUgenesMI))]
mi_nosig_genes <- mi_nosig_genes[sapply(mi_nosig_genes, function(x){
  !!sum(c("m", "n", "j") %in% filter(gr, gene_id == x)$class_code)
}, simplify = TRUE)]
pdf("plots/cov_notSigMI_mnj.pdf")
for(i in mi_nosig_genes){
  plot_cov_genes2(i, anno_col = "transcript_id_status")
}
dev.off()



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

pdf("plots/test.pdf", width = 10)
# par(mfrow=c(1, 2))
plot_cov_genes2("TMEM80_1", cov_data = NULL, alignment_track = TRUE)
# plot_cov_genes2("TMEM80_1", cov_data = "long")
# make_gen_heatmap("TMEM80_1", c("isDTU", 'isDTE', "class_code"))
dev.off()

test_rg <- gr %>% filter(gene_id=="TMEM80_1", type=="gene")
test_alTrack_long <- AlignmentsTrack(design$Bam[1], isPaired = FALSE,
                                     chromosome = as.character(seqnames(test_rg)),
                                     start = start(test_rg),
                                     end = end(test_rg),
                                # range = gr %>% filter(gene_id=="TMEM80_1", type=="gene"),
                                options(ucscChromosomeNames=FALSE))
pdf("plots/testAlTr.pdf", width = 10)
plotTracks(test_alTrack_long,  from=start(test_rg),to=end(test_rg))
dev.off()