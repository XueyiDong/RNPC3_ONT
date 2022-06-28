library(superintronic)
suppressPackageStartupMessages(library(plyranges))
library(ggplot2)

# prepare annotation----
gff <- "../bambu/out/extended_annotations.genes.edited.gtf"
gr <- read_gff(gff, genome_info = "CHM13v2")
# parts <- collect_parts(gr)
# saveRDS(parts, "complete_annotation.rds")
parts <- readRDS("complete_annotation.rds")

# design ----
bam_dir = "../aligned_minimap2"
design <- data.frame(
  Sample = paste0(rep(paste(rep(c("NT", "targeted"), c(4, 4)), rep(1:4, 2), sep = "_"), rep(2, 8)), rep(c("_fail", "_pass"), 8)),
  Bam = file.path(bam_dir, list.files(bam_dir, ".bam$")),
  Replicate = rep(rep(1:4, 2), rep(2, 8)),
  siRNA = rep(c("NT", "RNPC3"), c(8,8))
)
# compute coverage ----
# cvg <- compute_coverage_long(design, source="bam")
# saveRDS(cvg, "complete_coverage.rds")
cvg <- readRDS("complete_coverage.rds")


# transform intron coord into new version----
# write minor intron information into bed format
tmap <- read.delim("../bambu/out/bambu_comp.extended_annotations.gtf.tmap")
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

add_mi <- function(gene_name, gene_track){
  out <- gene_track
  mi_range <- mi %>% filter(name == gene_name)
  if(length(mi_range) > 0){
    # cat("Adding", length(mi_range),  "minor introns.\n")
    for(j in 1:length(mi_range)){
      out <- out +
        geom_segment(aes(x = start(mi_range)[j], xend=end(mi_range)[j], y=0.5, yend = 0.5),
                     size = 10, alpha = 0.5)
    }
  }
  return(out)
}

# plot for each feature ----
plot_cov_genes <- function(gene){
  if (gene %in% parts$gene_id){
    # cat("Making plot for gene", gene, ".\n")
    features <- parts %>% filter(gene_id == gene)
    cvg_over_features <- cvg %>% 
      select(-Bam) %>% 
      join_parts(features)
    if(length(unique(BiocGenerics::strand(unnest_parts(features)))) > 1) {
      cat("Gene", gene, "strand not unique.\n")
    } else{
      gene_track <- view_segments(unnest_parts(features),
                                  colour = feature_type)
      # mark retained introm for mi gene
      gene_name <- tmap$ref_gene_id[match(gene, tmap$qry_gene_id)]
      gene_track <- add_mi(gene_name, gene_track)
      p <- cvg_over_features %>% 
        mutate(strand = feature_strand) %>% 
        view_coverage(score = score, 
                      colour = feature_type, 
                      facets = vars(siRNA)) + 
        scale_color_brewer(palette = "Dark2") +
        guides(colour = FALSE) +
        labs(title = gene)
      p / gene_track + patchwork::plot_layout(heights = c(3, 1))
    }
  } else {
    cat("Gene", gene, "not found.\n")
  }
}

# make plot ---
genes <- list()
for(i in c("DTEgenesMI", "DTEgenesNotMI", "DTUgenesMI", "DTUgenesNotMI")){
  cat("Working on", i, ".\n")
  genes[[i]] <- readRDS(paste0(i, ".RDS"))
  pdf(paste0("plots/cov_", i, ".pdf"))
  for(j in genes[[i]]){
    print(plot_cov_genes(j))
  }
  dev.off()
}

pdf("plots/cov_allMI.pdf")
for(i in unique(na.omit(tmap$qry_gene_id[match(mi_gene$Gene.Name, tmap$ref_gene_id)]))){
  print(plot_cov_genes(i))
}
dev.off()

# Visualize coverage using Gviz ----
library(Gviz)
plot_cov_genes2 <- function(gene){
    if (gene %in% parts$gene_id){
      # cat("Making plot for gene", gene, ".\n")
      features <- parts %>% filter(gene_id == gene)
      cvg_over_features <- cvg %>% 
        select(-Bam) %>% 
        join_parts(features)
      if(length(unique(BiocGenerics::strand(unnest_parts(features)))) > 1) {
        cat("Gene", gene, "strand not unique.\n")
      } else{
        # axis
        axisTrack <- GenomeAxisTrack()
        # transcripts annotation
        anno <- gr %>% filter(gene_id == gene, type=="exon")
        aTrack <- AnnotationTrack(anno, 
                                  group = anno$transcript_id_status,
                                  options(ucscChromosomeNames=FALSE),
                                  name = "Isoforms")
        # coverage histogram
        # calculate mean coverage score by group
        data <- cvg_over_features %>% plyranges::group_by(siRNA) %>% disjoin_ranges_directed(score = mean(score))
        dTrack <- lapply(unique(data$siRNA), function(x){
          DataTrack(range = data %>% filter(siRNA == x),
                    options(ucscChromosomeNames=FALSE), 
                    name = x,
                    data = "score")
        })
        # add highlight for minor intron region
        ref_gene <- tmap$ref_gene_id[match(gene, tmap$qry_gene_id)][1]
        ht <- HighlightTrack(trackList = append(dTrack, aTrack),
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

# make plot using gviz ---
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
  plot_cov_genes2(i)
}
dev.off()

# test----
plot_cov_genes2("UBL5_1")
plot_cov_genes2("MYCBP_1")
plot_cov_genes2("KRTCAP2_1")
pdf("test.pdf")
for(i in c("UBL5_1", "SPCS2_1")){
  plot_cov_genes2(i)
}
dev.off()