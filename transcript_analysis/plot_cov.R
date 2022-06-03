library(superintronic)
suppressPackageStartupMessages(library(plyranges))
library(ggplot2)

# prepare annotation----
# gff <- "../bambu/out/extended_annotations.genes.gtf"
# gr <- read_gff(gff, genome_info = "CHM13v2")
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

# plot for each feature ----
plot_cov_genes <- function(genes){
  plot_list = list(genes)
  for (i in genes){
    features <- parts %>% filter(gene_id == i)
    cvg_over_features <- cvg %>% 
      select(-Bam) %>% 
      join_parts(features)
    gene_track <- view_segments(unnest_parts(features), 
                                colour = feature_type)
    p <- cvg_over_features %>% 
      mutate(strand = feature_strand) %>% 
      view_coverage(score = score, 
                    colour = feature_type, 
                    facets = vars(siRNA)) + 
      scale_color_brewer(palette = "Dark2") +
      guides(colour = FALSE) +
      labs(title = i)
    plot_list[[i]] = p / gene_track + patchwork::plot_layout(heights = c(3, 1))
  }
  return(plot_list)
}

# make plot ---
for(i in c("DTEgenesMI", "DTEgenesNotMI", "DTUgenesMI", "DTUgenesNotMI")){
  genes <- readRDS(paste0(i, ".RDS"))
  plot_list <- plot_cov_genes(genes)
  pdf(paste0("cov_", i, ".pdf"))
  for(j in genes){
    print(plot_list[[j]])
  }
  dev.off()
}

