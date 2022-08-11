library(ggplot2)

qcdata <- readRDS("./rds/summaryInfo.RDS")
qcdata.filt <- qcdata[qcdata$Barcode != "other", ]
maxLength = max(qcdata.filt$Read_length)
maxLength
qcdata.filt$LengthGroup <- Hmisc::cut2(qcdata.filt$Read_length, cuts = c(0, 500, 1000,
                                                     2000, maxLength))
qcdata.filt$LengthGroup <- gsub(" ", "", qcdata.filt$LengthGroup)
qcdata.filt$LengthGroup <- gsub(",", ", ", qcdata.filt$LengthGroup)
qcdata.filt$LengthGroup <- factor(qcdata.filt$LengthGroup, levels = c(
  "[0, 500)", "[500, 1000)", "[1000, 2000)", paste0("[2000, ", maxLength, "]")
))

qLim <- quantile(qcdata$Qscore, 0.9999)
qLim

# stat_box_data <- function(y, upper_limit = max(qcdata$Qscore) * 1.15) {
stat_box_data <- function(y, upper_limit = qLim * 1.15) {
  return( 
    data.frame(
      y = 0.95 * upper_limit,
      label = length(y)
    )
  )
}

pdf("plots/lengthQualityViolin.pdf", height = 5, width = 8)
ggplot(qcdata.filt, aes(x=LengthGroup, y=Qscore, fill=LengthGroup, colour = LengthGroup)) +
  geom_violin(alpha = 0.4) +
  geom_boxplot(width = 0.2, outlier.colour = NA, alpha = 0) +
  theme_bw() +
  ylim(c(NA, qLim * 1.15)) +
  labs(x = "Read length")+
  stat_summary(
    fun.data = stat_box_data, 
    geom = "text", 
    hjust = 0.5,
    vjust = 0.9
  )  +
  theme(text = element_text(size=20), legend.position = "none",
        axis.text.x = element_text(angle = 30, hjust = 1))
dev.off()

# pdf("plots/sampleQualityBox.pdf", height = 4, width = 8)
# ggplot(qcdata, aes(x=Barcode, y=Qscore, fill=Barcode, colour = Barcode)) +
#   geom_boxplot(alpha = 0.4) +
#   theme_bw()  +
#   labs(x = "Read length")+
#   stat_summary(
#     fun.data = stat_box_data, 
#     geom = "text", 
#     hjust = 0.5,
#     vjust = 0.9
#   ) +
#   theme(text = element_text(size=20), legend.position = "none", axis.text.x = element_text(angle = 90))
# dev.off()
