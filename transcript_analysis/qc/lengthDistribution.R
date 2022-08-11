library(ggplot2)
library(MetBrewer)

qcdata <- readRDS("./rds/summaryInfo.RDS")


# dim(qcdata)
# table(qcdata$Barcode)
# summary(qcdata$Read_length)
# summary(qcdata$Read_length[qcdata$Pass_filtering == "TRUE"])
# summary(qcdata$Read_length[qcdata$Pass_filtering == "FALSE"])
# quantile(qcdata$Read_length, 0.999)
# summary(qcdata$Qscore)
# quantile(qcdata$Qscore, 0.999)

lenLim <- quantile(qcdata$Read_length, 0.9999)
lenLim
pdf("plots/LengthDistSample.pdf", height = 5, width = 8)
ggplot(qcdata, aes(x=Read_length, colour=Barcode, linetype=Barcode)) +
  geom_density() +
  theme_bw() +
  theme(text = element_text(size = 20)) +
  scale_x_continuous(trans="log10", limits = c(NA, lenLim)) +
  labs(x= "ONT read length", colour="sample", linetype="sample") +
  scale_colour_manual(values = c(met.brewer("Tiepolo", 8), "grey70")) +
  scale_linetype_manual(values = c(rep("solid", 8), "dashed"))
# facet by whether passed filtering
ggplot(qcdata, aes(x=Read_length, colour=Barcode, linetype=Barcode)) +
  geom_density() +
  facet_grid(rows=vars(Pass_filtering)) +
  theme_bw() +
  theme(text = element_text(size = 20)) +
  scale_x_continuous(trans="log10", limits = c(NA, lenLim)) +
  labs(x= "ONT read length", colour="sample", linetype="sample") +
  scale_colour_manual(values = c(met.brewer("Tiepolo", 8), "grey70")) +
  scale_linetype_manual(values = c(rep("solid", 8), "dashed"))
dev.off()

# undemultiplexed reads are filtered out
pdf("plots/LengthDist.pdf", height = 5, width = 8)
ggplot(qcdata[qcdata$Barcode != "other", ], aes(x=Read_length)) +
  geom_density(colour="black") +
  geom_histogram(aes(y=..density..), fill="#438DAC", colour="#438DAC", alpha = .3, bins=50) +
  theme_bw() +
  theme(text = element_text(size = 20)) +
  scale_x_continuous(trans="log10", limits = c(NA, lenLim)) +
  labs(x= "ONT read length")
dev.off()

pdf("plots/LengthDistFilter.pdf", height = 5, width = 8)
ggplot(qcdata[qcdata$Barcode != "other", ], aes(x=Read_length)) +
  geom_density(colour="black") +
  geom_histogram(aes(y=..density..), fill="#438DAC", colour="#438DAC", alpha = .3, bins=50) +
  facet_grid(rows=vars(Pass_filtering)) +
  theme_bw() +
  theme(text = element_text(size = 20)) +
  scale_x_continuous(trans="log10", limits = c(NA, lenLim)) +
  labs(x= "ONT read length")
dev.off()
