library("ggplot2")
library("scales")

metrics.file <- "#metricsSummary"

metrics <- read.delim(metrics.file, as.is=T, row.names=1, header=FALSE, skip=1)
titles <- read.delim(metrics.file, as.is=T, row.names=1, header=FALSE, nrows=1)
colnames(metrics) <- as.vector(as.matrix(titles))

width <- length(rownames(metrics)) * 15
if (width < 500) width <- 500

#plot bases genomic distribution
bases.alignment.ext <- "basesAlignment"

names <- rep(rownames(metrics), 5)
status <- c(rep("MRNA",length(rownames(metrics))),rep("INTRONIC",length(rownames(metrics))),rep("INTERGENIC",length(rownames(metrics))),rep("RIBOSOMAL",length(rownames(metrics))),rep("UNALIGNED",length(rownames(metrics))))
mrna.bases <- metrics$CODING_BASES + metrics$UTR_BASES
intronic.bases <- metrics$INTRONIC_BASES
intergenic.bases <- metrics$INTERGENIC_BASES
ribosomal.bases <- metrics$RIBOSOMAL_BASES
unaligned.bases <- metrics$PF_BASES - metrics$PF_ALIGNED_BASES
bases <- c(mrna.bases,intronic.bases,intergenic.bases,ribosomal.bases,unaligned.bases)
millions <- rep(1000000,length(bases))
bases <- bases/millions
bases.alignment.df <- data.frame(names=factor(names,levels=rownames(metrics)), status=factor(status,levels=c("MRNA","INTRONIC","INTERGENIC","RIBOSOMAL","UNALIGNED")), bases=bases)

png(paste(metrics.file, ".", bases.alignment.ext, ".png", sep=""), width=width+120, height=500)

ggplot(bases.alignment.df,aes(names,bases,fill=status))+geom_bar(stat="identity")+scale_fill_hue(guide=guide_legend(title=NULL),labels=c("mRNA","intron","intergenic region","rRNA","unaligned"))+ggtitle("GENOMIC DISTRIBUTION OF ALIGNED BASES")+scale_x_discrete("SAMPLE ID")+scale_y_continuous("BASES, M",expand = c(0, 0),labels = comma)+theme(plot.title=element_text(face="bold",vjust=2),axis.title.x=element_text(vjust=-0.5),axis.text.x=element_text(angle=90,vjust=0.5),axis.title.y=element_text(vjust=0.3),panel.background=element_rect(fill="white"),panel.grid.major.x=element_line(colour="NA"),panel.grid.major.y=element_line(colour="light grey"))

dev.off()

#plot median 5' to 3' bias
median.bias.ext <- "median5primeTo3primeBias"

png(paste(metrics.file, ".", median.bias.ext, ".png", sep=""), width=width, height=500)

ggplot(metrics,aes(rownames(metrics),MEDIAN_5PRIME_TO_3PRIME_BIAS))+geom_bar(stat="identity",fill="#DD8888")+ggtitle("MEDIAN RATIO OF COVEARGE AT 5' END TO 3' END\n for 1000 most hightly expressed transcripts")+scale_x_discrete("SAMPLE ID")+scale_y_continuous("RATIO",expand = c(0, 0))+theme(plot.title=element_text(face="bold",vjust=2),axis.title.x=element_text(vjust=-0.5),axis.text.x=element_text(angle=90,vjust=0.5),axis.title.y=element_text(vjust=0.3),panel.background=element_rect(fill="white"),panel.grid.major.x=element_line(colour="NA"),panel.grid.major.y=element_line(colour="light grey"))

dev.off()

