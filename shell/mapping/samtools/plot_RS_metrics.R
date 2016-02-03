library("ggplot2")
library("scales")
library(RColorBrewer)

metrics.file <- "#metricsSummary"
counts.file <- "#countsFile"
ref.file <- "#refFile"

#metrics.file <- "/project/tgu/results/standfield_claudication/mergetag/2015-07-27/multisample/standfield_claudication.2015-07-27.RnaSeqMetrics"
#counts.file <- "/project/tgu/results/standfield_claudication/mergetag/2015-07-27/multisample/standfield_claudication.2015-07-27.readCounts"

metrics <- read.delim(metrics.file, as.is=T, row.names=1, header=FALSE, skip=1)
titles <- read.delim(metrics.file, as.is=T, row.names=1, header=FALSE, nrows=1)
colnames(metrics) <- as.vector(as.matrix(titles))

width <- length(rownames(metrics)) * 15
if (width < 500) width <- 500

#plot bases genomic distribution

names <- rep(rownames(metrics), 5)
status <- c(rep("RNA",length(rownames(metrics))),rep("INTRONIC",length(rownames(metrics))),rep("INTERGENIC",length(rownames(metrics))),rep("RIBOSOMAL",length(rownames(metrics))),rep("UNALIGNED",length(rownames(metrics))))
mrna.bases <- metrics$CODING_BASES + metrics$UTR_BASES
intronic.bases <- metrics$INTRONIC_BASES
intergenic.bases <- metrics$INTERGENIC_BASES
ribosomal.bases <- metrics$RIBOSOMAL_BASES
unaligned.bases <- metrics$PF_BASES - metrics$PF_ALIGNED_BASES
bases <- c(mrna.bases,intronic.bases,intergenic.bases,ribosomal.bases,unaligned.bases)
millions <- rep(1000000,length(bases))
bases <- bases/millions
bases.alignment.df <- data.frame(names=factor(names,levels=rownames(metrics)), status=factor(status,levels=c("RNA","INTRONIC","INTERGENIC","RIBOSOMAL","UNALIGNED")), bases=bases)
colors.manual <- c("cornflowerblue","orchid3","springgreen3","yellow3","salmon")

png(paste(metrics.file, ".basesAlignment.png", sep=""), width=width+120, height=500)

ggplot(bases.alignment.df,aes(names,bases,fill=status))+geom_bar(stat="identity")+scale_fill_hue(guide=guide_legend(title=NULL),labels=c("RNA","intron","intergenic region","rRNA","unaligned"))+ggtitle("GENOMIC DISTRIBUTION OF ALIGNED BASES")+scale_x_discrete("SAMPLE ID")+scale_y_continuous("BASES [M]",expand = c(0, 0),labels = comma)+theme(plot.title=element_text(face="bold",vjust=2),axis.title.x=element_text(vjust=-0.5),axis.text.x=element_text(angle=90,vjust=0.5),axis.title.y=element_text(vjust=0.3),panel.background=element_rect(fill="white"),panel.grid.major.x=element_line(colour="NA"),panel.grid.major.y=element_line(colour="light grey"))+scale_fill_manual(values=colors.manual)

dev.off()

#plot median 5' to 3' bias

png(paste(metrics.file, ".median5primeTo3primeBias.png", sep=""), width=width, height=500)

ggplot(metrics,aes(rownames(metrics),MEDIAN_5PRIME_TO_3PRIME_BIAS))+geom_bar(stat="identity",fill="#DD8888")+ggtitle("MEDIAN RATIO OF COVEARGE AT 5' END TO 3' END\n for 1000 most hightly expressed transcripts")+scale_x_discrete("SAMPLE ID")+scale_y_continuous("RATIO",expand = c(0, 0),limits=c(0,1))+theme(plot.title=element_text(face="bold",vjust=2),axis.title.x=element_text(vjust=-0.5),axis.text.x=element_text(angle=90,vjust=0.5),axis.title.y=element_text(vjust=0.3),panel.background=element_rect(fill="white"),panel.grid.major.x=element_line(colour="NA"),panel.grid.major.y=element_line(colour="light grey"))

dev.off()

#plot read counts per chromosome

counts <- read.delim(counts.file, as.is=T, header=FALSE, skip=1)
titles <- read.delim(counts.file, as.is=T, header=FALSE, nrows=1)
colnames(counts) <- as.vector(as.matrix(titles))

sample.name <- counts$sample
chrom.name <- counts$chrom 
chrom.counts <- counts$read_counts
thousands <- rep(1000,length(chrom.counts))
chrom.counts <- chrom.counts/thousands
counts.df <- data.frame(sample.name=factor(sample.name,levels=unique(sample.name)), chrom.name=factor(chrom.name,levels=unique(chrom.name)), chrom.counts=chrom.counts)

colors.count <- length(unique(chrom.name))
getPalette = colorRampPalette(brewer.pal(9, "Paired"))

png(paste(counts.file, ".png", sep=""), width=width+75, height=500)

ggplot(counts.df,aes(sample.name,chrom.counts,fill=chrom.name))+geom_bar(stat="identity")+scale_fill_manual(guide=guide_legend(title=NULL),labels=unique(chrom.name),values=getPalette(colors.count))+ggtitle("CHROMOSOMAL DISTRIBUTION OF ALIGNED READS")+scale_x_discrete("SAMPLE ID")+scale_y_continuous("READS [K]",expand = c(0, 0),labels = comma)+theme(plot.title=element_text(face="bold",vjust=2),axis.title.x=element_text(vjust=-0.5),axis.text.x=element_text(angle=90,vjust=0.5),axis.title.y=element_text(vjust=0.3),panel.background=element_rect(fill="white"),panel.grid.major.x=element_line(colour="NA"),panel.grid.major.y=element_line(colour="light grey"))

dev.off()

#plot normalized read counts per chromosome

if (ref.file != "/project/tgu/resources/reference/hsapiens/GRCh37/fasta/GRCh37.fa") { stop() }

counts <- read.delim(counts.file, as.is=T, header=FALSE, skip=1)
titles <- read.delim(counts.file, as.is=T, header=FALSE, nrows=1)
colnames(counts) <- as.vector(as.matrix(titles))

lengths.file <- "/project/tgu/resources/annotations/Homo_sapiens.GRCh37.75.transcript_per_chrom_length"
transcript.length <- read.delim(lengths.file, row.names=1, as.is=T, header=FALSE, skip=1)

sample.name <- counts$sample
chrom.name <- counts$chrom 
chrom.counts <- counts$read_counts
chrom.counts <- chrom.counts/transcript.length[chrom.name,]
counts.df <- data.frame(sample.name=factor(sample.name,levels=unique(sample.name)), chrom.name=factor(chrom.name,levels=unique(chrom.name)), chrom.counts=chrom.counts)

png(paste(counts.file, ".normalized.png", sep=""), width=width+75, height=500)

ggplot(counts.df,aes(sample.name,chrom.counts,fill=chrom.name))+geom_bar(stat="identity")+scale_fill_manual(guide=guide_legend(title=NULL),labels=unique(chrom.name),values=getPalette(colors.count))+ggtitle("CHROMOSOMAL DISTRIBUTION OF ALIGNED READS\nnormalized by the total length of transcripts annotated on the chromosome")+scale_x_discrete("SAMPLE ID")+scale_y_continuous("READS [N] / TRANSCRIPT LENGTH [BP]",expand = c(0, 0),labels = comma)+theme(plot.title=element_text(face="bold",vjust=2),axis.title.x=element_text(vjust=-0.5),axis.text.x=element_text(angle=90,vjust=0.5),axis.title.y=element_text(vjust=0.3),panel.background=element_rect(fill="white"),panel.grid.major.x=element_line(colour="NA"),panel.grid.major.y=element_line(colour="light grey"))

dev.off()

counts.df <- counts.df[counts.df$chrom.name != "MT",]
counts.df <- counts.df[counts.df$chrom.name != "GL",]

png(paste(counts.file, ".normalizedNoMTGL.png", sep=""), width=width+75, height=500)

ggplot(counts.df,aes(sample.name,chrom.counts,fill=chrom.name))+geom_bar(stat="identity")+scale_fill_manual(guide=guide_legend(title=NULL),labels=unique(chrom.name),values=getPalette(colors.count))+ggtitle("CHROMOSOMAL DISTRIBUTION OF ALIGNED READS\nnormalized by the total length of transcripts annotated on the chromosome")+scale_x_discrete("SAMPLE ID")+scale_y_continuous("READS [N] / TRANSCRIPT LENGTH [BP]",expand = c(0, 0),labels = comma)+theme(plot.title=element_text(face="bold",vjust=2),axis.title.x=element_text(vjust=-0.5),axis.text.x=element_text(angle=90,vjust=0.5),axis.title.y=element_text(vjust=0.3),panel.background=element_rect(fill="white"),panel.grid.major.x=element_line(colour="NA"),panel.grid.major.y=element_line(colour="light grey"))

dev.off()



