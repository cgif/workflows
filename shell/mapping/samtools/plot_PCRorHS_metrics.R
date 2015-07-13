library("ggplot2")
library("scales")

metrics.file <- "#metricsSummary"
experiment <- "#experimentType"

metrics <- read.delim(metrics.file, as.is=T, row.names=1)
width <- length(rownames(metrics)) * 15

#plot uniq and aligned reads
reads.aligned.ext <- "readsAlligned"

names <- rep(rownames(metrics), 3)
status <- c(rep("PF_READS",length(rownames(metrics))),rep("PF_UNIQUE_READS",length(rownames(metrics))),rep("PF_UQ_READS_ALIGNED",length(rownames(metrics))))
reads <- c(metrics$PF_READS, metrics$PF_UNIQUE_READS, metrics$PF_UQ_READS_ALIGNED)
reads.aligned.df <- data.frame(names=factor(names,levels=rownames(metrics)), status=factor(status,levels=c("PF_READS","PF_UNIQUE_READS","PF_UQ_READS_ALIGNED")), reads=reads)

png(paste(metrics.file, ".", reads.aligned.ext, ".png", sep=""), width=width, height=500)

ggplot(reads.aligned.df,aes(names,reads,fill=status))+geom_bar(stat="identity",position=position_dodge())+scale_fill_hue(h.start=240,guide=guide_legend(title=NULL),labels=c("total","unique","unique and aligned"))+ggtitle("UNIQUE AND ALIGNED READS")+scale_x_discrete("")+scale_y_continuous("",expand = c(0, 0),labels = comma)+theme(plot.title=element_text(face="bold",vjust=1),axis.text.x=element_text(angle=90,vjust=0.5),panel.background=element_rect(fill="white"),panel.grid.major.x=element_line(colour="NA"),panel.grid.major.y=element_line(colour="light grey"))

dev.off()

#plot reads ON and OFF bait
bases.on.bait.ext <- "basesOnBait"

names <- rep(rownames(metrics), 3)

if (experiment == "TP") {

	status <- c(rep("ON_AMPLICON_BASES",length(rownames(metrics))),rep("NEAR_AMPLICON_BASES",length(rownames(metrics))),rep("OFF_AMPLICON_BASES",length(rownames(metrics))))
	bases <- c(metrics$ON_AMPLICON_BASES, metrics$NEAR_AMPLICON_BASES, metrics$OFF_AMPLICON_BASES)
	bases.on.bait.df <- data.frame(names=factor(names,levels=rownames(metrics)), status=factor(status,levels=c("ON_AMPLICON_BASES","NEAR_AMPLICON_BASES","OFF_AMPLICON_BASES")), bases=bases)
	title <- "BASES ON AMPLICON"
	labels <- c("on amplicon","near amplicon","off amplicon")
} 

if (experiment == "HS") {

	status <- c(rep("ON_BAIT_BASES",length(rownames(metrics))),rep("NEAR_BAIT_BASES",length(rownames(metrics))),rep("OFF_BAIT_BASES",length(rownames(metrics))))
	bases <- c(metrics$ON_BAIT_BASES, metrics$NEAR_BAIT_BASES, metrics$OFF_BAIT_BASES)
	bases.on.bait.df <- data.frame(names=factor(names,levels=rownames(metrics)), status=factor(status,levels=c("ON_BAIT_BASES","NEAR_BAIT_BASES","OFF_BAIT_BASES")), bases=bases)
	title <- "BASES ON BAIT"
	labels <- c("on bait","near bait","off bait")
}

png(paste(metrics.file, ".", bases.on.bait.ext, ".png", sep=""), width=width, height=500)

ggplot(bases.on.bait.df,aes(names,bases,fill=status))+geom_bar(stat="identity")+scale_fill_hue(h.start=240,guide=guide_legend(title=NULL),labels=labels)+ggtitle(title)+scale_x_discrete("")+scale_y_continuous("",expand = c(0, 0),labels = comma)+theme(plot.title=element_text(face="bold",vjust=1),axis.text.x=element_text(angle=90,vjust=0.5),panel.background=element_rect(fill="white"),panel.grid.major.x=element_line(colour="NA"),panel.grid.major.y=element_line(colour="light grey"))

dev.off()

#plot target coverage
mean.target.coverage.ext <- "meanTargetCoverage"

png(paste(metrics.file, ".", mean.target.coverage.ext, ".png", sep=""), width=width, height=500)

ggplot(metrics,aes(rownames(metrics),MEAN_TARGET_COVERAGE))+geom_bar(stat="identity",fill="#DD8888")+ggtitle("MEAN TARGET COVERAGE")+scale_x_discrete("")+scale_y_continuous("",expand = c(0, 10))+theme(plot.title=element_text(face="bold",vjust=1),axis.text.x=element_text(angle=90,vjust=0.5),panel.background=element_rect(fill="white"),panel.grid.major.x=element_line(colour="NA"),panel.grid.major.y=element_line(colour="light grey"))

dev.off()

#plot fold enrichment
fold.enrichment.ext <- "foldEnrichment"

png(paste(metrics.file, ".", fold.enrichment.ext, ".png", sep=""), width=width, height=500)

ggplot(metrics,aes(rownames(metrics),FOLD_ENRICHMENT))+geom_bar(stat="identity",fill="#DD8888")+ggtitle("FOLD ENRICHMENT")+scale_x_discrete("")+scale_y_continuous("",expand = c(0, 10))+theme(plot.title=element_text(face="bold",vjust=1),axis.text.x=element_text(angle=90,vjust=0.5),panel.background=element_rect(fill="white"),panel.grid.major.x=element_line(colour="NA"),panel.grid.major.y=element_line(colour="light grey"))

dev.off()

#plot cumulative coverage
cumulative.coverage.ext <- "cumulativeCoverage"

names <- rep(rownames(metrics), 4)
coverage <- c(rep("x2",length(rownames(metrics))),rep("x10",length(rownames(metrics))),rep("x20",length(rownames(metrics))),rep("x30",length(rownames(metrics))))
targets.covered <- c(metrics$PCT_TARGET_BASES_2X, metrics$PCT_TARGET_BASES_10X, metrics$PCT_TARGET_BASES_20X, metrics$PCT_TARGET_BASES_30X)
coverage.df <- data.frame(names=factor(names,levels=rownames(metrics)), coverage=factor(coverage,levels=c("x2","x10","x20","x30")), targets.covered=targets.covered)

png(paste(metrics.file, ".", cumulative.coverage.ext, ".png", sep=""), width=width, height=500)

ggplot(coverage.df,aes(coverage,targets.covered,colour=names,group=names))+geom_line()+guides(col = guide_legend(title=NULL,nrow = 25))+ggtitle("FRACTION OF TARGETS WITH X COVERAGE")+scale_x_discrete("COVERAGE")+scale_y_continuous("",expand = c(0, 0),lim=c(0,1))+theme(plot.title=element_text(face="bold",vjust=1),panel.background=element_rect(fill="white"),panel.grid.major.x=element_line(colour="NA"),panel.grid.major.y=element_line(colour="light grey"))

dev.off()


