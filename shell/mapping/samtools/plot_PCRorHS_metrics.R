library("ggplot2")
library("scales")

metrics.file <- "#metricsSummary"
experiment <- "#experimentType"

metrics <- read.delim(metrics.file, as.is=T, row.names=1)
width <- length(rownames(metrics)) * 15

#plot uniq and aligned reads
reads.aligned.ext <- "readsAlligned"

names <- rep(rownames(metrics), 3)
status <- c(rep("UQ_ALIGNED",length(rownames(metrics))),rep("AMB_ALIGNED",length(rownames(metrics))),rep("DUPLICATES",length(rownames(metrics))))
duplicates <- metrics$PF_READS - metrics$PF_UNIQUE_READS
amb.aligned <- metrics$PF_UNIQUE_READS - metrics$PF_UQ_READS_ALIGNED
uq.aligned <- metrics$PF_UQ_READS_ALIGNED
reads <- c(uq.aligned, amb.aligned, duplicates)
thousands <- rep(1000,length(reads))
reads <- reads/thousands
reads.aligned.df <- data.frame(names=factor(names,levels=rownames(metrics)), status=factor(status,levels=c("UQ_ALIGNED","AMB_ALIGNED","DUPLICATES")), reads=reads)

if (experiment == "TP") {

labels=c("uniquely mapped reads", "non-uniquely mapped reads", "unmapped reads")
width.add=170

} 

if (experiment == "HS") {

labels=c("non-duplicate, uniquely mapped reads", "non-duplicate, non-uniquely mapped reads", "unmapped + mapped duplicate reads")
width.add=230

}

png(paste(metrics.file, ".", reads.aligned.ext, ".png", sep=""), width=width+width.add, height=500)

ggplot(reads.aligned.df,aes(names,reads,fill=status))+geom_bar(stat="identity")+scale_fill_hue(h.start=240,guide=guide_legend(title=NULL),labels=labels)+ggtitle("UNIQUE AND ALIGNED READS")+scale_x_discrete("SAMPLE ID")+scale_y_continuous("READS [K]",expand = c(0, 0),labels = comma)+theme(plot.title=element_text(face="bold",vjust=2),axis.title.x=element_text(vjust=-0.5),axis.text.x=element_text(angle=90,vjust=0.5),axis.title.y=element_text(vjust=0.5),panel.background=element_rect(fill="white"),panel.grid.major.x=element_line(colour="NA"),panel.grid.major.y=element_line(colour="light grey"))

dev.off()

#plot reads ON and OFF bait
bases.on.bait.ext <- "basesOnBait"

names <- rep(rownames(metrics), 3)

if (experiment == "TP") {

	status <- c(rep("ON_AMPLICON_BASES",length(rownames(metrics))),rep("NEAR_AMPLICON_BASES",length(rownames(metrics))),rep("OFF_AMPLICON_BASES",length(rownames(metrics))))
	bases <- c(metrics$ON_AMPLICON_BASES, metrics$NEAR_AMPLICON_BASES, metrics$OFF_AMPLICON_BASES)
	millions <- rep(1000000,length(bases))
	bases <- bases/millions
	bases.on.bait.df <- data.frame(names=factor(names,levels=rownames(metrics)), status=factor(status,levels=c("ON_AMPLICON_BASES","NEAR_AMPLICON_BASES","OFF_AMPLICON_BASES")), bases=bases)
	title <- "BASES ON AMPLICON"
	labels <- c("on amplicon","near amplicon","off amplicon")
} 

if (experiment == "HS") {

	status <- c(rep("ON_BAIT_BASES",length(rownames(metrics))),rep("NEAR_BAIT_BASES",length(rownames(metrics))),rep("OFF_BAIT_BASES",length(rownames(metrics))))
	bases <- c(metrics$ON_BAIT_BASES, metrics$NEAR_BAIT_BASES, metrics$OFF_BAIT_BASES)
	millions <- rep(1000000,length(bases))
	bases <- bases/millions
	bases.on.bait.df <- data.frame(names=factor(names,levels=rownames(metrics)), status=factor(status,levels=c("ON_BAIT_BASES","NEAR_BAIT_BASES","OFF_BAIT_BASES")), bases=bases)
	title <- "BASES ON BAIT"
	labels <- c("on bait","near bait","off bait")
}

png(paste(metrics.file, ".", bases.on.bait.ext, ".png", sep=""), width=width+100, height=500)

ggplot(bases.on.bait.df,aes(names,bases,fill=status))+geom_bar(stat="identity")+scale_fill_hue(h.start=240,guide=guide_legend(title=NULL),labels=labels)+ggtitle(title)+scale_x_discrete("SAMPLE ID")+scale_y_continuous("BASES [M]",expand = c(0, 0),labels = comma)+theme(plot.title=element_text(face="bold",vjust=2),axis.title.x=element_text(vjust=-0.5),axis.text.x=element_text(angle=90,vjust=0.5),axis.title.y=element_text(vjust=0.5),panel.background=element_rect(fill="white"),panel.grid.major.x=element_line(colour="NA"),panel.grid.major.y=element_line(colour="light grey"))

dev.off()

#plot target coverage
mean.target.coverage.ext <- "meanTargetCoverage"

png(paste(metrics.file, ".", mean.target.coverage.ext, ".png", sep=""), width=width, height=500)

ggplot(metrics,aes(rownames(metrics),MEAN_TARGET_COVERAGE))+geom_bar(stat="identity",fill="#DD8888")+ggtitle("MEAN COVERAGE FOR TARGETS with at least x2 coverage at one base")+scale_x_discrete("SAMPLE ID")+scale_y_continuous("COVERAGE [X FOLDS]",expand = c(0, 10))+theme(plot.title=element_text(face="bold",vjust=2),axis.title.x=element_text(vjust=-0.5),axis.text.x=element_text(angle=90,vjust=0.5),axis.title.y=element_text(vjust=0.5),panel.background=element_rect(fill="white"),panel.grid.major.x=element_line(colour="NA"),panel.grid.major.y=element_line(colour="light grey"))

dev.off()

#plot fold enrichment
fold.enrichment.ext <- "foldEnrichment"

png(paste(metrics.file, ".", fold.enrichment.ext, ".png", sep=""), width=width, height=500)

ggplot(metrics,aes(rownames(metrics),FOLD_ENRICHMENT))+geom_bar(stat="identity",fill="#DD8888")+ggtitle("FOLD ENRICHMENT")+scale_x_discrete("SAMPLE ID")+scale_y_continuous("ENRICHMENT [X FOLDS]",expand = c(0, 10))+theme(plot.title=element_text(face="bold",vjust=2),axis.title.x=element_text(vjust=-0.5),axis.text.x=element_text(angle=90,vjust=0.5),axis.title.y=element_text(vjust=0.5),panel.background=element_rect(fill="white"),panel.grid.major.x=element_line(colour="NA"),panel.grid.major.y=element_line(colour="light grey"))

dev.off()

#plot cumulative coverage
cumulative.coverage.ext <- "cumulativeCoverage"

names <- rep(rownames(metrics), 4)
coverage <- c(rep(2,length(rownames(metrics))),rep(10,length(rownames(metrics))),rep(20,length(rownames(metrics))),rep(30,length(rownames(metrics))))
targets.covered <- c(metrics$PCT_TARGET_BASES_2X, metrics$PCT_TARGET_BASES_10X, metrics$PCT_TARGET_BASES_20X, metrics$PCT_TARGET_BASES_30X)
percent <- rep(100,length(targets.covered))
targets.covered <- targets.covered*percent
coverage.df <- data.frame(names=factor(names,levels=rownames(metrics)), coverage=coverage, targets.covered=targets.covered)

png(paste(metrics.file, ".", cumulative.coverage.ext, ".png", sep=""), width=width, height=500)

ggplot(coverage.df,aes(coverage,targets.covered,colour=names,group=names))+geom_point(shape=4,size=4)+geom_line(linetype=1,alpha=I(1/2))+guides(col = guide_legend(title=NULL,nrow = 25))+ggtitle("PERCENTAGE OF TARGET BASES ACHIEVING X COVERAGE OR GREATER")+scale_x_continuous("COVERAGE, X",breaks=c(10,20,30))+scale_y_continuous("BASES [%]",expand = c(0, 0),lim=c(0,100))+theme(plot.title=element_text(face="bold",vjust=2),axis.title.x=element_text(vjust=-0.5),axis.title.y=element_text(vjust=0.5),panel.background=element_rect(fill="white"),panel.grid.major.x=element_line(colour="grey95"),panel.grid.major.y=element_line(colour="grey90"))

dev.off()

