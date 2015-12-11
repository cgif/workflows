library("ggplot2")
library("reshape")


metrics.dir <- "#resultsDir"
sample.interval.summary <- "#sampleIntervalSummary"
non_overlapping.interval.summary <- "#nonOverlappingIntervalSummary"

sample.status <- "#sampleStatus"
date <- "#mergetagDate"
project <- "#projectName"

coverage <- read.delim(paste(metrics.dir, "/", sample.interval.summary, sep=""), as.is=T, row.names=1)

summary.samples <- data.frame(row.names = row.names(coverage))
summary.samples$median.cov <- apply(coverage,1,median)

#now get status of the samples

samples.status <- read.delim(sample.status, as.is=T, row.names=1, header = FALSE, col.names = c("sample", "status", "batch"))
summary.samples <- merge(summary.samples, samples.status, by = "row.names")
names(summary.samples)[names(summary.samples) == 'Row.names'] <- 'SampleID'

#generate list of samples to include and exclude
samples.include <- NULL
samples.exclude <- NULL

#for each batch
for (cur.batch in unique(samples.status$batch) ) {

	#get maximum median coverage for blank control in the batch
	max.blank <- NULL
	max.blank <- max(summary.samples[ which(summary.samples$status == 0 & summary.samples$batch == cur.batch), ][,"median.cov"])
	cut.off <- NULL
	cut.off <- 10*max.blank

	if (!exists('samples.include')) {
		samples.include <- summary.samples[ which(summary.samples$median.cov > cut.off & summary.samples$batch == cur.batch), ]
	} else {
		samples.include <- rbind(samples.include, summary.samples[ which(summary.samples$median.cov > cut.off & summary.samples$batch == cur.batch), ])
	}

	if (!exists('samples.exclude')) {
		samples.exclude <- summary.samples[ which(summary.samples$median.cov <= cut.off & summary.samples$batch == cur.batch), ]
	}else{
		samples.exclude <- rbind(samples.exclude, summary.samples[ which(summary.samples$median.cov <= cut.off & summary.samples$batch == cur.batch), ])
	}

}

#generate list of samples to include and exclude
samples.include <- samples.include[order(samples.include$batch),]
write.table(x = samples.include, file = paste(metrics.dir, "/", project, ".samples_include.tsv", sep = ""), row.names = FALSE, quote = FALSE)

samples.exclude <- samples.exclude[order(samples.exclude$batch),]
write.table(x = samples.exclude, file = paste(metrics.dir, "/", project, ".samples_exclude.tsv", sep = ""), row.names = FALSE, quote = FALSE)

samples.include$date <- date
samples.include$project <- project
write.table(samples.include[,c("SampleID", "date", "project")], file = paste(metrics.dir, "/gatk2_samples.", project, ".tsv", sep = ""), row.names = FALSE, col.names = FALSE, quote = FALSE)

#recalculate median amplicon coverage using only included samples and get regions to exclude
coverage.included <- subset(coverage, row.names = samples.include$SampleID)

coverage.included <- coverage[which(rownames(coverage) %in% samples.include$SampleID),]
target.coverage <- data.frame(row.names = colnames(coverage.included))
target.coverage$median.cov <- apply(coverage.included,2,median)
cut.off.target <- 5
targets.exclude <- subset(target.coverage, median.cov < cut.off.target)

coverage.non_overlapping <- read.delim(paste(metrics.dir, "/", non_overlapping.interval.summary, sep=""), as.is=T, row.names=1)
coverage.non_overlapping.included <- subset(coverage.non_overlapping, row.names = samples.include$SampleID)
target.non_overlapping.coverage <- data.frame(row.names = colnames(coverage.non_overlapping.included))
target.non_overlapping.coverage$median.cov <- apply(coverage.non_overlapping.included,2,median)
targets.non_overlapping.exclude <- subset(target.non_overlapping.coverage, median.cov < cut.off.target)

targets.exclude <- rbind(targets.exclude, targets.non_overlapping.exclude)

write.table(targets.exclude, file = paste(metrics.dir, "/", project, ".targets.exclude.tsv", sep = ""), col.names=FALSE, quote = FALSE)


#############################################

#plot target coverage

plot.coverage <- function (metrics.dir, sample.interval.summary) {
    coverage <- read.delim(paste(metrics.dir, "/", sample.interval.summary, sep=""), as.is=T, row.names=1, check.names=FALSE)
    
    #order samples by median amplicon coverage and add median values
    sample.order <- names(sort(apply(coverage,1,median), decreasing=T))
    target.order <- names(sort(apply(coverage,2,median), decreasing=T))

    coverage.n <- as.matrix(coverage[sample.order,])
    coverage.sample.median <- data.frame(row.names = row.names(coverage.n))
    coverage.sample.median$sample.median.cov <- apply(coverage.n,1,median)

    coverage.n.m <- cbind(coverage.sample.median, coverage.n)
    coverage.target.median <- t(as.data.frame(apply(coverage.n.m,2,median)))
    rownames(coverage.target.median) <- "target.median.cov"

    coverage.n.m <- rbind(coverage.n.m, coverage.target.median)
    coverage.n.m <- as.matrix(coverage.n.m)

    coverage.n.melt <- melt.array(coverage.n.m, varnames = names(dimnames(coverage.n.m)))

    #categorise coverage values
    coverage.n.melt.categorised <- cut(coverage.n.melt[,3],breaks = c(0,10,100,500,1000,5000,Inf),right = FALSE)
    coverage.m.melt <- cbind(coverage.n.melt, coverage.n.melt.categorised)
    colnames(coverage.m.melt) <- c("sample", "amplicon", "coverage.nominal", "coverage")

    width <- length(colnames(coverage)) * 10
    height <- length(rownames(coverage)) * 10
    if(height < 500){
	height <- 500
    }

    target.order <- c("sample.median.cov", target.order)
    sample.order <- c("target.median.cov", sample.order) 

    png(paste(metrics.dir, "/", sample.interval.summary, ".png", sep=""), width, height)
    png.plot <- ggplot(data = coverage.m.melt, aes(x = amplicon, y = sample)) + geom_tile(aes(fill = coverage), colour = "white") + scale_fill_brewer(palette = "YlGnBu") + theme(axis.text.x=element_text(angle=90,vjust=0.5))	+ scale_y_discrete(limits=sample.order) + scale_x_discrete(limits=target.order)
    print(png.plot)
    dev.off()
}

plot.coverage(metrics.dir = metrics.dir, sample.interval.summary = sample.interval.summary )
plot.coverage(metrics.dir, non_overlapping.interval.summary)




