library("ggplot2")
library("reshape")

#for testing
#results.dir <- "/groupvol/cgi/results/aitman_eds/gatk2/2013-12-20/multisample/metrics"
#sample.interval.summary <- "aitman_eds.2013-12-20.sample_interval_summary"
#sample.cumulative.coverage.proportions <- "aitman_eds.2013-12-20.sample_cumulative_coverage_proportions"

results.dir <- "#resultsDir"
sample.interval.summary <- "#sampleIntervalSummary"
sample.cumulative.coverage.proportions <- "#sampleCumulativeCoverageProportions"

metrics.dir <- paste(results.dir, "/metrics", sep="")

#plot target coverage
coverage <- read.delim(paste(metrics.dir, "/", sample.interval.summary, sep=""), as.is=T, row.names=1)

#order samples by median amplicon coverage
sample.order <- names(sort(apply(coverage,1,median), decreasing=T))
target.order <- names(sort(apply(coverage,2,median), decreasing=T))

coverage.n <- as.matrix(coverage[sample.order,])

coverage.n.melt <- melt.array(coverage.n, varnames = names(dimnames(coverage)))

#categorise coverage values
coverage.n.melt.categorised <- cut(coverage.n.melt[,3],breaks = c(0,100,500,1000,5000,Inf),right = FALSE)
coverage.m.melt <- cbind(coverage.n.melt, coverage.n.melt.categorised)
colnames(coverage.m.melt) <- c("sample", "amplicon", "coverage.nominal", "coverage")

width <- length(colnames(coverage)) * 10
height <- length(rownames(coverage)) * 10
if(height < 500){
	height <- 500
}

png(paste(metrics.dir, "/", sample.interval.summary, ".png", sep=""), width, height)
ggplot(data = coverage.m.melt, aes(x = amplicon, y = sample)) + geom_tile(aes(fill = coverage), colour = "white") + scale_fill_brewer(palette = "YlGnBu") + theme(axis.text.x=element_text(angle=90))	+ scale_y_discrete(limits=sample.order) + scale_x_discrete(limits=target.order)
dev.off()


#plot cumulative coverage %
coverage.cum <- read.delim(paste(metrics.dir, "/", sample.cumulative.coverage.proportions, sep=""), as.is=T, row.names=1)

idx <- sort(coverage.cum[,"gte_100"], decreasing=T, index.return=T)
sample.order.coverage.cum <- rownames(coverage.cum[idx[[2]],])

coverage.cum.df <- data.frame(sample=sample.order.coverage.cum, gte0=coverage.cum[sample.order.coverage.cum,"gte_0"]*100, gte100=coverage.cum[sample.order.coverage.cum,"gte_100"]*100, gte500=coverage.cum[sample.order.coverage.cum,"gte_500"]*100)

width <- 500
png(paste(metrics.dir, "/", sample.cumulative.coverage.proportions, ".png", sep=""), width, height)
ggplot(coverage.cum.df, aes(y=sample, x=gte100)) + geom_point() + scale_y_discrete(limits=sample.order.coverage.cum) + xlab("% bases covered >= 100x")
dev.off()

