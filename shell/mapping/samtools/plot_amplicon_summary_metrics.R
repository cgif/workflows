library("ggplot2")
library("reshape")


metrics.dir <- "#resultsDir"
sample.interval.summary <- "#sampleIntervalSummary"

#plot target coverage
coverage <- read.delim(paste(metrics.dir, "/", sample.interval.summary, sep=""), as.is=T, row.names=1)

#order samples by median amplicon coverage
sample.order <- names(sort(apply(coverage,1,median), decreasing=T))
target.order <- names(sort(apply(coverage,2,median), decreasing=T))

coverage.n <- as.matrix(coverage[sample.order,])

coverage.n.melt <- melt.array(coverage.n, varnames = names(dimnames(coverage)))

#categorise coverage values
coverage.n.melt.categorised <- cut(coverage.n.melt[,3],breaks = c(0,10,100,500,1000,5000,Inf),right = FALSE)
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


