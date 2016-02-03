library(gplots)
library(RColorBrewer)

input <- #input
output <- #output

#for testing
#input <- "/project/tgu/results/mergetag/project/aitman_fh/multisample/2013-06-12/aitman_fh.2013-06-12.perTargetCoverage.mean_coverage"
#output <- "/project/tgu/results/mergetag/project/aitman_fh/multisample/2013-06-12/aitman_fh.2013-06-12.perTargetCoverage.mean_coverage.png"

target.names <- as.vector(read.delim(pipe(paste("head -n 1",input)), as.is=T, header=F))
target.names <- target.names[c(2:length(target.names))]

coverage <- read.table(input, header=T, sep="\t", row.names=1)
coverage_matrix <- as.matrix(coverage)

#shorten target names
for(i in 1:length(target.names)){

      tokens <- strsplit(as.character(target.names[i]), "\\[")
      name <- tokens[[1]][1]
      tokens <- strsplit(name, "\\|")
      name <- tokens[[1]][1]
      target.names[i] <- name
      
}

png(filename=output, width=1200, height=length(coverage[,1])*10)

heatmap.2(coverage_matrix, scale="none", Rowv=NA, Colv=NA,
							lhei = c(0.1,0.9),
								labCol=target.names,
											breaks=seq(0,10000,1000),
                dendrogram = "none",  ## to suppress warnings
                margins=c(15,10), cexRow=1.0, cexCol=1.0, key=TRUE, keysize=0.5,
                trace="none", xlab="target", ylab="sample")

dev.off()
