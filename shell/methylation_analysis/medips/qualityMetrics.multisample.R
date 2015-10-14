library(MEDIPS)
library(BSgenome.Hsapiens.UCSC.hg19)

#reference genome
BSgenome="BSgenome.Hsapiens.UCSC.hg19"
#remove duplicate reads
uniq=TRUE
#extend reads
extend=#extendReads
#do not shift reads
shift=0
#set window size to 100nt
ws=#windowSize
#paired reads 
paired=#pairedReads

tmp.dir <- "#tmpDir"
sample.info.file <- "#sampleInfo"
medips.dir <- "#medipsDir"

#read sample info
sample.info <- read.delim(sample.info.file, header=F, skip=1, col.names=c("name", "category"), as.is=T)

#make list of Input sets
list <- NULL

#calculate pairwise Pearson correlations between provided sets

for ( sample in sample.info$name ) { 

	list <- c(list, MEDIPS.createSet(
							file = paste(tmp.dir, "/", sample , ".bam", sep=""), 
							BSgenome = BSgenome, 
							uniq = uniq, 
							extend = extend, 
							shift = shift, 
							window_size = ws, 
							paired = paired, 
							sample_name = sample))

}

cor.matrix <- MEDIPS.correlation(MSets = list)
write.table(cor.matrix, file = paste(medips.dir, "/", "correlation_matrix.tsv", sep=""), sep = "\t", row.names = sample.info$name, col.names = sample.info$name)

metrics <- read.delim(file = paste(medips.dir, "/", "correlation_matrix.tsv", sep=""), as.is=T, row.names=1, header=FALSE, skip=1)
colnames(metrics) <- rownames(metrics)
hc <- hclust(dist(1-metrics),method="average")

png(file = paste(medips.dir, "/", "correlation_plot.png", sep=""))
	plot(hc)
dev.off()
