library("DESeq2")
library("pheatmap")
library("ggplot2")
library("RColorBrewer")

results.dir = "#resultsDir"
counts.table = "#countsTable"
metadata.file = "#metadataFile"
paired = "#paired"
htseq = "#htSeq"

#########################################################################
#read count data from matrix and create dds
#########################################################################

if (htseq == "FALSE") {

	countData<-read.table(counts.table, header = TRUE, row.names = 1)
	colData<-read.table(metadata.file, header = TRUE, row.names = 1)
	dds<-DESeqDataSetFromMatrix(countData = countData, colData = colData, design = ~condition)

} else if (htseq == "TRUE") {

    sampleTable = read.table(metadata.file, header=TRUE)
    dds = DESeqDataSetFromHTSeqCount(sampleTable, directory = counts.table, design = ~condition)

}

#pre-filtering removing genes/rows with 0 or 1 reads
dds <- dds[ rowSums(counts(dds)) > 1, ]

#relevel
#dds$condition <- relevel(dds$condition, ref="untreated")
#dds$condition <- droplevels(dds$condition)


#differential expression analysis using 'DESeq' function
dds <- DESeq(dds, betaPrior=FALSE)
#get normalized counts and write to file
normalizedCounts<-as.data.frame(counts(dds, normalized = TRUE))
write.table(normalizedCounts, file = paste(results.dir, "normalizedCounts.tsv", sep="/"), sep = "\t", row.names = TRUE, col.names = NA, quote = FALSE)
res <- results(dds)
#res<-results(dds,addMLE=TRUE)

#plot fold change
pdf(file = paste( results.dir, "FC_plot.pdf", sep="/" ))
plotMA(res, main="DESeq2", ylim=c(-4,4))
#plotMA(res, MLE=TRUE, main="DESeq2", ylim=c(-2,2))
dev.off()

#export results
resOrdered <- res[order(res$pvalue),]
resSig <- subset(resOrdered, pvalue < 0.05)
write.table(resOrdered, file = paste(results.dir, "resultsOrdered.tsv", sep="/"), sep = "\t",col.names = NA, quote = FALSE)
write.table(resSig, file = paste(results.dir, "resultsSig.tsv", sep="/"), sep = "\t",col.names = NA, quote = FALSE)

#data transformation and visualization
#rld <- rlog(dds)
#vsd <- varianceStabilizingTransformation(dds)

#heatmap of count matrix

#select <- order(rowMeans(counts(dds,normalized=TRUE)),decreasing=TRUE)[1:1000]
#df <- as.data.frame(colData(dds)["condition"])
#pdf("countmatrix_heatmap.pdf")
#pheatmap(assay(vsd)[select,], cluster_rows=FALSE, show_rownames=FALSE,
#cluster_cols=FALSE, annotation_col=df)
#dev.off()

#heatmap of the sample-to-sample distances

#sampleDists <- dist(t(assay(vsd)))
#sampleDistMatrix <- as.matrix(sampleDists)
#rownames(sampleDistMatrix) <- paste(vsd$condition)
#colnames(sampleDistMatrix) <- NULL
#colors <- colorRampPalette( rev(brewer.pal(9, "Blues")) )(255)
#pdf("sample_to_sample_distance_heatmap.pdf")
#pheatmap(sampleDistMatrix,
#clustering_distance_rows=sampleDists,
#clustering_distance_cols=sampleDists,
#col=colors)
#dev.off()

#principal component plot of the samples

#data <- plotPCA(vsd, intgroup=c("condition"), returnData=TRUE, ntop = 1000)
#percentVar <- round(100 * attr(data, "percentVar"))
#pdf("pca_samples.pdf")
#ggplot(data, aes(PC1, PC2, color=condition)) +
#geom_point(size=3) +
#xlab(paste0("PC1: ",percentVar[1],"% variance")) +
#ylab(paste0("PC2: ",percentVar[2],"% variance"))
#dev.off()
