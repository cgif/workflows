library("DESeq2")
library("pheatmap")
library("ggplot2")
library("RColorBrewer")

results.dir = "#resultsDir"
counts.table = "#countsTable"
metadata.file = "#metadataFile"
paired = "#paired"

#########################################################################
#read count data from matrix and create dds
#########################################################################

#sampleTable = read.table(metadata.file, header=TRUE)

#ddsHTSeq <- DESeqDataSetFromHTSeqCount(sampleTable = sampleTable,
#directory = directory,
#design= ~ condition)

dds<-DESeqDataSetFromMatrix(countData = counts.table, colData = metadata.file, design = ~condition)

#pre-filtering removing genes/rows with 0 or 1 reads
dds <- dds[ rowSums(counts(dds)) > 1, ]
#dds$condition <- relevel(dds$condition, ref="untreated")
#dds$condition <- droplevels(dds$condition)


#differential expression analysis using 'DESeq' function, contrasts?
dds <- DESeq(dds)
res <- results(dds)

#plot fold change
pdf("MAplot.pdf")
plotMA(res, main="DESeq2", ylim=c(-2,2))
dev.off()

#export results
resOrdered <- res[order(res$padj),]
resSig <- subset(resOrdered, padj < 0.1)
write.csv(as.data.frame(resSig), file="results.txt")

#data transformation and visualization
#rld <- rlog(dds)
vsd <- varianceStabilizingTransformation(dds)

#heatmap of count matrix

select <- order(rowMeans(counts(dds,normalized=TRUE)),decreasing=TRUE)[1:20]
df <- as.data.frame(colData(dds)["condition"])
pdf("countmatrix_heatmap.pdf")
pheatmap(assay(vsd)[select,], cluster_rows=FALSE, show_rownames=FALSE,
cluster_cols=FALSE, annotation_col=df)
dev.off()

#heatmap of the sample-to-sample distances
sampleDists <- dist(t(assay(vsd)))
sampleDistMatrix <- as.matrix(sampleDists)
rownames(sampleDistMatrix) <- paste(vsd$condition)
colnames(sampleDistMatrix) <- NULL
colors <- colorRampPalette( rev(brewer.pal(9, "Blues")) )(255)
pdf("sample_to_sample_distance_heatmap.pdf")
pheatmap(sampleDistMatrix,
clustering_distance_rows=sampleDists,
clustering_distance_cols=sampleDists,
col=colors)
dev.off()

#principal component plot of the samples

print(plotPCA(vsd, intgroup=c("condition"), ntop = 1000))

data <- plotPCA(vsd, intgroup=c("condition"), returnData=TRUE)
percentVar <- round(100 * attr(data, "percentVar"))
pdf("pca_samples.pdf")
ggplot(data, aes(PC1, PC2, color=condition)) +
geom_point(size=3) +
xlab(paste0("PC1: ",percentVar[1],"% variance")) +
ylab(paste0("PC2: ",percentVar[2],"% variance"))
dev.off()
