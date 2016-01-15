library( "DESeq" )
require( "EDASeq" )
library("RColorBrewer")
library("gplots")
library("genefilter")

results.dir = "#resultsDir"
counts.table = "#countsTable"
design.file = "#designFile"
htseq = "#htSeq"
disp.mode = "#dispertionMode"
filtering.stats= "#filteringStats"
filtering.cutoff= #filteringCutoff
GCcont.file = "#CGcontent"
length.file = "#geneLength"

#########################################################################
#generate new data set
#########################################################################
if (htseq == "F") {

    countTable = read.table( counts.table, header=TRUE, row.names=1 )
    head(countTable)
    conditionTable = read.table( design.file, header=TRUE, row.names=1 )
    conditionTable
    conditionTableSync = data.frame( row.names = colnames(countTable), condition = conditionTable[colnames(countTable),] )
    conditionTableSync
    cds = newCountDataSet( countTable, conditionTableSync$condition )
    cds

} else if (htseq == "T") {

    sampleTable = read.table(design.file, header=TRUE)
    cds = newCountDataSetFromHTSeqCount(sampleTable, directory = counts.table)
    cds

}

#remove lines with all 0 counts
use = (rowSums ( counts ( cds )) > 0)
cds = cds[ use, ]
cds

geneLevelData = counts(cds)

#filter genes according to filtering statistics and cutoff 
if (filtering.stats == "mean") {
    rs = rowMeans(geneLevelData)
} else if (filtering.stats == "min") {
    rs = rowMin(geneLevelData)
} else if (filtering.stats == "max") {
    rs = rowMax(geneLevelData)
} 

use = (rs > quantile(rs, probs=filtering.cutoff))
geneLevelData = geneLevelData[ use, ]

conditionData = data.frame( row.names = colnames(geneLevelData), conditions = conditions(cds)[colnames(geneLevelData)] )

gcData = read.table( GCcont.file, row.names=1 )
lengthData = read.table( length.file, row.names=1 )
features = data.frame( row.names = rownames(geneLevelData), gc = gcData[rownames(geneLevelData),], length = lengthData[rownames(geneLevelData),] )

data = newSeqExpressionSet(exprs = as.matrix(geneLevelData), featureData = features, phenoData = conditionData )

#####################################################################################
#run within sample normalization over GC-content and gene length
#####################################################################################

data.GCnorm = withinLaneNormalization(data,"gc",which="full")
data.GCnorm.lengthNorm  = withinLaneNormalization(data.GCnorm,"length",which="full")

#####################################################################################
#run between sample normalization over sequencing depth
#####################################################################################

data.GCnorm.lengthNorm.depthNorm = betweenLaneNormalization(data.GCnorm.lengthNorm, which="full")

#####################################################################################
#plot bias before and after normalization
#####################################################################################

#GC-content bias
png( file = paste( results.dir, "norm_GC_bias.png", sep="/" ) )
    biasPlot(data.GCnorm.lengthNorm.depthNorm, "gc", log=TRUE, xlim=c(20,80), ylim=c(0,10), xlab = "Gene GC-content, %", ylab = "Gene-level counts (log)", main = "GC-content bias before and after normalization", legend = F)
    legend("topleft", legend = c("control","condition"), fill = c("black","red"))
    par(new=TRUE)
    biasPlot(data, "gc", log=TRUE, xlim=c(20,80), ylim=c(0,10), xlab = "", ylab = "", main = "", legend = F, lty = 3)
dev.off()

#gene length bias
png( file = paste( results.dir, "norm_length_bias.png", sep="/" ) )
    biasPlot(data.GCnorm.lengthNorm.depthNorm, "length", log=TRUE, xlim=c(0,30000), ylim=c(0,10), xlab = "Gene length, bp", ylab = "Gene-level counts (log)", main = "Gene length bias before and after normalization", legend = F)
    legend("topleft", legend = c("control","condition"), fill = c("black","red"))
    par(new=TRUE)
    biasPlot(data, "length", log=TRUE, xlim=c(0,30000), ylim=c(0,10), xlab = "", ylab = "", main = "", legend = F, lty = 3)
dev.off()

#sequencing depth bias
png( file = paste( results.dir, "norm_read_depth_bias.png", sep="/" ) )
    par(mfrow=c(1,2), oma=c(0,0,1,0))
    boxplot(data, lof = F, main = "Before normalization")
    boxplot(data.GCnorm.lengthNorm.depthNorm, main = "After normalization")
dev.off()

#########################################################################
#estimate size factor and dispersion
#########################################################################

cds = as(data.GCnorm.lengthNorm.depthNorm,"CountDataSet")
sizeFactors(cds) = c(1)
use = (rowSums (counts(cds)) > 0)
cds = cds[use,]

cds = estimateSizeFactors( cds )

if (disp.mode == "maximum") {
    cds = estimateDispersions( cds, method="pooled", sharingMode="maximum", fitType="parametric")
} else if (disp.mode == "gene") {
    cds = estimateDispersions( cds,  method="pooled", sharingMode="gene-est-only", fitType="parametric")
} else if (disp.mode == "fit") {
    cds = estimateDispersions( cds, method="blind", sharingMode="fit-only", fitType="parametric")
}

#plot dispersion
png( file = paste( results.dir, "dispersion_plot.png", sep="/" ) )
    plotDispEsts( cds )
dev.off()

#########################################################################
#calculate differential expression 
#########################################################################

res = nbinomTest( cds, "1", "2" )
write.table( res[ order(res$pval), ], file = paste( results.dir, "nbinom.filt.norm", sep="/"), sep = "\t" )

resSig = res[ res$padj < 0.1, ]
write.table( resSig[ order(resSig$pval), ], file = paste( results.dir, "nbinom.filt.norm.sig", sep="/"), sep = "\t" )

#plot fold change
png( file = paste( results.dir, "FC_plot.filt.png", sep="/" ) )
    plotMA( res )
dev.off()

#########################################################################
#sample clustering and PCA analysis
#########################################################################
#sample clustering following variance stabilizing transformation for the most highly expressed 10000 genes
cdsBlind = as(data.GCnorm.lengthNorm.depthNorm,"CountDataSet")
sizeFactors(cdsBlind) = c(1)

if (disp.mode == "maximum") {
    cdsBlind = estimateDispersions( cdsBlind, method="blind", sharingMode="maximum", fitType="parametric")
} else if (disp.mode == "gene") {
    cdsBlind = estimateDispersions( cdsBlind,  method="blind", sharingMode="gene-est-only", fitType="parametric")
} else if (disp.mode == "fit") {
    cdsBlind = estimateDispersions( cdsBlind, method="blind", sharingMode="fit-only", fitType="parametric")
}

useBlind = (rowSums ( counts ( cdsBlind )) > 0)
cdsBlind = cdsBlind[ useBlind, ]
vsd = varianceStabilizingTransformation( cdsBlind )

colnames(exprs(vsd)) =  with(pData(vsd), paste(colnames(exprs(vsd)), condition, sep=":"))
hmcol = colorRampPalette(brewer.pal(9, "GnBu"))(100)
select = order(rowMeans(counts(cds)), decreasing=TRUE)[1:1000]

png( file = paste( results.dir, "vst_heatmap.png", sep="/" ) )
    heatmap.2(exprs(vsd)[select,], col = hmcol, trace="none", margin=c(10, 6))
dev.off()

#plot Eucludian distances between samples
dists = dist( t( exprs(vsd) ) )
mat = as.matrix( dists )

png( file = paste( results.dir, "distances_heatmap.png", sep="/" ) )
    heatmap.2(mat, trace="none", col = hmcol, margin=c(13, 13))
dev.off()

#plot PCA results
png( file = paste( results.dir, "PCA.png", sep="/" ) )
    print(plotPCA(vsd, intgroup=c("condition"), ntop = 1000))
dev.off()

#######################################################################################
#boxplot p-value for DE analysis binned by GC-content or gene length
#######################################################################################
#cap -log(p-value) to 20
log = -log10(res$pval)
log[log > 20] = 20
png( file = paste( results.dir, "DEanalysis_gc_bias.png", sep="/" ) )
    biasBoxplot(log, fData(data)[res$id,]$gc, xlab = "Gene GC-content, %", ylab = "-log10(p-value)", main = "DE analysis: p-value binned by gene GC-content")
dev.off()

png( file = paste( results.dir, "DEanalysis_length_bias.png", sep="/" ) )
    biasBoxplot(log, fData(data)[res$id,]$length, xlab = "Gene length, bp", ylab = "-log10(p-value)", main = "DE analysis: p-value binned by gene length")
dev.off()

######################################################################################
#plot persentage of DE genes vs GC content or gene length
######################################################################################
hist.gc =  hist(fData(data)$gc, plot = F)
hist.gc.DEseq = hist(fData(data)[resSig$id,]$gc, breaks = hist.gc$breaks, plot = F)

png( file = paste( results.dir, "DEgenes_gc_bias.png", sep="/" ) )
    par(mar=c(c(5, 4, 4, 4) + 0.1))
    plot(hist.gc$mids, (hist.gc.DEseq$counts/hist.gc$counts)*100, type = "l", xlim = c(20,80), ylim=c(0,5), xlab = "Gene GC-content, %", ylab = "DE genes, %", main = "Significant DE genes binned by gene GC-content")
    legend("topleft", legend = c("DE genes, %","DE genes, #","GC bias in NGS data"), fill = c("black","red","blue"))
    par(new=TRUE)
    plot(hist.gc$mids, hist.gc.DEseq$counts, type = "h", xlim = c(20,80), axes = FALSE, xlab = "",  ylab = "", col = "red")
    mtext("DE genes, #", side = 4, col = "red", line = 3) 
    axis(4, col = "red", col.axis = "red")
    par(new=TRUE)
    biasPlot(data, "gc", log=TRUE, xlim=c(20,80), ylim=c(0,10), axes=FALSE, xlab = "", ylab = "", legend = F, lty = 3, col = "blue")
dev.off()

hist.length = hist(fData(data)$length, breaks = 100, plot = F)
hist.length.DEseq = hist(fData(data)[resSig$id,]$length, breaks = hist.length$breaks, plot = F)

png( file = paste( results.dir, "DEgenes_length_bias.png", sep="/" ) )
    par(mar=c(c(5, 4, 4, 4) + 0.1))
    plot(hist.length$mids, (hist.length.DEseq$counts/hist.length$counts)*100, xlim=c(0,30000), ylim=c(0,5), type = "l", xlab = "Gene length, bp", ylab = "DE genes, %", main = "Significant DE genes binned by gene length")
    legend("topleft", legend = c("DE genes, %","DE genes, #","gene length bias in NGS data"), fill = c("black","red","blue"))
    par(new=TRUE)
    plot(hist.length$mids, hist.length.DEseq$counts, type = "h", xlim = c(0,30000), axes = FALSE, xlab = "",  ylab = "", col = "red")
    mtext("DE genes, #", side = 4, col = "red", line = 3) 
    axis(4, col = "red", col.axis = "red")
    par(new=TRUE)
    biasPlot(data, "length", log=TRUE, xlim=c(0,30000), ylim=c(0,10), axes=FALSE, xlab = "", ylab = "", legend = F, lty = 3, col = "blue")
dev.off()

