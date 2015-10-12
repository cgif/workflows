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

    sampleTable = read.table(design.file)
    cds = newCountDataSetFromHTSeqCount(sampleTable, directory = counts.table)
    cds

}

#remove lines with all 0 counts
use = (rowSums ( counts ( cds )) > 0)
cds = cds[ use, ]
cds

#########################################################################
#estimate size factor and dispersion
#########################################################################
cds = estimateSizeFactors( cds )
write.table(counts( cds, normalized=TRUE ) , file = paste( results.dir, "counts.norm.txt", sep="/"), sep = "\t")

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
write.table( res[ order(res$pval), ], file = paste( results.dir, "nbinom.txt", sep="/"), sep = "\t" )

#select set of differentially expressed genes at FDR = 10% 
resSig = res[ res$padj < 0.1, ]
write.table( resSig[ order(resSig$pval), ], file = paste( results.dir, "nbinom.sig.txt", sep="/"), sep = "\t" )

#plot fold change
png( file = paste( results.dir, "FC_plot.png", sep="/" ) )
    plotMA( res )
dev.off()

#plot p-val histogram
png( file = paste( results.dir, "pval_hist.png", sep="/" ) )
    hist( res$pval, breaks=100 )
dev.off()

#######################################################################
#pick up filtering parameters
#######################################################################
filterChoices = data.frame(
    'mean' = rowMeans(counts(cds)),
    'geneID' = as.numeric(sub("[:alpha:]*", "", rownames(res))),
    'min' = rowMin(counts(cds)),
    'max' = rowMax(counts(cds)),
    'sd' = rowSds(counts(cds))
)
theta = seq(from=0, to=0.8, by=0.02)
rejChoices = sapply(filterChoices, function(f) filtered_R(alpha=0.1, filter=f, test=res$pval, theta=theta, method="BH"))
myColours = brewer.pal(ncol(filterChoices), "Set1")
png( file = paste( results.dir, "filtering_param_plot.png", sep="/" ) )
    matplot(theta, rejChoices, type="l", lty=1, col=myColours, lwd=2, xlab="filtering quantile", ylab="number of significant tests (differentially expressed genes)")
    legend("bottomleft", legend=colnames(filterChoices), fill=myColours)
    grid(22)
dev.off()

#########################################################################
#sample clustering and PCA analysis
#########################################################################
#sample clustering following variance stabilizing transformation for the most highly expressed 10000 genes
cdsBlind = newCountDataSetFromHTSeqCount(sampleTable, directory = counts.table)
cdsBlind = estimateSizeFactors( cdsBlind )

if (disp.mode == "maximum") {
    cdsBlind = estimateDispersions( cdsBlind, method="blind", sharingMode="maximum", fitType="parametric")
} else if (disp.mode == "gene") {
    cdsBlind = estimateDispersions( cdsBlind,  method="blind", sharingMode="gene-est-only", fitType="parametric")
} else if (disp.mode == "fit") {
    cdsBlind = estimateDispersions( cdsBlind, method="blind", sharingMode="fit-only", fitType="parametric")
}

useBlind = (rowSums ( counts ( cdsBlind )) > 0)
cdsBlind = cdsBlind[ use, ]
vsd = varianceStabilizingTransformation( cdsBlind )

colnames(exprs(vsd)) =  with(pData(vsd), paste(colnames(exprs(vsd)), condition, sep=":"))
hmcol = colorRampPalette(brewer.pal(9, "GnBu"))(100)
select = order(rowMeans(counts(cds)), decreasing=TRUE)[1:10000]

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
    print(plotPCA(vsd, intgroup=c("condition"), ntop = 10000))
dev.off()

#########################################################################
#plot GC content and length bias for gene-level counts and genome-wide distribution of genes
#########################################################################
geneLevelData = counts(cds)

conditionData = data.frame( row.names = colnames(geneLevelData), conditions = conditions(cds)[colnames(geneLevelData)] )

gcData = read.table( GCcont.file, row.names=1 )
lengthData = read.table( length.file, row.names=1 )
features = data.frame( row.names = rownames(geneLevelData), gc = gcData[rownames(geneLevelData),], length = lengthData[rownames(geneLevelData),] )

data = newSeqExpressionSet(exprs = as.matrix(geneLevelData), featureData = features, phenoData = conditionData )

gc.hist = hist(fData(data)$gc, plot = F)

png( file = paste( results.dir, "input_gc_bias.png", sep="/" ) )
    par(mar=c(c(5, 4, 4, 4) + 0.1))
    biasPlot(data, "gc", log=TRUE, xlim=c(20,80), ylim=c(0,10), xlab = "Gene GC-content, %", ylab = "Gene-level counts (log)", main = "Lowess regression of the gene-level counts (log) on GC-content", legend = F)
    legend("topleft", legend = c("control","condition"), fill = c("black","red"))
    par(new=TRUE)
    plot(gc.hist$mids, log(gc.hist$counts), type='h', col = "grey", xlim=c(20,80), xlab = "", ylab = "", axes = F)
    mtext("Genome-wide number of genes (log) binned by GC-content", col = "grey", side = 4, line = 3)
    axis(4, col = "grey", col.axis = "grey")
dev.off()

length.hist = hist(fData(data)$length, breaks=100, plot = F)

png( file = paste( results.dir, "input_length_bias.png", sep="/" ) )
    par(mar=c(c(5, 4, 4, 4) + 0.1))
    biasPlot(data, "length", log=TRUE, xlim=c(0,30000), ylim=c(0,10), xlab = "Gene length, bp", ylab = "Gene-level counts (log)", main = "Lowess regression of the gene-level counts (log) on gene length", legend = F)
    legend("topleft", legend = c("control","condition","counts ~ gene length"), fill = c("black","red","blue"))
    lines(seq(0,30000),log(0.1*seq(0,30000)),col="blue",lty=2)
    lines(seq(0,30000),log(0.05*seq(0,30000)),col="blue",lty=2)
    lines(seq(0,30000),log(0.01*seq(0,30000)),col="blue",lty=2)
    par(new=TRUE)
    plot(length.hist$mids, log(length.hist$counts), type='h', col = "grey", xlim=c(0,30000), xlab = "", ylab = "", axes = F)
    mtext("Genome-wide number of genes (log) binned by length", col = "grey", side = 4, line = 3)
    axis(4, col = "grey", col.axis = "grey")
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
#plot GC content and length bias for DE genes
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
