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

    sampleTable = read.table(design.file)
    cds = newCountDataSetFromHTSeqCount(sampleTable, directory = counts.table)
    cds

}

#remove lines with all 0 counts
use = (rowSums ( counts ( cds )) > 0)
cds = cds[ use, ]
cds

########################################################################
#estimate size factor, dispersion and run binomial test for original data
########################################################################
cds = estimateSizeFactors( cds )

if (disp.mode == "maximum") {
    cds = estimateDispersions( cds, method="pooled", sharingMode="maximum", fitType="parametric")
} else if (disp.mode == "gene") {
    cds = estimateDispersions( cds,  method="pooled", sharingMode="gene-est-only", fitType="parametric")
} else if (disp.mode == "fit") {
    cds = estimateDispersions( cds, method="blind", sharingMode="fit-only", fitType="parametric")
}

res = nbinomTest( cds, "1", "2" )

#########################################################################
#filter genes according to filtering statistics and cutoff 
#########################################################################
if (filtering.stats == "mean") {
    rs = rowMeans(counts(cds))
} else if (filtering.stats == "min") {
    rs = rowMin(counts(cds))
} else if (filtering.stats == "max") {
    rs = rowMax(counts(cds))
} else if (filtering.stats == "sd") {
    rs = rowSds(counts(cds))
}

#plot counts rank vs p-value to assess filtering statistics
#cap -log(p-value) to 20
log = -log10(res$pval)
log[log > 20] = 20
png( file = paste( results.dir, "rank_vs_pval_plot.png", sep="/" ) )
    plot(rank(rs)/length(rs), log, ylab = "-log10(p-value)")
dev.off()

use = (rs > quantile(rs, probs=filtering.cutoff))
cdsFilt = cds[ use, ]
cdsFilt

########################################################################
#estimate dispersion and run binomial test for filtered data
########################################################################

if (disp.mode == "maximum") {
    cdsFilt = estimateDispersions( cdsFilt, method="pooled", sharingMode="maximum", fitType="parametric")
} else if (disp.mode == "gene") {
    cdsFilt = estimateDispersions( cdsFilt,  method="pooled", sharingMode="gene-est-only", fitType="parametric")
} else if (disp.mode == "fit") {
    cdsFilt = estimateDispersions( cdsFilt, method="blind", sharingMode="fit-only", fitType="parametric")
}

write.table( fData( cdsFilt ), file = paste( results.dir, "dispersion.filt.txt", sep="/"), sep = "\t" )

#plot dispersion
png( file = paste( results.dir, "dispersion_plot.png", sep="/" ) )
    plotDispEsts( cdsFilt )
dev.off()

#calculate differential expression
resFilt = nbinomTest( cdsFilt, "1", "2" )
write.table( resFilt[ order(resFilt$pval), ], file = paste( results.dir, "nbinom.filt.txt", sep="/"), sep = "\t" )

#select set of differentially expressed genes at FDR = 10% 
resFiltSig = resFilt[ resFilt$padj < 0.1, ]
write.table( resFiltSig[ order(resFiltSig$pval), ], file = paste( results.dir, "nbinom.filt.sig.txt", sep="/"), sep = "\t" )

#plot fold change
png( file = paste( results.dir, "FC_plot.png", sep="/" ) )
    plotMA( resFilt )
dev.off()

#make a table comparing the power of tests on original and filtered data
resFiltForComparison = rep(+Inf, length(res$padj))
resFiltForComparison[use] = resFilt$padj
table.filt = addmargins(table(`no filtering` =  res$padj < 0.1, `with filtering` = resFiltForComparison < 0.1 ))
write.table(table.filt, file = paste( results.dir, "original_vs_filt.tsv", sep="/"), sep = "\t", col.names = NA, row.names = c("filt FALSE", "filt TRUE", "filt sum"))

#########################################################################
#plot GC content and length bias for gene-level counts and genome-wide distribution of genes
#########################################################################
geneLevelData = counts(cdsFilt)

conditionData = data.frame( row.names = colnames(geneLevelData), conditions = conditions(cdsFilt)[colnames(geneLevelData)] )

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
    plot(length.hist$mids, log(length.hist$counts), type='h',  col = "grey", xlim=c(0,30000), xlab = "", ylab = "", axes = F)
    mtext("Genome-wide number of genes (log) binned by length", col = "grey", side = 4, line = 3)
    axis(4, col = "grey", col.axis = "grey")
dev.off()

#######################################################################################
#boxplot p-value for DE analysis binned by GC-content or gene length
#######################################################################################
#cap -log(p-value) to 20
log = -log10(resFilt$pval)
log[log > 20] = 20
png( file = paste( results.dir, "DEanalysis_gc_bias.png", sep="/" ) )
    biasBoxplot(log, fData(data)[resFilt$id,]$gc, xlab = "Gene GC-content, %", ylab = "-log10(p-value)", main = "DE analysis: p-value binned by gene GC-content")
dev.off()

png( file = paste( results.dir, "DEanalysis_length_bias.png", sep="/" ) )
    biasBoxplot(log, fData(data)[resFilt$id,]$length, xlab = "Gene length, bp", ylab = "-log10(p-value)", main = "DE analysis: p-value binned by gene length")
dev.off()

######################################################################################
#plot GC content and length bias for DE genes
######################################################################################
hist.gc =  hist(fData(data)$gc, plot = F)
hist.gc.DEseq = hist(fData(data)[resFiltSig$id,]$gc, breaks = hist.gc$breaks, plot = F)

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
hist.length.DEseq = hist(fData(data)[resFiltSig$id,]$length, breaks = hist.length$breaks, plot = F)

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
