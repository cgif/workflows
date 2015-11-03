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

