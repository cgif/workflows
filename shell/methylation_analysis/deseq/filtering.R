library( "DESeq" )
require( "EDASeq" )
library("RColorBrewer")
library("gplots")
library("genefilter")

counts.table = "#countsTable"
design.file = "#designFile"
filtering.stats= "#filteringStats"
filtering.cutoff= #filteringCutoff
results.dir = "#resultsDir"

#########################################################################
#generate new data set
#########################################################################
sampleTable = read.table(design.file, header = TRUE)
cds = newCountDataSetFromHTSeqCount(sampleTable, directory = counts.table)
cds

#remove lines with all 0 counts
use = (rowSums ( counts ( cds )) > 0)
cds = cds[ use, ]
cds

########################################################################
#estimate size factor, dispersion and run binomial test for original data
########################################################################
cds = estimateSizeFactors( cds )
write.table(counts( cds, normalized=TRUE ) , file = paste( results.dir, "counts.norm.tsv", sep="/"), sep = "\t")

cds = estimateDispersions( cds,  method="pooled", sharingMode="gene-est-only", fitType="local")

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

cdsFilt = estimateDispersions( cdsFilt,  method="pooled", sharingMode="gene-est-only", fitType="local")

#plot dispersion
png( file = paste( results.dir, "dispersion_plot.png", sep="/" ) )
    plotDispEsts( cdsFilt )
dev.off()

#calculate differential expression
resFilt = nbinomTest( cdsFilt, "1", "2" )
write.table( resFilt[ order(resFilt$pval), ], file = paste( results.dir, "nbinom.tsv", sep="/"), sep = "\t" )

#select set of differentially expressed genes at FDR = 10% 
resFiltSig = resFilt[ resFilt$padj < 0.1, ]
write.table( resFiltSig[ order(resFiltSig$pval), ], file = paste( results.dir, "nbinom.sig.tsv", sep="/"), sep = "\t" )

#plot fold change
png( file = paste( results.dir, "FC_plot.png", sep="/" ) )
    plotMA( resFilt )
dev.off()

#make a table comparing the power of tests on original and filtered data
resFiltForComparison = rep(+Inf, length(res$padj))
resFiltForComparison[use] = resFilt$padj
table.filt = addmargins(table(`no filtering` =  res$padj < 0.1, `with filtering` = resFiltForComparison < 0.1 ))
write.table(table.filt, file = paste( results.dir, "original_vs_filt.tsv", sep="/"), sep = "\t", col.names = NA, row.names = c("filt FALSE", "filt TRUE", "filt sum"))

