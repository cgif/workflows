library( "DESeq" )
require( "EDASeq" )
library("RColorBrewer")
library("gplots")
library("genefilter")

counts.table = "#countsTable"
design.file = "#designFile"
results.dir="#resultsDir"

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

#########################################################################
#estimate size factor and dispersion
#########################################################################
cds = estimateSizeFactors( cds )
write.table(counts( cds, normalized=TRUE ) , file = paste( results.dir, "counts.norm.tsv", sep="/"), sep = "\t")

cds = estimateDispersions( cds,  method="pooled", sharingMode="gene-est-only", fitType="parametric")

#plot dispersion
png( file = paste( results.dir, "dispersion_plot.png", sep="/" ) )
    plotDispEsts( cds )
dev.off()

#########################################################################
#calculate differential expression 
#########################################################################
res = nbinomTest( cds, "1", "2" )
write.table( res[ order(res$pval), ], file = paste( results.dir, "nbinom.tsv", sep="/"), sep = "\t" )

#select set of differentially expressed genes at FDR = 10% 
resSig = res[ res$padj < 0.1, ]
write.table( resSig[ order(resSig$pval), ], file = paste( results.dir, "nbinom.sig.tsv", sep="/"), sep = "\t" )

#plot fold change
png( file = paste( results.dir, "FC_plot.png", sep="/" ) )
    plotMA( res )
dev.off()

#######################################################################
#pick up filtering parameters
#######################################################################
filterChoices = data.frame(
    'mean' = rowMeans(counts(cds)),
    'geneID' = as.numeric(sub("[:alpha:]*", "", rownames(res))),
    'min' = rowMin(counts(cds)),
    'max' = rowMax(counts(cds))
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

cdsBlind = estimateDispersions( cdsBlind,  method="blind", sharingMode="gene-est-only", fitType="parametric")

useBlind = (rowSums ( counts ( cdsBlind )) > 0)
cdsBlind = cdsBlind[ useBlind, ]
vsd = varianceStabilizingTransformation( cdsBlind )

colnames(exprs(vsd)) =  with(pData(vsd), paste(colnames(exprs(vsd)), condition, sep=":"))
hmcol = colorRampPalette(brewer.pal(9, "GnBu"))(100)
select = order(rowMeans(counts(cds)), decreasing=TRUE)[1:1000]

png( file = paste( results.dir, "vst_heatmap.png", sep="/" ) )
    heatmap.2(exprs(vsd)[select,], col = hmcol, trace="none", margin=c(10, 6))
dev.off()

#plot PCA results
png( file = paste( results.dir, "PCA.png", sep="/" ) )
    print(plotPCA(vsd, intgroup=c("condition"), ntop = 1000))
dev.off()

#plot Eucludian distances between samples
dists = dist( t( exprs(vsd) ) )
mat = as.matrix( dists )

png( file = paste( results.dir, "distances_heatmap.png", sep="/" ) )
    heatmap.2(mat, trace="none", col = hmcol, margin=c(13, 13))
dev.off()



