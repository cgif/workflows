library(MEDIPS)
library(BSgenome.Hsapiens.UCSC.hg19)

#reference genome
BSgenome="BSgenome.Hsapiens.UCSC.hg19"
#remove duplicate reads
uniq=TRUE
#extend reads to 300 nt
extend=200
#do not shift reads
shift=0
#set window size to 100nt
ws=100

inBam <- "#bamFile"
sample <- "#sample"
dir <- "#medipsSampleDir"

#estimating the saturation and reproducibility for obtaining full genome coverage profles
sr = MEDIPS.saturation(file = inBam, 
			BSgenome = BSgenome,						
			uniq = uniq, 
			extend = extend,
			shift = shift, 
			window_size = ws, 						
			nit = 10, 
			nrit = 3, 
			empty_bins = TRUE,
			rank = FALSE, 
			paired = TRUE)

png(paste(dir, "/", sample, ".sat_plot.png", sep=""))
	MEDIPS.plotSaturation(saturationObj = sr, main = "Saturation plot")
dev.off()

cat("#total number of reads ", sr$numberReads, "\n")
cat("#estimated correlation for artifically doubled set of reads ", sr$maxEstCor, "\n")

#analyzing the coverage of genome wide CpG sequence patterns
cr = MEDIPS.seqCoverage(file = inBam, 
			BSgenome = BSgenome, 
			uniq = uniq, 
			extend = extend,
			shift = shift, 
			pattern = "CG",
			paired = TRUE)

png(paste(dir, "/", sample, ".seqcov_pie.png", sep=""))
	MEDIPS.plotSeqCoverage(seqCoverageObj = cr, type = "pie", cov.level = c(0,1, 2, 3, 4, 5), main = "Coverage of genome wide CpG sequence patterns")
dev.off()

CpG.no.cov = (sum(cr$cov.res == 0)/length(cr$cov.res))*100
cat("#number of CpG that are not covered ", CpG.no.cov, "\n")
cat("#number of reads which do not cover CpG ", cr$numberReadsWO, "\n")

#test the enrichment of CpGs within the genomic regions covered by the given set of short reads compared to the full reference genome
er = MEDIPS.CpGenrich(file = inBam, 
			BSgenome = BSgenome,
			uniq = uniq, 
			extend = extend,
			shift = shift, 
			paired = TRUE)

cat("#observed/expected ratio of CpGs within the covered regions / same ratio within the genome", er$enrichment.score.GoGe, "\n")





