library(MEDIPS)
library(BSgenome.Hsapiens.UCSC.hg19)

#reference genome
BSgenome="BSgenome.Hsapiens.UCSC.hg19"
#remove duplicate reads
uniq=TRUE
#extend reads to 300 nt
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
condition <- "#condition"

#read sample info
sample.info <- read.delim(sample.info.file, header=F, skip=1, col.names=c("name", "condition"), as.is=T)

#calculate methylation profile
list <- NULL

for ( sample in sample.info$name[sample.info$condition == condition] ) { 

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

profile <- MEDIPS.meth(
							MSet1 = list,  
							MeDIP = FALSE, 
							CNV = FALSE)

write.table(profile, file = paste(medips.dir, "/WIG/", condition, "_profile.txt", sep=""), sep = "\t")





