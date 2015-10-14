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
condition1 <- "#condition1"
condition2 <- "#condition2"

#read sample info
sample.info <- read.delim(sample.info.file, header=F, skip=1, col.names=c("name", "condition"), as.is=T)

#make list of Input sets
input.list <- NULL

for ( sample in sample.info$name[sample.info$condition == "input"] ) { 

	input.list <- c(input.list, MEDIPS.createSet(
							file = paste(tmp.dir, "/", sample , ".bam", sep=""), 
							BSgenome = BSgenome, 
							uniq = uniq, 
							extend = extend, 
							shift = shift, 
							window_size = ws, 
							paired = paired, 
							sample_name = sample))

}

#make list of sets for each category
list1 <- NULL
for ( sample in sample.info$name[sample.info$condition == condition1] ) { 

	list1 <- c(list1, MEDIPS.createSet(
							file = paste(tmp.dir, "/", sample , ".bam", sep=""), 
							BSgenome = BSgenome, 
							uniq = uniq, 
							extend = extend, 
							shift = shift, 
							window_size = ws, 
							paired = paired, 
							sample_name = sample))

}

list2 <- NULL
for ( sample in sample.info$name[sample.info$condition == condition2] ) { 

	list2 <- c(list2, MEDIPS.createSet(
							file = paste(tmp.dir, "/", sample , ".bam", sep=""), 
							BSgenome = BSgenome, 
							uniq = uniq, 
							extend = extend, 
							shift = shift, 
							window_size = ws, 
							paired = paired, 
							sample_name = sample))

}

#calculate differential methylation profile
DM.profile <- MEDIPS.meth(
							MSet1 = list1,  
							MSet2 = list2,  
							ISet1 = input.list,
							ISet2 = input.list,
							MeDIP = FALSE, 
							p.adj = bonferroni,
							diff.method = edgeR,
							CNV = FALSE)

write.table(DM.profile, file = paste(medips.dir, "/DM/", condition1, "_", condition2, "_profile.txt", sep=""), sep = "\t")
