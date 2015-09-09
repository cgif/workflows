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
ws=100
#paired reads 
paired=#pairedReads

tmp.dir <- "#tmpDir"
sample.info.file <- "#sampleInfo"
medips.dir <- "#medipsDir"

#sample.info.file <- "/project/tgu/rawdata/kemp_copd2/docs/kemp_copd2.sample_groups.tsv"

#read sample info
sample.info <- read.delim(sample.info.file, header=F, skip=1, col.names=c("name", "category"), as.is=T)

#make list of Input sets
input.list <- NULL
all.list <- NULL

for ( sample in sample.info$name[sample.info$category == "input"] ) { 

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

all.list <- input.list

#make list of sets for each category and calculate methylation profile
categories <- levels(factor(sample.info$category[sample.info$category != "input"]))
for ( category in categories ) {

	category.list <- NULL

	for ( sample in sample.info$name[sample.info$category == category] ) { 

		category.list <- c(category.list, MEDIPS.createSet(
							file = paste(tmp.dir, "/", sample , ".bam", sep=""), 
							BSgenome = BSgenome, 
							uniq = uniq, 
							extend = extend, 
							shift = shift, 
							window_size = ws, 
							paired = paired, 
							sample_name = sample))

	}

	category.profile <- MEDIPS.meth(
							MSet1 = category.list,  
							ISet1 = input.list,
							MeDIP = FALSE, 
							CNV = FALSE)
	write.table(category.profile, file = paste(medips.dir, "/", category, "_profile.txt", sep=""), sep = "\t")

	all.list <- c(all.list, category.list)

}


#calculate pairwise Pearson correlations between provided sets

png(paste(medips.dir, "/", "correlation_matrix.png", sep=""))
	MEDIPS.correlation(MSets = all.list, plot = T)
dev.off()

cor.matrix <- MEDIPS.correlation(MSets = all.list)
write(cor.matrix, file = paste(medips.dir, "/", "correlation_matrix.tsv", sep=""), sep = "\t")



