# R functions to parse bam files to get exon counts and to call CNVs

################################################################
# get exon counts, required to run once for each batch of samples
################################################################

# target.bed should contain exons or target regions. If target regions are long/continuous, they must be broken into the fragments with ~300bp gaps (to avoid the same read counting in two regions)

# bam.files have to have duplicates removed

# exon.counts.file - output  Rdata file to write the results to be used in calling CNVs. 

get.exon.counts <- function (target.bed, 
                             bam.files, 
                             exon.counts.file) {

        library(ExomeDepth)

    	bam.files <- scan(bam.files, what = 'character')
    	exon.counts <- getBamCounts(bed.file = target.bed, bam.files = bam.files)
	save(exon.counts, file = exon.counts.file)

}


###########################################################
# get CNV calls, to be run for each sample
###########################################################

call.cnvs <- function (exon.counts.file,
                       target.bed,
		       test.sample,			# test sample name as in exon.counts.file object
		       ref.samples.file,		# file containing a list of reference samples, one per line
		       annotations.file,		# file with all required annotations listed
		       all.exons.output,		# output Rdata file saving the ExomeDepth object
		       cnv.calls.file,			# output text file containing CNV calls
		       summary.file,			# output summary file
			bam.suffix) {			# suffix for bam files, default '.bam'

        library(ExomeDepth)
        load(exon.counts.file)

	#convert exon counts object into dataframe 
	exon.counts.dafr <- as(exon.counts[, colnames(exon.counts)], 'data.frame')
	print(head(exon.counts.dafr))

	#remove bam file suffix from column names
	names(exon.counts.dafr) <- gsub(names(exon.counts.dafr), pattern = paste(bam.suffix), replacement = '')

        #make sure there are no chr prefix in chr name
        exon.counts.dafr$chromosome <- gsub(as.character(exon.counts.dafr$space),
                                           pattern = 'chr',
                                           replacement ='')

        
	print(head(exon.counts.dafr))

        #subset autosomes and X chromosome

        exon.counts.dafr.autosomes <- subset(exon.counts.dafr, chromosome %in% c(1:22))
	print(tail(exon.counts.dafr.autosomes))
	print(nrow(exon.counts.dafr.autosomes))
        exon.counts.dafr.X <- subset(exon.counts.dafr, chromosome == "X")
	print(tail(exon.counts.dafr.X))
	print(nrow(exon.counts.dafr.X))

        #analyse autosomes

        #analyse X chromosome

	#read a test sample
	my.test <- paste('exon.counts.dafr$', test.sample, sep = "")
	my.test <- eval(parse(text=my.test)) 

        #select samples for a reference set
	my.ref.samples <- scan(ref.samples.file, what = 'character')
	my.reference.set <- as.matrix(exon.counts.dafr[, my.ref.samples])

	my.choice <- select.reference.set(test.counts = my.test,
		     		          reference.counts = my.reference.set,
	                                  bin.length = (exon.counts.dafr$end - exon.counts.dafr$start)/1000,
					  formula = 'cbind(test, reference) ~ 1')

	#print choice of reference samples
	message("size of reference set: ", length(my.choice$reference.choice))

        print(my.choice$reference.choice)
	cat('reference set samples: ', my.choice$reference.choice, "\n")
        print(my.choice$summary.stats)

	# construct a reference set
        my.matrix <- as.matrix( exon.counts.dafr[, my.choice$reference.choice] )
	my.reference.selected <- apply(X = my.matrix,
				       MAR = 1,
				       FUN = sum )

	# call CNVs
	all.exons <- new('ExomeDepth',
			 test = my.test,
			 reference = my.reference.selected,
			 formula = 'cbind(test, reference) ~ 1')

	all.exons <- CallCNVs(x = all.exons,
			      transition.probability = 10^-3,
			      chromosome = exon.counts.dafr$space,
			      start = exon.counts.dafr$start,
			      end = exon.counts.dafr$end,
			      name = exon.counts.dafr$names)

	print(head(all.exons@CNV.calls))
	save(all.exons, file = all.exons.output)


	#annotate CNVs with target names
        targets <- read.table(target.bed, col.names = c("chromosome","start","end","name"))
        targets.GRanges <- GRanges(seqnames = targets$chromosome,
                                   IRanges(start = targets$start,end = targets$end),
                                   names = targets$name)

        all.exons <- AnnotateExtra(x = all.exons,
                                   reference.annotation = targets.GRanges,
                                   min.overlap = 0.0001,
                                   column.name = 'geneID')

	#annotate CNVs with common CNVs
	data(Conrad.hg19)
        all.exons <- AnnotateExtra(x = all.exons,
                                   reference.annotation = Conrad.hg19.common.CNVs,
                                   min.overlap = 0.5,
                                   column.name = 'Conrad.hg19')

	#annotate CNVs with data from public databases
        annotations <- read.table(annotations.file, col.names = c("db","cutoff","path"), as.is=TRUE)
	print(head(annotations))

        lines <- length(annotations$db)
	for (i in 1:lines){
	        print(annotations$db[i])
    	        features <- read.table(annotations$path[i], col.names = c("chromosome","start","end","name"))
        	features.GRanges <- GRanges(seqnames = features$chromosome,
                                            IRanges(start = features$start,end = features$end),
                                            names = features$name)

     		all.exons <- AnnotateExtra(x = all.exons,
                                           reference.annotation = features.GRanges,
                                           min.overlap = annotations$cutoff[i],
                                           column.name = annotations$db[i])
	}

	# rank CNVs by BF:
	cnv.calls <- (all.exons@CNV.calls[ order ( all.exons@CNV.calls$BF, decreasing = TRUE),])

	cnv.calls$size <- cnv.calls$end - cnv.calls$start
	cnv.calls$sample <- test.sample

	# reorder for clarity
	my.block <- c('sample', 'id', 'type', 'nexons', 'size', 'chromosome', 'start', 'end', 'BF', 'reads.expected', 'reads.observed', 'reads.ratio') 
	my.names <-  c(subset(my.block, my.block %in% names(cnv.calls)), subset(names(cnv.calls), ! (names(cnv.calls) %in% my.block) ) )

#        write.table(x = all.exons@CNV.calls, file = cnv.calls.file, quote = FALSE, row.names = FALSE)
	write.table(x = cnv.calls[,my.names], file = cnv.calls.file, quote = FALSE, row.names = FALSE, sep = "\t")

	# now collect summary statistics for this sample
	summary = data.frame(id = test.sample)

	cor.test.reference <- cor(all.exons@test, all.exons@reference)
	message('Once more: Correlation between reference and tests count is ', signif(cor.test.reference, 5))
	summary$cor.test.reference[1] <- signif(cor.test.reference, 5)
	summary$size.ref.set[1] <- length(my.choice$reference.choice)

	summary$total.cnvs[1] <- nrow(all.exons@CNV.calls)
	summary$del.count[1] <- sum(all.exons@CNV.calls$type == 'deletion')
	summary$dup.count[1] <- sum(all.exons@CNV.calls$type == 'duplication')
	summary$proportion.del[1] <- format(summary$del.count/summary$total.cnvs, digits = 3)
	summary$proportion.dup[1] <- format(summary$dup.count/summary$total.cnvs, digits = 3)
	summary$count.nonConrad[1] <- sum(is.na(all.exons@CNV.calls$Conrad.hg19))

	summary$ref.set[1] <- gsub(", ",",",toString(my.choice$reference.choice))


#	print(summary)
#	summary.file <- paste(cnv.calls.file, ".summary", sep = '')
	write.table(x = summary, file = summary.file, quote = FALSE, row.names = FALSE, sep = "\t")

#	cor.test.reference <- cor(all.exons@test, all.exons@reference)
#	message('Once more: Correlation between reference and tests count is ', signif(cor.test.reference, 5))
#	summary$cor.test.reference[1] <- signif(cor.test.reference, 5)

	print(summary)

}
