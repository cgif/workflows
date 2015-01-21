library(edgeR)

targets.file <- "#targetsFile"
samples.file <- "#samplesFile"       #file containing 3 columns, sample ID, patient ID and status, with header 
results.dir <- "#resultsDir"
file.prefix <- "#prefix"
min.cpm <- minCPM
group.size <- groupSize


prefix <- paste(results.dir, file.prefix, sep = "/")

targets <- read.delim(file = targets.file, stringsAsFactors = FALSE)
targets

#read in data files
d <- readDGE(targets, header = FALSE)
samples <- read.table(file = samples.file, header = TRUE, stringsAsFactors = FALSE)
colnames(d) <- samples[,1]
d$samples
dim(d)

#keep tags that achieve two counts per million for at least 3 libraries
keep <- rowSums(cpm(d) > min.cpm) >= group.size
d <- d[keep,]
dim(d)

#reset library size
d$samples$lib.size <- colSums(d$counts)
d$samples

#normalization for RNA composition using TMM 
d <- calcNormFactors(d)
d$samples

#plot MDS
png(file = paste(prefix, "MDS.png", sep = "."))
plotMDS(d)
dev.off()

#plot MDS using log-counts-per-million
y <- cpm(d, prior.count=2, log=TRUE)
png(file = paste(prefix, "MDS_cpm.png", sep = "."))
plotMDS(y)
dev.off()

#make design matrix
Patient <- factor(samples[,2])
Status <- factor(samples[,3])
data.frame(Sample=colnames(d), Patient, Status)
design <- model.matrix(~Patient+Status)
rownames(design) <- colnames(d)

#estimate dispersion
d <- estimateGLMCommonDisp(d, design, verbose = TRUE)
d <- estimateGLMTrendedDisp(d, design)
d <- estimateGLMTagwiseDisp(d, design)

#plot dispersion
png(file = paste(prefix, "dispersion.png", sep = "."))
plotBCV(d)
dev.off()

#differential expression
fit <- glmFit(d, design)
lrt <- glmLRT(fit)
topTags(lrt)

#full results ordered by p value
o <- order(lrt$table$PValue)
ordered <- lrt$table[o,]
ordered$FDR <- p.adjust(ordered$PValue, method = "BH")

#merge with counts per million for each sample
merged <- merge(ordered, cpm(d), by = 0)
merged <- merged[order(merged$PValue),]
colnames(merged)[colnames(merged)=="Row.names"] <- "miRNA"
write.table( x = merged, file = paste(prefix, "DE_results.txt", sep = "."), row.names=FALSE, sep = "\t")

#number of differentially expressed genes at 10% FDR
summary(de <- decideTestsDGE(lrt, p.value=0.1))

#plot log-fold changes against log-counts per million
detags <- rownames(d)[as.logical(de)]
png(file = paste(prefix, "logFC.png", sep = "."))
plotSmear(lrt, de.tags=detags)
abline(h=c(-1, 1), col="blue")
dev.off()

