#define variables

sample <- "#sampleName"
output.dir <- "#outputDir"
read.statistics.file <- "#readStatistics"
nohit.file <- "#nohitFile"
collapsed.reads <- "#collapsedReadsFile"
expression.file <- "#expressionFile"

lengths <- read.table(read.statistics.file, skip = 1, stringsAsFactors = FALSE)
colnames(lengths) <- c("length", "unique.reads", "total.reads")

#read length distribution
file <- paste(output.dir, sample, sep = "/")
png(file = paste(file, "read_lengths.png", sep = "."))
barplot(lengths$total.reads, names.arg = lengths$length, main = paste(sample, "read length distribution", sep = " "), xlab = "read length", ylab = "count")
dev.off()


summary <- data.frame(id = sample)

nohit <- read.table(nohit.file)
names(nohit) <- c("count", "tag")
nohit.count <- sum(nohit$count)

all <- read.table(collapsed.reads)
names(all) <- c("count","read")
summary$total.reads <- sum(all$count)

summary$aligned.reads <- summary$total.reads - nohit.count
summary$unaligned.reads <- nohit.count

expr <- read.table(expression.file)
names(expr) <- c("miRNA", "count")
expr <- expr[order(-expr$count),]

summary$unique.miRNA <- nrow(expr)

total.count <- sum(expr[,2])
summary$top3.prop <- signif((sum(expr[1:3,2])/total.count), 3)
summary$top10.prop <- signif((sum(expr[1:10,2])/total.count), 3)

file <- paste(output.dir, sample, sep = "/")
write.table(x = summary, file = paste(file, "summary.txt", sep = "."), row.names = FALSE, sep = "\t")

