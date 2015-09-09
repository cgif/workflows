input <- "#inputTable"
output <- paste(input, ".tr", sep="")

#read table
tbl <- read.delim(input, header=T, row.names=1)

#transpose table
tbl.tr <- t(tbl)

#write transposed table
#col.names=NA writes an empty column header for the row name column
write.table(tbl.tr, file=output, quote=F, sep="\t", col.names=NA)

