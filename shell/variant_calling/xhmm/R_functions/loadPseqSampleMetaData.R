#################################################
# For working with population samples:
#################################################
readPhenotypesFile <- function(file) {
	readFile = paste("grep -v '^##' ", file, sep="")
	con <- pipe(readFile)

	phenotypes = read.table(con, header=TRUE, check.names=FALSE, stringsAsFactors=FALSE, comment.char="", sep="\t", colClasses=c("#ID"="character", "FID"="character", "IID"="character"))

	return(phenotypes)
}

phenotypeDataToBinarySampleProperties <- function(phenotypes) {
	VAL_NAME_SEP = "_"

	EXCLUDE_COLS = c("#ID", "FID", "IID")
	propertyCols = setdiff(colnames(phenotypes), EXCLUDE_COLS)

	SAMPLES = as.character(phenotypes[, "#ID"])

	UNIQUE_VALS = c()
	for (col in propertyCols) {
		UNIQUE_VALS = c(UNIQUE_VALS, paste(col, VAL_NAME_SEP, unique(phenotypes[, col]), sep=""))
	}
	UNIQUE_VALS = sort(UNIQUE_VALS)

	# Add in features as binary traits:
	binaryProperties = list()
	for (col in propertyCols) {
		uniqueColVals = unique(phenotypes[, col])
		for (val in uniqueColVals) {
			valName = paste(col, VAL_NAME_SEP, val, sep="")

			binaryFeatures = 1 * (phenotypes[, col] == val)
			names(binaryFeatures) = SAMPLES
			binaryProperties[[valName]] = binaryFeatures
		}
	}

	return(binaryProperties)
}




#################################################
# For working with family pedigrees:
#################################################
readPedigreeFile <- function(file) {
	pedigree = read.table(file, header=TRUE, check.names=FALSE, stringsAsFactors=FALSE, comment.char="", sep="\t", colClasses=c("#ID"="character", "FID"="character", "IID"="character", "PAT"="character", "MAT"="character"))
	children = pedigree[, "#ID"]

	fathers = pedigree[, "PAT"]
	names(fathers) = children

	mothers = pedigree[, "MAT"]
	names(mothers) = children

	sex = as.numeric(pedigree[, "SEX"])
	names(sex) = children

	return(list(pedigree=pedigree, children=children, fathers=fathers, mothers=mothers, sex=sex))
}

pedigreeDataToBinarySampleProperties <- function(pedigreeData) {
	sex = PEDIGREE_DATA[["sex"]]

	children = PEDIGREE_DATA[["children"]]
	fathers = PEDIGREE_DATA[["fathers"]]
	mothers = PEDIGREE_DATA[["mothers"]]

	samples = c(children, fathers, mothers)

	# 1=male; 2=female
	sexMap = c(sex, rep(1, length(fathers)), rep(2, length(mothers)))
	names(sexMap) = samples

	childrenMap = c(rep(1, length(children)), rep(0, length(fathers)), rep(0, length(mothers)))
	names(childrenMap) = samples

	fathersMap = c(rep(0, length(children)), rep(1, length(fathers)), rep(0, length(mothers)))
	names(fathersMap) = samples

	mothersMap = c(rep(0, length(children)), rep(0, length(fathers)), rep(1, length(mothers)))
	names(mothersMap) = samples

	return(list(sex=sexMap, children=childrenMap, fathers=fathersMap, mothers=mothersMap))
}