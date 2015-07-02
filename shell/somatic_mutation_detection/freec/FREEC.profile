[general]

#Parameter chrLenFile is required.
chrLenFile=#chrLenFile

#Parameter ploidy is required.
ploidy=#polidy

sex=#sex

#Either coefficientOfVariation or window must be specified.
window=#window
step=#step
minCNAlength=1

#The closer "breakPointThreshold" is to 0, the more breakpoints will be called. 
breakPointThreshold=#breakPoint

#make a separate fragment of the ambiguous regions (poly-N or low mappability regions between 
#two different copy number values) and do not assign any copy number to this region at all 
breakPointType=4

contaminationAdjustment=TRUE

#simply model "sample RC ~ Control RC" 
forceGCcontentNormalization=0

#length of pre-telomeric and pre-centromeric regions: Control-FREEC will not output small CNAs and LOH found within these regions
telocentromeric=50000

maxThreads=6

outputDir=#resultsDir

samtools=/apps/samtools/#samtoolsVersion/bin/samtools


[sample]

#Either mateFile or mateCopyNumberFile must be specified.
mateFile=#tumorBam

#Parameters inputFormat and mateOrientation are required if mateFile is specified.
inputFormat=BAM

#use "mateOrientation=0" for sorted .SAM and .BAM
mateOrientation=0


[control]

mateFile=#normalBam

inputFormat=BAM

mateOrientation=0

