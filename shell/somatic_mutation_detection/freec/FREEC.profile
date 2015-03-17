[general]

#Parameter chrLenFile is required.
chrLenFile=#chrLenFile

#Parameter ploidy is required.
ploidy=2

#Either coefficientOfVariation or window must be specified.
window=400
step=200
minCNAlength=1

#The closer "breakPointThreshold" is to 0, the more breakpoints will be called. 
breakPointThreshold=0.6

#make a separate fragment of the ambiguous regions (poly-N or low mappability regions between 
#two different copy number values) and do not assign any copy number to this region at all 
breakPointType=4

contaminationAdjustment=TRUE

#simply model "sample RC ~ Control RC" 
forceGCcontentNormalization=0

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

