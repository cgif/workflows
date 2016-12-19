#!/bin/bash

#PBS -l walltime=72:00:00
#PBS -l select=1:ncpus=1:mem=20gb:tmpspace=#tmpSpacegb

#PBS -M igf@imperial.ac.uk
#PBS -m ea
#PBS -j oe

#PBS -q pqcgi

# load modules

module load R/#Rversion

NOW="date +%Y-%m-%d%t%T%t"

# define variables
R_FUNCTIONS=#Rfunctions
R_SCRIPT=#Rscript
TARGET=#target
BAM_LIST=#BamList
RESULTS_DIR=#resultsFolder
SUMMARY_SCRIPT_PATH=#summaryScriptPath

echo "`$NOW`creating R script for counting reads per exon"
echo "`$NOW`R script: $R_SCRIPT"
echo "`$NOW`input: $BAM_LIST"
echo "`$NOW`target: $TARGET"
echo "`$NOW`results: $RESULTS_DIR"

echo "`$NOW`copying bam files to $TMPDIR"

while read BAM_PATH; do
	SAMPLE=`basename $BAM_PATH .bam`
	echo "`$NOW`copying $BAM_PATH to temporary space"
	cp $BAM_PATH $TMPDIR/$SAMPLE.bam
	cp $BAM_PATH.bai $TMPDIR/$SAMPLE.bam.bai
	echo "$TMPDIR/$SAMPLE.bam" >> $TMPDIR/bam.list
done < $BAM_LIST

echo "`$NOW`copying target file to $TMPDIR"

cp $TARGET $TMPDIR/target.bed

echo "
source('$R_FUNCTIONS')

target.bed <- '$TMPDIR/target.bed'
bam.files <- '$TMPDIR/bam.list'
exon.counts.file <- '$TMPDIR/exon.counts.Rdata'

get.exon.counts(target.bed = target.bed,
                bam.files = bam.files,
		exon.counts.file = exon.counts.file)

sessionInfo()
	
" > $R_SCRIPT

chmod 770 $R_SCRIPT

echo "`${NOW}`parsing bam files to get exon counts"
R CMD BATCH --no-save --no-restore $R_SCRIPT ${R_SCRIPT}.log

chmod 660 ${R_SCRIPT}.log

cp $TMPDIR/*.Rdata $RESULTS_DIR
chmod 660 $RESULTS_DIR/*

#run summary script
perl $SUMMARY_SCRIPT_PATH

ls -al

echo "`${NOW}`done"

