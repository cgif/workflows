#!/bin/bash

#PBS -l walltime=72:00:00
#PBS -l select=1:ncpus=1:mem=10gb

#PBS -M cgi@imperial.ac.uk
#PBS -m bea
#PBS -j oe

#PBS -q pqcgi

# load modules

module load R/3.1.0
module load samtools/0.1.18

NOW="date +%Y-%m-%d%t%T%t"

# define variables
R_FUNCTIONS=#Rfunctions
R_SCRIPT=#Rscript
TARGET=#target
BAM_LIST=#BamList
RESULTS_DIR=#resultsFolder
CHROMOSOMES=#chromosomes

R_SCRIPT=${R_SCRIPT}.${CHROMOSOMES}.R

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

if [[ $CHROMOSOMES == "autosomes" ]]; then
	grep -Pv '^[X|Y]\s' $TARGET > $TMPDIR/target.bed
elif [[ $CHROMOSOMES == "X_chromosome" ]]; then
	grep -P '^X\s' $TARGET > $TMPDIR/target.bed
else
	echo "illegal value for chromosomes $CHROMOSOMES"
	exit 1
fi

echo "
source('$R_FUNCTIONS')

target.bed <- '$TMPDIR/target.bed'
bam.files <- '$TMPDIR/bam.list'
exon.counts.file <- '$TMPDIR/exon.counts.${CHROMOSOMES}.Rdata'

get.exon.counts(target.bed = target.bed,
                bam.files = bam.files,
		exon.counts.file = exon.counts.file)
	
" > $R_SCRIPT

chmod 770 $R_SCRIPT

echo "`${NOW}`parsing bam files to get exon counts"
R CMD BATCH --no-save --no-restore $R_SCRIPT ${R_SCRIPT}.log

chmod 660 ${R_SCRIPT}.log

cp $TMPDIR/*.Rdata $RESULTS_DIR
chmod 660 $RESULTS_DIR/*

ls -al

echo "`${NOW}`done"

