#!/bin/bash

#PBS -l walltime=1:00:00
#PBS -l select=1:ncpus=1:mem=5gb

#PBS -M cgi@imperial.ac.uk
#PBS -m bea
#PBS -j oe

#PBS -q pqcgi

# load modules

module load R/3.1.0

NOW="date +%Y-%m-%d%t%T%t"

# define variables
R_FUNCTIONS=#Rfunctions
R_SCRIPT=#Rscript
RESULTS_DIR=#resultsFolder
TARGET=#target
PED_FILE=#pedFile
CHROMOSOMES=#chromosomes
ANNOTATIONS=#annotations
BAM_SUFFIX='#bamSuffix'

SAMPLE=`basename $RESULTS_DIR`
PROJECT_RESULTS_DIR=`dirname $RESULTS_DIR`

# select reference sets
FAMILY=`grep $SAMPLE $PED_FILE|cut -f 1|uniq`
 
if [[ $CHROMOSOMES == "autosomes" ]]; then

	#select unrelated samples for reference set
	grep -Pv "$FAMILY\t" $PED_FILE | grep -v "#" | cut -f 2 > $TMPDIR/ref.${CHROMOSOMES}.sample.list

elif [[ $CHROMOSOMES == "X_chromosome" ]]; then

	#select unrelated same sex samples for reference set
	SEX=`grep -P "$FAMILY\t$SAMPLE" $PED_FILE|cut -f 5` 
	grep -Pv "$FAMILY\t" $PED_FILE | grep -v "#" | cut -f 2,5 | grep -P "\t$SEX" | cut -f 1 > $TMPDIR/ref.${CHROMOSOMES}.sample.list

else
	echo "illegal value for chromosomes $CHROMOSOMES"
	exit 1
fi

echo "`$NOW`copying exon counts files to $TMPDIR"

cp $PROJECT_RESULTS_DIR/multisample/exon.counts.${CHROMOSOMES}.Rdata $TMPDIR/exon.counts.${CHROMOSOMES}.Rdata

echo "`$NOW`copying target and annotations list files to $TMPDIR"

if [[ $CHROMOSOMES == "autosomes" ]]; then
	grep -Pv '^[X|Y]\s' $TARGET > $TMPDIR/target.bed
elif [[ $CHROMOSOMES == "X_chromosome" ]]; then
	grep -P '^X\s' $TARGET > $TMPDIR/target.bed
else
	echo "illegal value for chromosomes $CHROMOSOMES"
	exit 1
fi

cp $ANNOTATIONS $TMPDIR/annotations.list

echo "`$NOW`writing R script"

R_SCRIPT=${R_SCRIPT}.${CHROMOSOMES}.R
chmod 770 $R_SCRIPT

echo "
source('$R_FUNCTIONS')

exon.counts.file <- '$TMPDIR/exon.counts.${CHROMOSOMES}.Rdata'
target.bed <- '$TMPDIR/target.bed'
test.sample <- '$SAMPLE'
ref.samples.file <- '$TMPDIR/ref.${CHROMOSOMES}.sample.list'
annotations.file <- '$TMPDIR/annotations.list'
all.exons.output <- '$TMPDIR/${SAMPLE}.all.exons.${CHROMOSOMES}.Rdata'
cnv.calls.file <- '$TMPDIR/${SAMPLE}.cnv.calls.${CHROMOSOMES}.tsv'
summary.file <- '$TMPDIR/${SAMPLE}.cnv.calls.${CHROMOSOMES}.summary.tsv'
bam.suffix <- '$BAM_SUFFIX'

call.cnvs(exon.counts.file = exon.counts.file,
          target.bed = target.bed,
	  test.sample = test.sample,
	  ref.samples.file = ref.samples.file,
          annotations.file = annotations.file,
	  all.exons.output = all.exons.output,
	  cnv.calls.file = cnv.calls.file,
	  summary.file = summary.file,
	  bam.suffix = bam.suffix)

" > $R_SCRIPT

echo "`${NOW}`calling CNVs"

R CMD BATCH --no-save --no-restore $R_SCRIPT ${R_SCRIPT}.log

echo "`${NOW}`done"

cp $TMPDIR/$SAMPLE*Rdata $RESULTS_DIR/
cp $TMPDIR/$SAMPLE*tsv $RESULTS_DIR/

chmod 660 $RESULTS_DIR/*
chmod 660 ${R_SCRIPT}.log
