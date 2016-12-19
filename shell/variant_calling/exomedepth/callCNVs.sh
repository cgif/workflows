#!/bin/bash

#PBS -l walltime=1:00:00
#PBS -l select=1:ncpus=1:mem=5gb

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
RESULTS_DIR=#resultsFolder
TARGET=#target
PED_FILE=#pedFile
ANNOTATIONS=#annotations
BAM_SUFFIX='#bamSuffix'
SUMMARY_SCRIPT_PATH=#summaryScriptPath


SAMPLE=`basename $RESULTS_DIR`
PROJECT_RESULTS_DIR=`dirname $RESULTS_DIR`

# select reference sets
FAMILY=`grep $SAMPLE $PED_FILE|cut -f 1|uniq`
 
#for autosomes, select unrelated samples for reference set
grep -Pv "$FAMILY\t" $PED_FILE | grep -v "#" | cut -f 2 > $TMPDIR/ref.autosomes.sample.list

#for X chromosome, select unrelated same sex samples for reference set
SEX=`grep -P "$FAMILY\t$SAMPLE" $PED_FILE|cut -f 5` 
grep -Pv "$FAMILY\t" $PED_FILE | grep -v "#" | cut -f 2,5 | grep -P "\t$SEX" | cut -f 1 > $TMPDIR/ref.X_chromosome.sample.list

echo "`$NOW`copying exon counts file to $TMPDIR"

cp $PROJECT_RESULTS_DIR/multisample/exon.counts.Rdata $TMPDIR

echo "`$NOW`copying target and annotations list files to $TMPDIR"

cp $TARGET $TMPDIR/target.bed
cp $ANNOTATIONS $TMPDIR/annotations.list

echo "`$NOW`writing R script"

echo "
source('$R_FUNCTIONS')

exon.counts.file <- '$TMPDIR/exon.counts.Rdata'
exon.counts.autosomes.file <- '$TMPDIR/exon.counts.autosomes.tsv'
exon.counts.X.file <- '$TMPDIR/exon.counts.X_chromosome.tsv'
bam.suffix <- '$BAM_SUFFIX'

make.dafr(exon.counts.file = exon.counts.file,
		exon.counts.autosomes.file = exon.counts.autosomes.file,
		exon.counts.X.file = exon.counts.X.file,
		bam.suffix)

target.bed <- '$TMPDIR/target.bed'
test.sample <- '$SAMPLE'
annotations.file <- '$TMPDIR/annotations.list'
prefix <- '$TMPDIR/${SAMPLE}'
cnv.output <- '$TMPDIR/${SAMPLE}.cnv.calls'

#for autosomes:
exon.counts.dafr.file <- '$TMPDIR/exon.counts.autosomes.tsv'
ref.samples.file <- '$TMPDIR/ref.autosomes.sample.list'
analysis.subset <- 'autosomes'

call.cnvs(exon.counts.dafr.file = exon.counts.dafr.file,
		test.sample = test.sample,
		ref.samples.file = ref.samples.file,
		annotations.file = annotations.file,
		analysis.subset = analysis.subset,
		prefix = prefix,
		target.bed = target.bed)

#for X chromosome:
exon.counts.dafr.file <- '$TMPDIR/exon.counts.X_chromosome.tsv'
ref.samples.file <- '$TMPDIR/ref.X_chromosome.sample.list'
analysis.subset <- 'X_chromosome'

call.cnvs(exon.counts.dafr.file = exon.counts.dafr.file,
		test.sample = test.sample,
		ref.samples.file = ref.samples.file,
		annotations.file = annotations.file,
		analysis.subset = analysis.subset,
		prefix = prefix,
		target.bed = target.bed)

sessionInfo()

" > $R_SCRIPT

chmod 770 $R_SCRIPT
echo "`${NOW}`calling CNVs"

R CMD BATCH --no-save --no-restore $R_SCRIPT ${R_SCRIPT}.log

echo "`${NOW}`merging outputs for $SAMPLE"

AUTO=$TMPDIR/${SAMPLE}.cnv.calls.autosomes.tsv
CHRX=$TMPDIR/${SAMPLE}.cnv.calls.X_chromosome.tsv
ALL=$TMPDIR/${SAMPLE}.cnv.calls.all.tsv

cp $AUTO $ALL
tail -n +2 $CHRX >> $ALL

AUTO_SUM=$TMPDIR/${SAMPLE}.cnv.calls.autosomes.summary.tsv
CHRX_SUM=$TMPDIR/${SAMPLE}.cnv.calls.X_chromosome.summary.tsv
ALL_SUM=$TMPDIR/${SAMPLE}.cnv.calls.all.summary.tsv

cp $AUTO_SUM $ALL_SUM
tail -n +2 $CHRX_SUM >> $ALL_SUM

#run summary script
perl $SUMMARY_SCRIPT_PATH

cp $TMPDIR/$SAMPLE*Rdata $RESULTS_DIR/
cp $TMPDIR/$SAMPLE*tsv $RESULTS_DIR/

chmod 660 $RESULTS_DIR/*
chmod 660 ${R_SCRIPT}.log

echo "`${NOW}`done"
