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
CHROM_X=#chromX
ANNOTATIONS=#annotations

SAMPLE=`basename $RESULTS_DIR`
PROJECT_RESULTS_DIR=`dirname $RESULTS_DIR`

if [[ "$CHROM_X" == "T" ]] 
then

    cut -f 2 $PED_FILE | grep -v "#" > $TMPDIR/all.sample.list

    #select unrelated proband samples for reference set
    FAMILY=`grep $SAMPLE $PED_FILE|cut -f 1|uniq` 
    cat $PED_FILE | awk '$7 == "1"' | grep -Pv "$FAMILY\t" | cut -f 2 > $TMPDIR/ref.autosomes.sample.list

    #chromosome X targets will be analysed against same sex reference set
    SEX=`grep -P "$FAMILY\t$SAMPLE" $PED_FILE|cut -f 5` 
    cat $PED_FILE | awk '$7 == "1"' | grep -Pv "$FAMILY\t" | cut -f 2,5 | grep -P "\t$SEX" | cut -f 1 > $TMPDIR/ref.chrX.sample.list

    cp $PROJECT_RESULTS_DIR/multisample/exon.counts.autosomes.Rdata $TMPDIR/exon.counts.autosomes.Rdata
    cp $PROJECT_RESULTS_DIR/multisample/exon.counts.chrX.Rdata $TMPDIR/exon.counts.chrX.Rdata

    grep -Pv '^[X|Y]\s' $TARGET > $TMPDIR/target.autosomes.bed
    grep -P '^X\s' $TARGET > $TMPDIR/target.chrX.bed

    cp $ANNOTATIONS $TMPDIR/annotations.list

    echo "
source('$R_FUNCTIONS')

exon.counts.file <- '$TMPDIR/exon.counts.autosomes.Rdata'
target.bed <- '$TMPDIR/target.autosomes.bed'
test.sample <- '$SAMPLE'
ref.samples.file <- '$TMPDIR/ref.autosomes.sample.list'
all.samples.file <- '$TMPDIR/all.sample.list'
annotations.file <- '$TMPDIR/annotations.list'
all.exons.output <- '$TMPDIR/${SAMPLE}.all.exons.autosomes.Rdata'
cnv.calls.file <- '$TMPDIR/${SAMPLE}.cnv.calls.autosomes.Rdata'

call.cnvs(exon.counts.file = exon.counts.file,
          target.bed = target.bed,
	  test.sample = test.sample,
	  ref.samples.file = ref.samples.file,
          all.samples.file = all.samples.file,
          annotations.file = annotations.file,
	  all.exons.output = all.exons.output,
	  cnv.calls.file = cnv.calls.file)


exon.counts.file <- '$TMPDIR/exon.counts.chrX.Rdata'
target.bed <- '$TMPDIR/target.chrX.bed'
ref.samples.file <- '$TMPDIR/ref.chrX.sample.list'
all.exons.output <- '$TMPDIR/${SAMPLE}.all.exons.chrX.Rdata'
cnv.calls.file <- '$TMPDIR/${SAMPLE}.cnv.calls.chrX.Rdata'

call.cnvs(exon.counts.file = exon.counts.file,
          target.bed = target.bed,
	  test.sample = test.sample,
	  ref.samples.file = ref.samples.file,
          all.samples.file = all.samples.file,
          annotations.file = annotations.file,
	  all.exons.output = all.exons.output,
	  cnv.calls.file = cnv.calls.file)

    " > $R_SCRIPT

else

    cut -f 2 $PED_FILE | grep -v "#" > $TMPDIR/all.sample.list

    #select unrelated AND unaffected proband samples for reference set
    FAMILY=`grep $SAMPLE $PED_FILE|cut -f 1|uniq` 
    cat $PED_FILE | awk '$7 == "1"' | grep -Pv "$FAMILY\t" | cut -f 2 > $TMPDIR/ref.sample.list

    cp $PROJECT_RESULTS_DIR/multisample/exon.counts.Rdata $TMPDIR/exon.counts.Rdata

    cp $TARGET $TMPDIR/target.bed

    cp $ANNOTATIONS $TMPDIR/annotations.list

    echo "
source('$R_FUNCTIONS')

exon.counts.file <- '$TMPDIR/exon.counts.Rdata'
target.bed <- '$TMPDIR/target.bed'
test.sample <- '$SAMPLE'
ref.samples.file <- '$TMPDIR/ref.sample.list'
all.samples.file <- '$TMPDIR/all.sample.list'
annotations.file <- '$TMPDIR/annotations.list'
all.exons.output <- '$TMPDIR/${SAMPLE}.all.exons.Rdata'
cnv.calls.file <- '$TMPDIR/${SAMPLE}.cnv.calls.Rdata'

call.cnvs(exon.counts.file = exon.counts.file,
          target.bed = target.bed,
	  test.sample = test.sample,
	  ref.samples.file = ref.samples.file,
          all.samples.file = all.samples.file,
          annotations.file = annotations.file,
	  all.exons.output = all.exons.output,
	  cnv.calls.file = cnv.calls.file)
    " > $R_SCRIPT

fi

echo "`${NOW}`calling CNVs"

R CMD BATCH --no-save --no-restore $R_SCRIPT ${R_SCRIPT}.log

echo "`${NOW}`done"

cp $TMPDIR/$SAMPLE*Rdata $RESULTS_DIR/
chmod 660 $RESULTS_DIR/*
