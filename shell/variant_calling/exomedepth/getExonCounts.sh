#!/bin/bash

#PBS -l walltime=72:00:00
#PBS -l select=1:ncpus=1:mem=10gb

#PBS -M cgi@imperial.ac.uk
#PBS -m bea
#PBS -j oe

#PBS -q pqcgi

# load modules

module load R/3.0.1
module load samtools/0.1.18

NOW="date +%Y-%m-%d%t%T%t"

# define variables
R_FUNCTIONS=#Rfunctions
R_SCRIPT=#Rscript
TARGET=#target
BAM_LIST=#BamList
RESULTS_DIR=#resultsFolder
CHROM_X=#chromX

echo "`$NOW`creating R script for counting reads per exon"
echo "`$NOW`R script: $R_SCRIPT"
echo "`$NOW`input: $BAM_LIST"
echo "`$NOW`target: $TARGET"
echo "`$NOW`analise X chromosome with same sex samples: $CHROM_X"
echo "`$NOW`results: $RESULTS_DIR"


if [ "$CHROM_X" == "T" ]
then

    grep -Pv '^[X|Y]\s' $TARGET > $TMPDIR/target.autosomes.bed
    grep -P '^X\s' $TARGET > $TMPDIR/target.chrX.bed

    cp $BAM_LIST $TMPDIR/bam.list
    echo -n "" > $TMPDIR/bam.autosomes.list
    echo -n "" > $TMPDIR/bam.chrX.list

    while read BAM_PATH
    do

        SAMPLE=`basename $BAM_PATH .bam`
	cp $BAM_PATH $TMPDIR/$SAMPLE.bam

	samtools view -b -h -L $TMPDIR/target.autosomes.bed -o $TMPDIR/$SAMPLE.autosomes.bam $TMPDIR/$SAMPLE.bam
	samtools index $TMPDIR/$SAMPLE.autosomes.bam
	echo "$TMPDIR/$SAMPLE.autosomes.bam" >> $TMPDIR/bam.autosomes.list

	samtools view -b -h -L $TMPDIR/target.chrX.bed -o $TMPDIR/$SAMPLE.chrX.bam $TMPDIR/$SAMPLE.bam
	samtools index $TMPDIR/$SAMPLE.chrX.bam
	echo "$TMPDIR/$SAMPLE.chrX.bam" >> $TMPDIR/bam.chrX.list

    done < $TMPDIR/bam.list

    echo "
source('$R_FUNCTIONS')

target.bed <- '$TMPDIR/target.autosomes.bed'
bam.files <- '$TMPDIR/bam.autosomes.list'
exon.counts.file <- '$TMPDIR/exon.counts.autosomes.Rdata'

get.exon.counts(target.bed = target.bed,
                bam.files = bam.files,
		exon.counts.file = exon.counts.file)

target.bed <- '$TMPDIR/target.chrX.bed'
bam.files <- '$TMPDIR/bam.chrX.list'
exon.counts.file <- '$TMPDIR/exon.counts.chrX.Rdata'

get.exon.counts(target.bed = target.bed,
                bam.files = bam.files,
		exon.counts.file = exon.counts.file)

    " > $R_SCRIPT

else

    cp $TARGET $TMPDIR/target.bed
    cp $BAM_LIST $TMPDIR/bam.list

    echo "
source('$R_FUNCTIONS')

target.bed <- '$TMPDIR/target.bed'
bam.files <- '$TMPDIR/bam.list'
exon.counts.file <- '$TMPDIR/exon.counts.Rdata'

get.exon.counts(target.bed = target.bed,
                bam.files = bam.files,
		exon.counts.file = exon.counts.file)

    " > $R_SCRIPT

fi

ls -l

echo "`${NOW}`parsing bam files to get exon counts"
R CMD BATCH --no-save --no-restore $R_SCRIPT ${R_SCRIPT}.log

cp $TMPDIR/*.Rdata $RESULTS_DIR
chmod 660 $RESULTS_DIR/*

echo "`${NOW}`done"

