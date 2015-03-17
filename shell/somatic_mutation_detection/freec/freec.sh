#!/bin/bash

## script to run FREEC

#PBS -l walltime=72:00:00
#PBS -l select=1:ncpus=6:mem=10gb:tmpspace=#tmpSpacegb

#PBS -m ea
#PBS -M cgi@imperial.ac.uk
#PBS -j oe

#PBS -q pqcgi

module load freec/#freecVersion
module load samtools/#SamtoolsVersion
module load R/#Rversion

CHROM_LENGTH=#ChrLenFile
FREEC_PROFILE=#freecProfile
SAMTOOLS_VERSION=#SamtoolsVersion
FREEC_GRAPH=#freecGraph
FREEC_SIGN=#freecSign

RESULTS_DIR=#ResultsDir
NORMAL_BAM=#NormalBam
TUMOR_BAM=#TumorBam

TMP_NORMAL_BAM=$TMPDIR/normal.bam
cp $NORMAL_BAM $TMP_NORMAL_BAM

#filter reads with non-uniq mapping
samtools view -h $TMPDIR/normal.bam|awk '!($5 == 0) {print}'|samtools view -hSb -  > $TMPDIR/normal.filtered.bam
ls -lh
mv $TMPDIR/normal.filtered.bam $TMPDIR/normal.bam
samtools index $TMPDIR/normal.bam

TMP_TUMOR_BAM=$TMPDIR/tumor.bam
cp $TUMOR_BAM $TMP_TUMOR_BAM

#filter reads with non-uniq mapping
samtools view -h $TMPDIR/tumor.bam|awk '!($5 == 0) {print}'|samtools view -hSb - > $TMPDIR/tumor.filtered.bam
ls -lh
mv $TMPDIR/tumor.filtered.bam $TMPDIR/tumor.bam
samtools index $TMPDIR/tumor.bam

TMP_CHROM_LENGTH=$TMPDIR/chrom.len
cp $CHROM_LENGTH $TMP_CHROM_LENGTH

TMP_FREEC_PROFILE=$TMPDIR/FREEC.profile
cp $FREEC_PROFILE $TMP_FREEC_PROFILE

sed -i -e "s/#resultsDir/${TMPDIR//\//\\/}/" $TMP_FREEC_PROFILE
sed -i -e "s/#normalBam/${TMP_NORMAL_BAM//\//\\/}/" $TMP_FREEC_PROFILE
sed -i -e "s/#tumorBam/${TMP_TUMOR_BAM//\//\\/}/" $TMP_FREEC_PROFILE

sed -i -e "s/#chrLenFile/${TMP_CHROM_LENGTH//\//\\/}/" $TMP_FREEC_PROFILE
sed -i -e "s/#samtoolsVersion/${SAMTOOLS_VERSION//\//\\/}/" $TMP_FREEC_PROFILE

freec -conf $TMP_FREEC_PROFILE

cat $FREEC_SIGN | R --slave --args $TMPDIR/tumor.bam_CNVs $TMPDIR/tumor.bam_ratio.txt

cat $FREEC_GRAPH | R --slave --args 2 $TMPDIR/tumor.bam_ratio.txt

cp $TMPDIR/tumor.bam_CNVs* $RESULTS_DIR
cp $TMPDIR/tumor.bam_ratio* $RESULTS_DIR

ls -lh
