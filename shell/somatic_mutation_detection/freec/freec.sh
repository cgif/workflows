#!/bin/bash

## script to run FREEC

#PBS -l walltime=72:00:00
#PBS -l select=1:ncpus=6:mem=10gb:tmpspace=#tmpSpacegb

#PBS -m ea
#PBS -M igf@imperial.ac.uk
#PBS -j oe

#PBS -q pqcgi

module load freec/#freecVersion
module load samtools/#SamtoolsVersion
module load R/#Rversion

SEX=#Sex
CHROM_LENGTH=#ChrLenFile
SAMTOOLS_VERSION=#SamtoolsVersion

FREEC_PROFILE=#freecProfile
FREEC_GRAPH=#freecGraph
FREEC_SIGN=#freecSign

RESULTS_DIR=#ResultsDir
NORMAL_BAM=#NormalBam
TUMOR_BAM=#TumorBam

WINDOW=#Window
STEP=#Step
PLOIDY=#Ploidy
BREAK_POINT=#BreakPoint

REF_FASTA=#RefFasta
TMP_REF_FASTA=$TMPDIR/ref.fasta
cp $REF_FASTA $TMP_REF_FASTA
cp $REF_FASTA.fai $TMP_REF_FASTA.fai

TMP_NORMAL_BAM=$TMPDIR/normal.bam
cp $NORMAL_BAM $TMP_NORMAL_BAM

#filter reads with non-uniq mapping
samtools view -h $TMPDIR/normal.bam|awk '!($5 == 0) {print}'|samtools view -hSb -  > $TMPDIR/normal.filtered.bam
mv $TMPDIR/normal.filtered.bam $TMPDIR/normal.bam 
ls -l 
samtools index $TMPDIR/normal.bam

TMP_TUMOR_BAM=$TMPDIR/tumor.bam
cp $TUMOR_BAM $TMP_TUMOR_BAM

#filter reads with non-uniq mapping
samtools view -h $TMPDIR/tumor.bam|awk '!($5 == 0) {print}'|samtools view -hSb - > $TMPDIR/tumor.filtered.bam
mv $TMPDIR/tumor.filtered.bam $TMPDIR/tumor.bam
ls -l 
samtools index $TMPDIR/tumor.bam

TMP_CHROM_LENGTH=$TMPDIR/chrom.len
cp $CHROM_LENGTH $TMP_CHROM_LENGTH

TMP_FREEC_PROFILE=$TMPDIR/FREEC.profile
cp $FREEC_PROFILE $TMP_FREEC_PROFILE

sed -i -e "s/#resultsDir/${TMPDIR//\//\\/}/" $TMP_FREEC_PROFILE
sed -i -e "s/#normalBam/${TMP_NORMAL_BAM//\//\\/}/" $TMP_FREEC_PROFILE
sed -i -e "s/#tumorBam/${TMP_TUMOR_BAM//\//\\/}/" $TMP_FREEC_PROFILE
sed -i -e "s/#sex/${SEX}/" $TMP_FREEC_PROFILE
sed -i -e "s/#chrLenFile/${TMP_CHROM_LENGTH//\//\\/}/" $TMP_FREEC_PROFILE
sed -i -e "s/#samtoolsVersion/${SAMTOOLS_VERSION//\//\\/}/" $TMP_FREEC_PROFILE
sed -i -e "s/#window/${WINDOW}/" $TMP_FREEC_PROFILE
sed -i -e "s/#step/${STEP}/" $TMP_FREEC_PROFILE
sed -i -e "s/#polidy/${PLOIDY}/" $TMP_FREEC_PROFILE
sed -i -e "s/#breakPoint/${BREAK_POINT}/" $TMP_FREEC_PROFILE

more $TMP_FREEC_PROFILE

freec -conf $TMP_FREEC_PROFILE

cat $FREEC_SIGN | R --slave --args $TMPDIR/tumor.bam_CNVs $TMPDIR/tumor.bam_ratio.txt

cat $FREEC_GRAPH | R --slave --args 2 $TMPDIR/tumor.bam_ratio.txt

SAMPLE=`basename $RESULTS_DIR`
cp $TMPDIR/tumor.bam_CNVs.p.value.txt $RESULTS_DIR/${SAMPLE}.CNVs.txt
cp $TMPDIR/tumor.bam_ratio.txt.png $RESULTS_DIR/${SAMPLE}.png

