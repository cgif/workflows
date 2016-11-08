#!/bin/bash

## script to run mutect

#PBS -l walltime=24:00:00
#PBS -l select=1:ncpus=1:mem=5gb:tmpspace=#tmpSpaceMbmb

#PBS -M igf@imperial.ac.uk
#PBS -m ea
#PBS -j oe

#PBS -q pqcgi

# load modules
module load mutect/#mutectVersion
module load java/#javaVersion

NOW="date +%Y-%m-%d%t%T%t"
JAVA_XMX=4G

# define variables
ANALYSIS_FILE=#analysisFile
REFERENCE_FASTA=#referenceFasta
REFERENCE_SEQ_DICT=#referenceSeqDict
MUTECT_COSMIC=#mutectCosmic
MUTECT_DBSNP=#mutectDBsnp
NORMAL_BAM=#normalBam
TUMOR_BAM=#tumorBam
INTERVALS_FILE=#intervalsFile
SUMMARY_SCRIPT_PATH=#summaryScriptPath

#copy input files to tmp dir
echo "`${NOW}` copying files to tmp directory..."
cp $REFERENCE_FASTA $TMPDIR/ref.fasta
cp $REFERENCE_FASTA.fai $TMPDIR/ref.fasta.fai
cp $REFERENCE_SEQ_DICT $TMPDIR/ref.dict
cp $MUTECT_COSMIC $TMPDIR/cosmic.vcf
cp $MUTECT_DBSNP $TMPDIR/dbsnp.vcf
cp $NORMAL_BAM $TMPDIR/normal.bam
cp $NORMAL_BAM.bai $TMPDIR/normal.bam.bai
cp $TUMOR_BAM $TMPDIR/tumor.bam
cp $TUMOR_BAM.bai $TMPDIR/tumor.bam.bai
cp $INTERVALS_FILE $TMPDIR/intervals.intervals

# make tmp folder for temporary java files
mkdir $TMPDIR/tmp

echo "`${NOW}` running mutect..."
java -Xmx$JAVA_XMX -XX:+UseSerialGC -jar -Djava.io.tmpdir=$TMPDIR/tmp $MUTECT_HOME/muTect-1.1.4.jar \
--analysis_type MuTect \
--reference_sequence $TMPDIR/ref.fasta \
--cosmic $TMPDIR/cosmic.vcf \
--dbsnp $TMPDIR/dbsnp.vcf \
--intervals $TMPDIR/intervals.intervals \
--input_file:normal $TMPDIR/normal.bam \
--input_file:tumor $TMPDIR/tumor.bam \
--out $TMPDIR/tmp.txt \
--vcf $TMPDIR/tmp.vcf \
-rf BadCigar

#make tsv file for input into Oncotator web server (requested by customers)

cat tmp.vcf | awk -F $'\t' 'BEGIN {OFS=FS} {print $1,$2,$2,$4,$5}' | sed 's/,\([ACGT]\)\t/\t\1/g' | grep -v "#" | awk 'NR==1{print "Chromosome\tStart_Position\tEnd_Position\tReference_Allele\tTumor_Seq_Allele}1' > tmp.tsv



echo "`${NOW}` copying files from tmp directory..."
cp $TMPDIR/tmp.txt $ANALYSIS_FILE.stats
cp $TMPDIR/tmp.vcf $ANALYSIS_FILE.vcf
cp $TMPDIR/tmp.tsv $ANALYSIS_FILE.tsv

perl $SUMMARY_SCRIPT_PATH
