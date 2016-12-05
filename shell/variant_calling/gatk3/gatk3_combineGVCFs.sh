#!/bin/bash

## script to combine genomic vcf files from multiple samples

#PBS -l walltime=72:00:00
#PBS -l select=1:ncpus=1:mem=16gb:tmpspace=#tmpSpacegb

#PBS -M igf@imperial.ac.uk
#PBS -m ea
#PBS -j oe

#PBS -q pqcgi

# load modules

module load gatk/#gatkVersion
module load java/#javaVersion

NOW="date +%Y-%m-%d%t%T%t"
JAVA_XMX=15G

# define variables

REFERENCE_FASTA=#referenceFasta
ANALYSIS_DIR=#analysisDir
PREFIX=#prefixName		### prefix for output file
FRAGMENT=#fragmentName

echo "`${NOW}`INFO copying reference fasta and index to tmp directory..."
cp $REFERENCE_FASTA $TMPDIR/reference.fa
cp $REFERENCE_FASTA.fai $TMPDIR/reference.fa.fai
REFERENCE_SEQ_DICT=`echo $REFERENCE_FASTA | perl -pe 's/\.fa/\.dict/'`
cp $REFERENCE_SEQ_DICT $TMPDIR/reference.dict
echo "`${NOW}`INFO copying chunk gVCF files to tmp directory..."

#get analysis directory of GATK3 run
ANALYSIS_DIR_RUN=`dirname $ANALYSIS_DIR`

#get sample names from analysis directory and copy
#respective genomic VCF files
INPUT_GVCF_ARGUMENT=""

for SAMPLE_NAME in `ls --color=never $ANALYSIS_DIR_RUN | grep -vE "multisample|aux"`; do

	GVCF="$ANALYSIS_DIR_RUN/$SAMPLE_NAME/chunks/$SAMPLE_NAME.$FRAGMENT.genomic.vcf"
	GVCF_NAME=`basename $GVCF`
	
	echo "`${NOW}`INFO $SCRIPT_CODE copying $GVCF to temporary space"
	cp $GVCF $TMPDIR

### unsure whether to copy indices
#	echo "`${NOW}`INFO $SCRIPT_CODE copying $GVCF.idx to temporary space"
#	cp $GVCF.idx $TMPDIR
	
	INPUT_GVCF_ARGUMENT="$INPUT_GVCF_ARGUMENT	-V $TMPDIR/$GVCF_NAME"

done;

echo -e "`${NOW}` size of $TMPDIR"
du -h $TMPDIR

# make tmp folder for temporary java files
mkdir $TMPDIR/tmp

echo "`${NOW}`INFO running CombineGVCFs"
java -Xmx$JAVA_XMX -XX:+UseSerialGC -Djava.io.tmpdir=$TMPDIR/tmp -jar $GATK_HOME/GenomeAnalysisTK.jar \
	-R $TMPDIR/reference.fa \
	-T CombineGVCFs \
	$INPUT_GVCF_ARGUMENT \
	-o $TMPDIR/$PREFIX.$FRAGMENT.genomic.vcf \

echo "`${NOW}`INFO copying combined gVCF file to output directory $ANALYSIS_DIR..."

cp $PREFIX.$FRAGMENT.genomic.vcf $ANALYSIS_DIR
cp $PREFIX.$FRAGMENT.genomic.vcf.idx $ANALYSIS_DIR

ls -al $TMPDIR

echo -e "`${NOW}` size of $TMPDIR"
du -h $TMPDIR

echo "`${NOW}`done"


