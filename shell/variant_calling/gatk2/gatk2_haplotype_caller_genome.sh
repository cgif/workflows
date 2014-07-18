#!/bin/bash

## script to run GATK Haplotype Caller for genome

#PBS -l walltime=72:00:00
#PBS -l select=1:ncpus=1:mem=15800mb

#PBS -M cgi@imperial.ac.uk
#PBS -m ea
#PBS -j oe

#PBS -q pqcgi

NOW="date +%Y-%m-%d%t%T%t"


# load modules
module load gatk/#gatkVersion
module load samtools/#samtoolsVersion
module load picard/#picardVersion
module load R/#rVersion
module load java/#javaVersion

JAVA_XMX=15G

SCRIPT_CODE="GATKHACA"



LOG_INFO="`${NOW}`INFO $SCRIPT_CODE"
LOG_ERR="`${NOW}`ERR $SCRIPT_CODE"
LOG_WARN="`${NOW}`WARN $SCRIPT_CODE"
LOG_DEBUG="`${NOW}`DEBUG $SCRIPT_CODE"

BUNDLE=/groupvol/cgi/resources/GATK_resource_bundle/2.3/b37

# define variables

INPUT_BAM="#inputBam"
DBSNP=#dbSnp
REFERENCE_FASTA=#referenceFasta
REFRENCE_SEQ_DICT=`echo $REFERENCE_FASTA | perl -pe 's/\.fa/\.dict/'`
ANALYSIS_DIR=#analysisDir
SAMPLE=#sampleName
FRAGMENT=#fragmentName
PED_FILE=#pedFile
FRAGMENT_FILE=#fragmentFile
DOWNSAMPLING=#downsamplingThreshold
TYPE=#seqType
RUN_LOG=#runLog

#copy input files to tmp dir

echo "`${NOW}`INFO $SCRIPT_CODE copying GATK resources to tmp directory..."
DBSNP_FILENAME=`basename $DBSNP`

echo "`${NOW}`INFO $SCRIPT_CODE $DBSNP"
cp $DBSNP $TMPDIR/$DBSNP_FILENAME
cp $DBSNP.idx $TMPDIR/$DBSNP_FILENAME.idx


cp $FRAGMENT_FILE $TMPDIR/fragment.intervals

echo "`${NOW}`INFO $SCRIPT_CODE copying reference fasta and indexto tmp directory..."
cp $REFERENCE_FASTA $TMPDIR/reference.fa
cp $REFERENCE_FASTA.fai $TMPDIR/reference.fa.fai
cp $REFRENCE_SEQ_DICT $TMPDIR/reference.dict

# make tmp folder for temporary java files
mkdir $TMPDIR/tmp

PED_ARG=""

##### check for ped file
if [[ "$PED_FILE" != "none" ]]; then
	PED_ARG="-ped $PED_FILE -pedValidationType SILENT"
fi

echo "`${NOW}`INFO $SCRIPT_CODE copying chunk BAM and index file to tmp directory..."

#get analysis directory of GATK2 run
ANALYSIS_DIR_RUN=`dirname $ANALYSIS_DIR`

INPUT_BAM_ARGUMENT=""
for SAMPLE_NAME in `ls --color=never $ANALYSIS_DIR_RUN | grep -v multisample`
do

	BAM="$ANALYSIS_DIR_RUN/$SAMPLE_NAME/recalibration/$SAMPLE_NAME.$FRAGMENT.realigned.recalibrated.reduced.bam"
	BAM_NAME=`basename $BAM`
	
	echo "`${NOW}`INFO $SCRIPT_CODE $BAM_NAME"
	cp $BAM $TMPDIR
	cp $BAM.bai $TMPDIR

	INPUT_BAM_ARGUMENT="$INPUT_BAM_ARGUMENT	-I $TMPDIR/$BAM_NAME"

done;

CUT_OFF=30
if [[ "$TYPE" != "TARGETED" ]]; then
	CUT_OFF=20
fi


# step 13b: haplotype caller
echo  "`${NOW}`INFO $SCRIPT_CODE running HaplotypeCaller"

java -Xmx$JAVA_XMX -Djava.io.tmpdir=$TMPDIR/tmp -jar $GATK_HOME/GenomeAnalysisTK.jar \
    -T HaplotypeCaller \
    -R $TMPDIR/reference.fa \
    $INPUT_BAM_ARGUMENT \
    --dbsnp $DBSNP_FILENAME \
    -G Standard \
    -A AlleleBalance \
    -stand_call_conf $CUT_OFF \
    -stand_emit_conf 10 \
    -dcov $DOWNSAMPLING \
    --bamOutput $SAMPLE.$FRAGMENT.haplotype.bam \
    -o $TMPDIR/$SAMPLE.$FRAGMENT.raw.vcf \
    -L $TMPDIR/fragment.intervals \
    $PED_ARG

echo "`${NOW}`INFO $SCRIPT_CODE copying raw VCF file to output directory $ANALYSIS_DIR/haplotypecaller..."
cp $SAMPLE.$FRAGMENT.raw.vcf $ANALYSIS_DIR/haplotypecaller/

echo "`${NOW}`INFO $SCRIPT_CODE copying assembled haplotype BAM to output directory $ANALYSIS_DIR/haplotypecaller..."
cp $SAMPLE.$FRAGMENT.haplotype.bam $ANALYSIS_DIR/haplotypecaller/

echo "`${NOW}`INFO $SCRIPT_CODE done"

#logging
STATUS=OK
if [[ ! -e $ANALYSIS_DIR/haplotypecaller/$SAMPLE.$FRAGMENT.raw.vcf ]]
then
	STATUS=FAILED
fi

echo -e "`${NOW}`$SCRIPT_CODE\tmultisample\t$FRAGMENT\traw_vcf\t$STATUS" >> $RUN_LOG

