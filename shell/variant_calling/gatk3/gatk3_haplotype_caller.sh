#!/bin/bash

## script to run GATK for counting covariates before base quality recalibration

#PBS -l walltime=48:00:00
#PBS -l select=1:ncpus=1:mem=7gb

#PBS -M cgi@imperial.ac.uk
#PBS -m ae
#PBS -j oe

#PBS -q pqcgi

NOW="date +%Y-%m-%d%t%T%t"

# load modules
module load gatk/#gatkVersion
module load samtools/#samtoolsVersion
module load java/#javaVersion

JAVA_XMX=6G
RUN_LOG=#runLog

SCRIPT_CODE="GATKHAPC"


LOG_INFO="`$NOW`INFO $SCRIPT_CODE"
LOG_ERR="`$NOW`ERROR $SCRIPT_CODE"
LOG_WARN="`$NOW`WARN $SCRIPT_CODE"
LOG_DEBUG="`$NOW`DEBUG $SCRIPT_CODE"


# define variables
REFERENCE_FASTA=#referenceFasta
SAMPLE=#sample
ANALYSIS_DIR=#analysisDir
RECALIBRATED_BAM=#recalibratedBam
FRAGMENT_FILE=#fragmentFile
FRAGMENT=#fragmentName
GENOMIC_VCF=#genomicVCF
SUMMARY_SCRIPT_PATH=#summaryScriptPath
DOWNSAMPLING=#downsamplingThreshold

mkdir $TMPDIR/tmp

echo "`${NOW}`INFO $SCRIPT_CODE copying recalibrated BAM and index files to tmp directory..."
cp $RECALIBRATED_BAM $TMPDIR/realigned.recalibrated.bam
cp $RECALIBRATED_BAM.bai $TMPDIR/realigned.recalibrated.bam.bai 

echo "`${NOW}`INFO $SCRIPT_CODE copying reference to tmp directory..."
cp $REFERENCE_FASTA $TMPDIR/reference.fa
cp $REFERENCE_FASTA.fai $TMPDIR/reference.fa.fai
REFRENCE_SEQ_DICT=`echo $REFERENCE_FASTA | perl -pe 's/\.fa/\.dict/'`
cp $REFRENCE_SEQ_DICT $TMPDIR/reference.dict

cp $FRAGMENT_FILE $TMPDIR/fragment.intervals

INTERVAL_ARG="-L $TMPDIR/fragment.intervals"


#########################################################
# run HaplotypeCaller on realigned recalibrated BAM file
#########################################################
## although we set up downsampling, it seems that HaplotypeCaller ignores it and uses hardcoded threshold of 1000
## http://gatkforums.broadinstitute.org/discussion/3989/downsampling-to-coverage-and-the-3-x-haplotypecaller
## -nct is not used, because it sometimes cause HaplotypeCaller to fail

java -Xmx$JAVA_XMX -XX:+UseSerialGC -Djava.io.tmpdir=$TMPDIR/tmp -jar $GATK_HOME/GenomeAnalysisTK.jar \
	-T HaplotypeCaller \
	-R $TMPDIR/reference.fa \
	-I $TMPDIR/realigned.recalibrated.bam \
	-ERC GVCF \
	--variant_index_type LINEAR \
	--variant_index_parameter 128000 \
	-o $TMPDIR/HCgenomic.vcf \
	-rf BadCigar \
	$INTERVAL_ARG

echo "`${NOW}`INFO $SCRIPT_CODE copying gVCF to output directory $ANALYSIS_DIR..."
cp $TMPDIR/HCgenomic.vcf $GENOMIC_VCF

#logging
STATUS=OK
if [[ ! -s $GENOMIC_VCF ]]
then
	STATUS=FAILED
fi

echo -e "`${NOW}`$SCRIPT_CODE\t$SAMPLE\t$FRAGMENT\tgenomic_vcf\t$STATUS" >> $RUN_LOG

#run summary script
perl $SUMMARY_SCRIPT_PATH

echo "`${NOW}`INFO $SCRIPT_CODE done"

