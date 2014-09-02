#!/bin/bash

## script to run GATK for counting covariates before base quality recalibration

#PBS -l walltime=48:00:00
#PBS -l select=1:ncpus=#cpuThreads:mem=7gb

#PBS -M cgi@imperial.ac.uk
#PBS -m ae
#PBS -j oe

#PBS -q pqcgi

NOW="date +%Y-%m-%d%t%T%t"

# load modules
module load gatk/#gatkVersion
module load samtools/#samtoolsVersion
module load java/#javaVersion

NXTGENUTILS_HOME=/groupvol/cgi/bin/nxtgen-utils-#nextGenUtilsVersion

JAVA_XMX=6G
NCT=#nctThreads
RUN_LOG=#runLog

SCRIPT_CODE="GATKPRHC"


LOG_INFO="`$NOW`INFO $SCRIPT_CODE"
LOG_ERR="`$NOW`ERROR $SCRIPT_CODE"
LOG_WARN="`$NOW`WARN $SCRIPT_CODE"
LOG_DEBUG="`$NOW`DEBUG $SCRIPT_CODE"


# define variables
REFERENCE_FASTA=#referenceFasta
SAMPLE=#sample
ANALYSIS_DIR=#analysisDir
REALIGNED_BAM_FILE=#realignedBamFile
RECALIBRATION_REPORT=#recalibrationReport
RECALIBRATED_OUTPUT_BAM=#recalibratedBamOutput
FRAGMENT_FILE=#fragmentFile
FRAGMENT=#fragmentName
GENOMIC_VCF=#genomicVCF
INCLUDES_UNMAPPED=#includesUnmapped
SUMMARY_SCRIPT_PATH=#summaryScriptPath
DOWNSAMPLING=#downsamplingThreshold

mkdir $TMPDIR/tmp

echo "`${NOW}`INFO $SCRIPT_CODE copying realigned BAM and index files to tmp directory..."
cp $REALIGNED_BAM_FILE $TMPDIR/tmp.bam
cp $REALIGNED_BAM_FILE.bai $TMPDIR/tmp.bam.bai 

echo "`${NOW}`INFO $SCRIPT_CODE copying merged recalibration report to tmp directory..."
cp $RECALIBRATION_REPORT $TMPDIR/merged_recal_data.grp

echo "`${NOW}`INFO $SCRIPT_CODE copying reference to tmp directory..."
cp $REFERENCE_FASTA $TMPDIR/reference.fa
cp $REFERENCE_FASTA.fai $TMPDIR/reference.fa.fai
REFRENCE_SEQ_DICT=`echo $REFERENCE_FASTA | perl -pe 's/\.fa/\.dict/'`
cp $REFRENCE_SEQ_DICT $TMPDIR/reference.dict

cp $FRAGMENT_FILE $TMPDIR/fragment.intervals

INTERVAL_ARG="-L $TMPDIR/fragment.intervals"
if [ "$INCLUDES_UNMAPPED" == "T" ]
then
	INTERVAL_ARG_PR="$INTERVAL_ARG -L unmapped"
else
	INTERVAL_ARG_PR=$INTERVAL_ARG
fi

#############################################
# print reads with recalibrated base calls
############################################
echo "`${NOW}`INFO $SCRIPT_CODE recalibrating chunk BAM..."
java -Xmx$JAVA_XMX -XX:+UseSerialGC -Djava.io.tmpdir=$TMPDIR/tmp -jar $GATK_HOME/GenomeAnalysisTK.jar \
   -nct $NCT \
   -T PrintReads \
   -R $TMPDIR/reference.fa \
   -I $TMPDIR/tmp.bam \
   -BQSR $TMPDIR/merged_recal_data.grp \
   -o $TMPDIR/realigned.recalibrated.bam \
   $INTERVAL_ARG_PR

# samtools index recalibrated BAM file
echo "`${NOW}`INFO $SCRIPT_CODE indexing recalibrated BAM..."
samtools index $TMPDIR/realigned.recalibrated.bam

#check that realigned and recalibrated bams have the same number of reads, 
# if so, copy to destination folder

READ_COUNT_INPUT=`samtools flagstat $TMPDIR/tmp.bam | head -n 1 | perl -e 'while(<>){ if(/(\d*?)\s\+\s(\d*?)\s/) { $retval=$1+$2; print "$retval\n"; }  }'`	
READ_COUNT_OUTPUT=`samtools flagstat $TMPDIR/realigned.recalibrated.bam | head -n 1 | perl -e 'while(<>){ if(/(\d*?)\s\+\s(\d*?)\s/) { $retval=$1+$2; print "$retval\n"; }  }'`
echo "Number of reads in realigned file: $READ_COUNT_INPUT"
echo "Number of reads in recalibrated file: $READ_COUNT_OUTPUT"

if [[ $READ_COUNT_INPUT -eq $READ_COUNT_OUTPUT ]]; then
	echo "`${NOW}`INFO $SCRIPT_CODE copying recalibrated BAM to output directory $ANALYSIS_DIR..."
	cp $TMPDIR/realigned.recalibrated.bam $RECALIBRATED_OUTPUT_BAM
	cp $TMPDIR/realigned.recalibrated.bam.bai $RECALIBRATED_OUTPUT_BAM.bai
else
	echo "ERROR: $SCRIPT_CODE realined and recalibrated bams have different numbers of reads"
fi

#logging and deleting internediate file
if [[ -s $RECALIBRATED_OUTPUT_BAM ]]; then
	STATUS=OK
	echo "deleting realigned bam file $REALIGNED_BAM_FILE"
	rm $REALIGNED_BAM_FILE $REALIGNED_BAM_FILE.bai
else
	STATUS=FAILED
	echo "recalibration failed, keeping realigned file for rerun $REALIGNED_BAM_FILE"
fi

echo -e "`${NOW}`$SCRIPT_CODE\t$SAMPLE\t$FRAGMENT\trecalibrated_bam\t$STATUS" >> $RUN_LOG

#########################################################
# run HaplotypeCaller on realigned recalibrated BAM file
#########################################################
## although we set up downsampling, it seems that HaplotypeCaller ignores it and uses hardcoded threshold of 1000
## http://gatkforums.broadinstitute.org/discussion/3989/downsampling-to-coverage-and-the-3-x-haplotypecaller
## -nct is not used, because it sometimes cause HaplotypeCaller to fail

java -Xmx$JAVA_XMX -XX:+UseSerialGC -Djava.io.tmpdir=$TMPDIR/tmp -jar $GATK_HOME/GenomeAnalysisTK.jar \
	-T HaplotypeCaller \
	-dcov $DOWNSAMPLING \
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

