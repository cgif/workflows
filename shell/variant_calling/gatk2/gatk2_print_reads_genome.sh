#!/bin/bash

## script to run GATK for counting covariates before base quality recalibration

#PBS -l walltime=48:00:00
#PBS -l select=1:ncpus=#cpuThreads:mem=50gb

#PBS -M igf@imperial.ac.uk
#PBS -m ea
#PBS -j oe

#PBS -q pqcgi

NOW="date +%Y-%m-%d%t%T%t"

# load modules
module load gatk/#gatkVersion
module load samtools/#samtoolsVersion
module load java/#javaVersion

NXTGENUTILS_HOME=/groupvol/cgi/bin/nxtgen-utils-#nextGenUtilsVersion

JAVA_XMX=49G
NCT=#cpuThreads
RUN_LOG=#runLog

SCRIPT_CODE="GATKPRRE"


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
REDUCED_OUTPUT_BAM=#reducedBamOutput
DOWNSAMPLING=#downsamplingThreshold
FRAGMENT_FILE=#fragmentFile
FRAGMENT=#fragmentName
INCLUDES_UNMAPPED=#includesUnmapped
SUMMARY_SCRIPT_PATH=#summaryScriptPath

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
	INTERVAL_ARG="$INTERVAL_ARG -L unmapped"
fi

# step 5: print reads with recalibrated base calls
echo "`${NOW}`INFO $SCRIPT_CODE recalibrating chunk BAM..."
java -Xmx$JAVA_XMX -XX:+UseSerialGC -Djava.io.tmpdir=$TMPDIR/tmp -jar $GATK_HOME/GenomeAnalysisTK.jar \
   -nct $NCT \
   -T PrintReads \
   -R $TMPDIR/reference.fa \
   -I $TMPDIR/tmp.bam \
   -BQSR $TMPDIR/merged_recal_data.grp \
   -o $TMPDIR/realigned.recalibrated.bam \
   $INTERVAL_ARG

# step 6: samtools index recalibrated BAM file
echo "`${NOW}`INFO $SCRIPT_CODE indexing recalibrated BAM..."
samtools index $TMPDIR/realigned.recalibrated.bam

echo "`${NOW}`INFO $SCRIPT_CODE copying recalibrated BAM to output directory $ANALYSIS_DIR..."
cp $TMPDIR/realigned.recalibrated.bam $RECALIBRATED_OUTPUT_BAM
cp $TMPDIR/realigned.recalibrated.bam.bai $RECALIBRATED_OUTPUT_BAM.bai

#logging
STATUS=OK
if [[ ! -e $RECALIBRATED_OUTPUT_BAM ]]
then
	STATUS=FAILED
fi

echo -e "`${NOW}`$SCRIPT_CODE\t$SAMPLE\t$FRAGMENT\trecalibrated_bam\t$STATUS" >> $RUN_LOG


# reduce reads
echo "`${NOW}`INFO $SCRIPT_CODE reducing recalibrated BAM..."

java -Xmx$JAVA_XMX -XX:+UseSerialGC -Djava.io.tmpdir=$TMPDIR/tmp -jar $GATK_HOME/GenomeAnalysisTK.jar \
   -R $TMPDIR/reference.fa \
   -T ReduceReads \
   -I $TMPDIR/realigned.recalibrated.bam \
   -o $TMPDIR/realigned.recalibrated.reduced.bam \
   -ds $DOWNSAMPLING \
   -dcov $DOWNSAMPLING

echo "`${NOW}`INFO $SCRIPT_CODE indexing reduced BAM"
samtools index $TMPDIR/realigned.recalibrated.reduced.bam

echo "`${NOW}`INFO $SCRIPT_CODE copying recalibrated BAM to output directory $ANALYSIS_DIR..."
cp $TMPDIR/realigned.recalibrated.reduced.bam $REDUCED_OUTPUT_BAM
cp $TMPDIR/realigned.recalibrated.reduced.bam.bai $REDUCED_OUTPUT_BAM.bai


#logging
STATUS=OK
if [[ ! -e $REDUCED_OUTPUT_BAM ]]
then
	STATUS=FAILED
fi

echo -e "`${NOW}`$SCRIPT_CODE\t$SAMPLE\t$FRAGMENT\treduced_bam\t$STATUS" >> $RUN_LOG

#run summary script
perl $SUMMARY_SCRIPT_PATH

echo "`${NOW}`INFO $SCRIPT_CODE done"

