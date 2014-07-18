#!/bin/bash

## script to run GATK

#PBS -l walltime=24:00:00
#PBS -l select=1:ncpus=2:mem=5gb

#PBS -M cgi@imperial.ac.uk
#PBS -m ea
#PBS -j oe

#PBS -q pqcgi

# load modules
module load gatk/#gatkVersion
module load samtools/#samtoolsVersion
module load picard/#picardVersion
module load R/#rVersion
module load java/#javaVersion

NOW="date +%Y-%m-%d%t%T%t"
JAVA_XMX=3800M

SCRIPT_CODE="GATKPRME"

RUN_LOG=#runLog

LOG_INFO="`$NOW`INFO $SCRIPT_CODE"
LOG_ERR="`$NOW`ERR $SCRIPT_CODE"
LOG_WARN="`$NOW`WARN $SCRIPT_CODE"
LOG_DEBUG="`$NOW`DEBUG $SCRIPT_CODE"

#set R library path
#(for GATK gsalib library required to 
#generate recalibration plots)
R_LIBS_USER=/groupvol/cgi/libs/R/2.15
export R_LIBS=${R_LIBS_USER}:${R_LIBS} 

# define variables

INPUT_BAM=#inputBam
INPUT_BAM_NAME=`basename $INPUT_BAM .bam`
REFERENCE_FASTA=#referenceFasta
FRAGMENT_FILE=#fragmentFile
FRAGMENT=#fragmentName
SAMPLE=#sample
ANALYSIS_DIR=#analysisDir
RESULTS_DIR=#resultsDir
PRE_RECALIBRATION_REPORT=#preRecalibrationReport

#GATK resources
INDELS_1000G=#indels1000G
INDELS_GOLDSTD=#indelsGoldStd
DBSNP=#dbSnp

echo "`${NOW}`copying GATK resources to tmp directory..."
INDELS_1000G_FILENAME=`basename $INDELS_1000G`
INDELS_GOLDSTD_FILENAME=`basename $INDELS_GOLDSTD`
DBSNP_FILENAME=`basename $DBSNP`

echo "`${NOW}`$INDELS_1000G"
cp $INDELS_1000G $TMPDIR/$INDELS_1000G_FILENAME
cp $INDELS_1000G.idx $TMPDIR/$INDELS_1000G_FILENAME.idx

echo "`${NOW}`$INDELS_GOLDSTD"
cp $INDELS_GOLDSTD $TMPDIR/$INDELS_GOLDSTD_FILENAME
cp $INDELS_GOLDSTD.idx $TMPDIR/$INDELS_GOLDSTD_FILENAME.idx

echo "`${NOW}`$DBSNP"
cp $DBSNP $TMPDIR/$DBSNP_FILENAME
cp $DBSNP.idx $TMPDIR/$DBSNP_FILENAME.idx

echo "`${NOW}`copying chunk BAM and index file to tmp directory..."
cp $INPUT_BAM $TMPDIR/chunk.bam
cp $INPUT_BAM.bai $TMPDIR/chunk.bam.bai

echo "`${NOW}`copying intervals file to tmp directory..."
cp $FRAGMENT_FILE $TMPDIR/fragment.intervals

echo "`${NOW}`copying reference fasta and indexto tmp directory..."
cp $REFERENCE_FASTA $TMPDIR/reference.fa
cp $REFERENCE_FASTA.fai $TMPDIR/reference.fa.fai
REFRENCE_SEQ_DICT=`echo $REFERENCE_FASTA | perl -pe 's/\.fa/\.dict/'`
cp $REFRENCE_SEQ_DICT $TMPDIR/reference.dict

# step 4: run GATK BaseRecalibrator to generate recalibration plots
# see http://gatkforums.broadinstitute.org/discussion/1539/baserecalibrator-plots
echo "`${NOW}`generating recalibration report and plot for realigned and recalibrated chunk BAM..."
java -Xmx$JAVA_XMX -jar $GATK_HOME/GenomeAnalysisTK.jar \
   -T BaseRecalibrator \
   -I $TMPDIR/chunk.bam \
   -R $TMPDIR/reference.fa \
   -knownSites $DBSNP_FILENAME \
   -knownSites $INDELS_1000G_FILENAME \
   -knownSites $INDELS_GOLDSTD_FILENAME \
   -BQSR $PRE_RECALIBRATION_REPORT \
   -o ${SAMPLE}.${FRAGMENT}.realigned.recalibrated.recal_data.grp \
   -L $TMPDIR/fragment.intervals \
   -rf BadCigar

echo "`${NOW}`copying recalibration report to $ANALYSIS_DIR/recalibration/reports/post..."
cp ${SAMPLE}.${FRAGMENT}.realigned.recalibrated.recal_data.grp $ANALYSIS_DIR/recalibration/reports/post/

#logging
STATUS=OK
if [[ ! -e $ANALYSIS_DIR/recalibration/reports/post/${SAMPLE}.${FRAGMENT}.realigned.recalibrated.recal_data.grp ]]
then
	STATUS=FAILED
fi

echo -e "`${NOW}`$SCRIPT_CODE\t$SAMPLE\tall\tpostrecalibration_report\t$STATUS" >> $RUN_LOG


echo "`${NOW}`done"

