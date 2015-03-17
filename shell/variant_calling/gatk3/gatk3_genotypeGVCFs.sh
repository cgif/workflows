#!/bin/bash

## script to run GATK unified genotyper for genome

#PBS -l walltime=24:00:00
#PBS -l select=1:ncpus=#dataThreads:mem=16gb

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
JAVA_XMX=15G

# define variables

DBSNP=#dbSnp
REFERENCE_FASTA=#referenceFasta
REFRENCE_SEQ_DICT=`echo $REFERENCE_FASTA | perl -pe 's/\.fa/\.dict/'`
ANALYSIS_DIR=#analysisDir
SAMPLE=#sampleName		### this is really a project name
FRAGMENT=#fragmentName
PED_FILE=#pedFile
GG_DATATHREADS=#dataThreads
FRAGMENT_FILE=#fragmentFile
TYPE=#seqType
RUN_LOG=#runLog
AUX_LIST=#auxList


SCRIPT_CODE="GATKGGVC"


LOG_INFO="`${NOW}`INFO $SCRIPT_CODE"
LOG_ERR="`$NOW`ERROR $SCRIPT_CODE"
LOG_WARN="`$NOW`WARN $SCRIPT_CODE"
LOG_DEBUG="`$NOW`DEBUG $SCRIPT_CODE"

#copy input files to tmp dir

echo "`${NOW}`INFO $SCRIPT_CODE copying GATK resources to tmp directory..."
DBSNP_FILENAME=`basename $DBSNP`

echo "`${NOW}`INFO $SCRIPT_CODE $DBSNP"
cp $DBSNP $TMPDIR/$DBSNP_FILENAME
cp $DBSNP.idx $TMPDIR/$DBSNP_FILENAME.idx


cp $FRAGMENT_FILE $TMPDIR/fragment.intervals

echo "`${NOW}`INFO $SCRIPT_CODE copying chunk gVCF file to tmp directory..."

#get analysis directory of GATK3 run
ANALYSIS_DIR_RUN=`dirname $ANALYSIS_DIR`

#get sample names from analysis directory and copy
#respective genomic VCF files
INPUT_GVCF_ARGUMENT=""

for SAMPLE_NAME in `ls --color=never $ANALYSIS_DIR_RUN | grep -vE "multisample|aux"`; do

	GVCF="$ANALYSIS_DIR_RUN/$SAMPLE_NAME/haplotypecaller/$SAMPLE_NAME.$FRAGMENT.genomic.vcf"
	GVCF_NAME=`basename $GVCF`
	
	echo "`${NOW}`INFO $SCRIPT_CODE copying $GVCF to temporary space"
	cp $GVCF $TMPDIR
	
	INPUT_GVCF_ARGUMENT="$INPUT_GVCF_ARGUMENT	-V $TMPDIR/$GVCF_NAME"

done;

if [[ $AUX_LIST != "" ]]; then 
	sort $AUX_LIST | uniq | while read SAMPLE GVCF; do
    		if [[ "$SAMPLE" != "" ]]; then

			GVCF_NAME=`basename $GVCF`
			echo "`${NOW}`INFO $SCRIPT_CODE copying $GVCF to temporary space"
			cp $GVCF $TMPDIR
			cp $GVCF.tbi $TMPDIR
			
			INPUT_GVCF_ARGUMENT="$INPUT_GVCF_ARGUMENT	-V $TMPDIR/$GVCF_NAME"

		fi
		### need to write into a file, because bash runs while loop in a subshell 
		### and variables are not available outside the loop
		echo "$INPUT_GVCF_ARGUMENT" > $TMPDIR/tmp.gvcf.arg  							
	done;
	INPUT_GVCF_ARGUMENT=`cat $TMPDIR/tmp.gvcf.arg`
fi


echo "`${NOW}`INFO $SCRIPT_CODE copying reference fasta and index to tmp directory..."
cp $REFERENCE_FASTA $TMPDIR/reference.fa
cp $REFERENCE_FASTA.fai $TMPDIR/reference.fa.fai
cp $REFRENCE_SEQ_DICT $TMPDIR/reference.dict


# make tmp folder for temporary java files
mkdir $TMPDIR/tmp


#PED_ARG="" 											# deactivated because GenotypeGVCF does not take into account family structure

##### check for ped file
#if [[ "$PED_FILE" != "none" ]]; then
#	PED_ARG="-ped $PED_FILE -pedValidationType SILENT"
#fi

############ genotype gVCFs #############################
# default downsampling to 1000 left 

CUT_OFF=30
#if [[ "$TYPE" != "TARGETED" ]]; then
#	CUT_OFF=20
#fi

## list all annotations which are required for VSQR in the next step, because by some reason some defaul annotations are not emmited unless specified in the command (QD and FS, for example)

echo "`${NOW}`INFO $SCRIPT_CODE running GenotypeGVCFs"
java -Xmx$JAVA_XMX -XX:+UseSerialGC -Djava.io.tmpdir=$TMPDIR/tmp -jar $GATK_HOME/GenomeAnalysisTK.jar \
	-nt $GG_DATATHREADS \
	-R $TMPDIR/reference.fa \
	-T GenotypeGVCFs \
	--dbsnp $DBSNP_FILENAME \
	$INPUT_GVCF_ARGUMENT \
	-o $TMPDIR/$SAMPLE.$FRAGMENT.raw.vcf \
	-stand_call_conf $CUT_OFF \
	-stand_emit_conf 10 \
	-A AlleleBalance \
	-A QualByDepth \
	-A FisherStrand \
	-A Coverage \
	-A ReadPosRankSumTest \
	-A MappingQualityRankSumTest \
	-A InbreedingCoeff \
	-A StrandOddsRatio \
	-L $TMPDIR/fragment.intervals

echo "`${NOW}`INFO $SCRIPT_CODE copying raw VCF file to output directory $ANALYSIS_DIR/genotypeGVCFs..."
cp $SAMPLE.$FRAGMENT.raw.vcf $ANALYSIS_DIR/genotypeGVCFs/

echo "`${NOW}`INFO $SCRIPT_CODE done"

#logging
STATUS=OK
if [[ ! -s $ANALYSIS_DIR/genotypeGVCFs/$SAMPLE.$FRAGMENT.raw.vcf ]]; then
	STATUS=FAILED
fi

echo -e "`${NOW}`$SCRIPT_CODE\tmultisample\t$FRAGMENT\traw_vcf\t$STATUS" >> $RUN_LOG

if [[ "$STATUS" == "FAILED" ]]; then
	exit 1;
fi

echo -e "`${NOW}` done"
