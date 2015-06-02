#!/bin/bash

## script to merge multisample GVCF chunks into one file

#PBS -l walltime=72:00:00
#PBS -l select=1:ncpus=1:mem=5gb:tmpspace=#tmpSpacegb

#PBS -M cgi@imperial.ac.uk
#PBS -m ea
#PBS -j oe

#PBS -q pqcgi


# load required modules
module load java/#javaVersion
module load gatk/#gatkVersion

NOW="date +%Y-%m-%d%t%T%t"
JAVA_XMX=4G

# define variables
REFERENCE_FASTA=#referenceFasta
REFRENCE_SEQ_DICT=`echo $REFERENCE_FASTA | perl -pe 's/\.fa/\.dict/'`
PREFIX=#prefixName
IN_GVCFS="#chunkGVCFs"
INPUT_DIR_GVCF=#inputDirGVCF
RESULTS_DIR=#resultsDir

#############################
## merge genomic VCF files ##
#############################

#copy reference to $TMP
echo "`${NOW}`INFO $SCRIPT_CODE copying reference fasta and index to tmp directory..."
cp $REFERENCE_FASTA $TMPDIR/reference.fa
cp $REFERENCE_FASTA.fai $TMPDIR/reference.fa.fai
cp $REFRENCE_SEQ_DICT $TMPDIR/reference.dict

# make tmp folder for temporary java files
mkdir $TMPDIR/tmp

OUT_GVCF=$RESULTS_DIR/$PREFIX.combined.genomic.vcf

# copy GVCF files to be merged to temp space
echo "`${NOW}`INFO copying GVCF files to $TMPDIR..."
TMP_IN_GVCF=""

for GVCF in $IN_GVCFS; do		

	GVCF_BASENAME=`basename $GVCF`
	echo "`${NOW}`INFO copying $GVCF_BASENAME to tmp space"
	cp $INPUT_DIR_GVCF/$GVCF_BASENAME $TMPDIR
	echo "`${NOW}`INFO copying $GVCF_BASENAME.idx to tmp space"	
	cp $INPUT_DIR_GVCF/$GVCF_BASENAME.idx $TMPDIR
	TMP_IN_GVCF="$TMP_IN_GVCF -V $TMPDIR/$GVCF_BASENAME"	

done

	# merge GVCF files

echo "`${NOW}`INFO $SCRIPT_CODE merging GVCF files..."

java -Xmx$JAVA_XMX -XX:+UseSerialGC -Djava.io.tmpdir=$TMPDIR/tmp -cp $GATK_HOME/GenomeAnalysisTK.jar org.broadinstitute.gatk.tools.CatVariants \
	-R $TMPDIR/reference.fa \
	$TMP_IN_GVCF \
	--assumeSorted \
	-out $TMPDIR/merged.combined.genomic.vcf

	# copy merged GVCF to destination folder
	
echo "`${NOW}`INFO copying merged GVCF to $OUT_GVCF ..."
cp $TMPDIR/merged.combined.genomic.vcf $OUT_GVCF							
cp $TMPDIR/merged.combined.genomic.vcf.idx $OUT_GVCF.idx
chmod 660 $OUT_GVCF*

ls -al $TMPDIR

du -h $TMPDIR

echo -e "`${NOW}`Done"
	

