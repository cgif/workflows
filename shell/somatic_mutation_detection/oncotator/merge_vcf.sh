#!/bin/bash

## script to run merge vcf files for the same sample

#PBS -l walltime=24:00:00
#PBS -l select=1:ncpus=1:mem=5gb

#PBS -M cgi@imperial.ac.uk
#PBS -m ea
#PBS -j oe

#PBS -q pqcgi

module load oncotator/#oncotatorVersion
module load python/#pythonVersion
module load samtools/#samtoolsVersion
module load tabix/#tabixVersion

module load java/#javaVersion
module load gatk/#gatkVersion
module load R/#rVersion

JAVA_XMX=4G
NT=2

INPUT_VCF_MUTECT=#pathVCFmutect
INPUT_VCF_SID=#pathVCFsid
ANALYSIS_DIR=#analysisDir
SAMPLE_NAME=`basename $INPUT_VCF_SID .vcf`
REFERENCE_FASTA=#referenceFasta
REFRENCE_SEQ_DICT=`echo $REFERENCE_FASTA | perl -pe 's/\.fa/\.dict/'`
ONCOTATOR_DB=#oncotatorDB

#copy reference to $TMP
cp $REFERENCE_FASTA $TMPDIR/reference.fa
cp $REFERENCE_FASTA.fai $TMPDIR/reference.fa.fai
cp $REFRENCE_SEQ_DICT $TMPDIR/reference.dict

mkdir $TMPDIR/oncotator_db
cp -R $ONCOTATOR_DB/* $TMPDIR/oncotator_db
ls -l $TMPDIR/oncotator_db

#copy vcf
cp $INPUT_VCF_MUTECT $TMPDIR/mutect.vcf
cp $INPUT_VCF_SID $TMPDIR/SID.vcf

# make tmp folder for temporary java files
mkdir $TMPDIR/tmp

# get number of variants in the input VCF file
VARIANT_COUNT_MUTECT=`grep -v '#' $TMPDIR/mutect.vcf | wc -l`
VARIANT_COUNT_SID=`grep -v '#' $TMPDIR/SID.vcf | wc -l`
VARIANT_COUNT_INPUT=$(($VARIANT_COUNT_MUTECT + $VARIANT_COUNT_SID))

echo "`$NOW`merging Mutect and SomaticIndelDetector VCF files"
java -Xmx$JAVA_XMX -XX:+UseSerialGC -Djava.io.tmpdir=$TMPDIR/tmp -jar $GATK_HOME/GenomeAnalysisTK.jar \
        -nt $NT \
        -R $TMPDIR/reference.fa \
        -T CombineVariants \
        -V $TMPDIR/mutect.vcf \
        -V $TMPDIR/SID.vcf \
	--assumeIdenticalSamples \
	-o $TMPDIR/merged.vcf

echo "`$NOW` run Oncotator on VCF files"
oncotator -i VCF --db-dir $TMPDIR/oncotator_db -o VCF $TMPDIR/merged.vcf $TMPDIR/annotated.vcf hg19

echo "`$NOW` copy VCF files"
# get number of reads in the output VCF file
VARIANT_COUNT_OUTPUT=`grep -v '#' $TMPDIR/annotated.vcf | wc -l`

if [[ $VARIANT_COUNT_INPUT -le $VARIANT_COUNT_OUTPUT ]]; then

        cp $TMPDIR/merged.vcf $ANALYSIS_DIR/${SAMPLE_NAME}.vcf
	chmod 660 $ANALYSIS_DIR/${SAMPLE_NAME}.vcf

	cat $TMPDIR/annotated.vcf |perl -e 'while (<>) { s/(^\S+\t\S+\t\S+\t\S+\t\S+\t\S+\t\S+\t)\"(\S+)\"(\t.*)/$1$2$3/g; s/^\<(\S+)\>(.*)/$1$2/; print; }' > $ANALYSIS_DIR/${SAMPLE_NAME}.annotated.vcf
	chmod 660 $ANALYSIS_DIR/${SAMPLE_NAME}.annotated.vcf

else

	echo "Output VCF does not contain the same number of variants as the input VCF files"
	exit 1

fi
