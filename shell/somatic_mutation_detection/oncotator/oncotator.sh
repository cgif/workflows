#!/bin/bash

## script to merge vcf files across all samples and to run Oncotator on the merged file

#PBS -l walltime=24:00:00
#PBS -l select=1:ncpus=1:mem=10gb

#PBS -M cgi@imperial.ac.uk
#PBS -m ea
#PBS -j oe

#PBS -q pqcgi

module load oncotator/#oncotatorVersion
module load python/#pythonVersion
module load samtools/#samtoolsVersion
module load tabix/#tabixVersion

module load gatk/#gatkVersion
module load java/#javaVersion
module load R/#rVersion

JAVA_XMX=4G
NT=2

IN_VCF="#multisampleVCF"
PROJECT=#project
RESULTS_DIR=#resultsDir
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

# make tmp folder for temporary java files
mkdir $TMPDIR/tmp

#merge vcf files across all samples
echo "`$NOW` mergeing VCF files across all samples"

VCF_FILES_COUNT=`echo $IN_VCF | perl -e '$in=<>; @tokens=split(/\s/,$in); $count=@tokens; print $count;'`

# if there is more than one input VCF file
if [ $VCF_FILES_COUNT -ge 2 ]; then

	# copy VCF files to be merged to temp space 
	TMP_IN_VCF=""
       
	for VCF_FILE in $IN_VCF; do		
	
		VCF_BASENAME=`basename $VCF_FILE`
		cp $VCF_FILE $TMPDIR/$VCF_BASENAME
		TMP_IN_VCF="$TMP_IN_VCF -V $TMPDIR/$VCF_BASENAME"

	done

	java -Xmx$JAVA_XMX -XX:+UseSerialGC -Djava.io.tmpdir=$TMPDIR/tmp -jar $GATK_HOME/GenomeAnalysisTK.jar \
		-nt $NT \
		-R $TMPDIR/reference.fa \
		-T CombineVariants \
		$TMP_IN_VCF \
		-o $TMPDIR/merged.vcf

elif [ $VCF_FILES_COUNT -eq 1 ]; then

	cp $IN_VCF $TMPDIR/merged.vcf

fi

cp $TMPDIR/merged.vcf $RESULTS_DIR/${PROJECT}.vcf
chmod 660 $RESULTS_DIR/${PROJECT}.vcf

echo "`$NOW` run Oncotator on VCF file"
oncotator -i VCF --db-dir $TMPDIR/oncotator_db -o VCF $TMPDIR/merged.vcf $TMPDIR/annotated.vcf hg19
cat $TMPDIR/annotated.vcf |perl -e 'while (<>) { s/(^\S+\t\S+\t\S+\t\S+\t\S+\t\S+\t\S+\t)\"(\S+)\"(\t.*)/$1$2$3/g; s/^\<(\S+)\>(.*)/$1$2/; print; }' > $RESULTS_DIR/${PROJECT}.annotated.vcf
chmod 660 $RESULTS_DIR/${PROJECT}.annotated.vcf

oncotator -i VCF --db-dir $TMPDIR/oncotator_db -o TCGAMAF --skip-no-alt $TMPDIR/merged.vcf $TMPDIR/annotated.maf hg19
cp $TMPDIR/annotated.maf $RESULTS_DIR/${PROJECT}.annotated.maf
chmod 660 $RESULTS_DIR/${PROJECT}.annotated.maf

