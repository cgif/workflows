#!/bin/bash

## script to merge mutect chunks

#PBS -l walltime=24:00:00
#PBS -l select=1:ncpus=1:mem=1gb

#PBS -M igf@imperial.ac.uk
#PBS -m ea
#PBS -j oe

#PBS -q pqcgi

NOW="date +%Y-%m-%d%t%T%t"

module load java/#javaVersion
module load gatk/#gatkVersion
JAVA_XMX=4G
NT=2

# define variables
RESULTS_DIR=#results_dir
SAMPLE=`basename $RESULTS_DIR`
RAW_STATS_PATH="#rawStatsFiles"
IN_VCF="#rawVcfFiles"
REFERENCE_FASTA=#referenceFasta
REFERENCE_SEQ_DICT=#referenceSeqDict

SUMMARY_SCRIPT_PATH=#summaryScriptPath
BASEDIR=#baseDir

#copy reference to $TMP
cp $REFERENCE_FASTA $TMPDIR/reference.fa
cp $REFERENCE_FASTA.fai $TMPDIR/reference.fa.fai
cp $REFERENCE_SEQ_DICT $TMPDIR/reference.dict

# make tmp folder for temporary java files
mkdir $TMPDIR/tmp

VCF_FILES_COUNT=`echo $IN_VCF | perl -e '$in=<>; @tokens=split(/\s/,$in); $count=@tokens; print $count;'`

# if there is more than one input VCF file
if [ $VCF_FILES_COUNT -ge 2 ]; then

	# copy VCF files to be merged to temp space and get total number of variants in all vcf files before merging 
	TMP_IN_VCF=""
	VARIANT_COUNT_INPUT=0
       
	for VCF_FILE in $IN_VCF; do		
	
		VCF_BASENAME=`basename $VCF_FILE`
		grep -v REJECT $VCF_FILE > $TMPDIR/$VCF_BASENAME
		TMP_IN_VCF="$TMP_IN_VCF -V $TMPDIR/$VCF_BASENAME"

		# get number of variants in the input VCF file
		VARIANT_COUNT=`grep -v '#' $TMPDIR/$VCF_BASENAME | wc -l`
		VARIANT_COUNT_INPUT=$(($VARIANT_COUNT_INPUT + $VARIANT_COUNT))
	done

	java -Xmx$JAVA_XMX -XX:+UseSerialGC -Djava.io.tmpdir=$TMPDIR/tmp -jar $GATK_HOME/GenomeAnalysisTK.jar \
		-nt $NT \
		-R $TMPDIR/reference.fa \
		-T CombineVariants \
		$TMP_IN_VCF \
	        --assumeIdenticalSamples \
		-o $TMPDIR/merged.vcf

	#make tsv file for input into Oncotator web server
	echo "`${NOW}` making tsv file..."
	perl $BASEDIR/../../helper/vcf_to_oncotator_tsv.pl $TMPDIR/merged.vcf $TMPDIR/tmp.tsv


	# get number of reads in the output GVCF file
	VARIANT_COUNT_OUTPUT=`grep -v '#' $TMPDIR/merged.vcf | wc -l`
	if [[ $VARIANT_COUNT_INPUT -eq $VARIANT_COUNT_OUTPUT ]]; then

		cp $TMPDIR/merged.vcf $RESULTS_DIR/$SAMPLE.mutect.vcf
		cp $TMPDIR/tmp.tsv $RESULTS_DIR/$SAMPLE.mutect.tsv
		chmod 660 $RESULTS_DIR/$SAMPLE.mutect.vcf
		chmod 660 $RESULTS_DIR/$SAMPLE.mutect.tsv

	else

		echo "Output VCF does not contain the same number of variants as the input VCF files"

        fi

fi

# if there is only one input VCF nothing to merge. Just copy VCF to results directory
if [ $VCF_FILES_COUNT -eq 1 ]; then

	VCF_BASENAME=`basename $IN_VCF`
	grep -v REJECT $IN_VCF > $TMPDIR/mutect.vcf

	#make tsv file for input into Oncotator web server
	echo "`${NOW}` making tsv file..."
	perl $BASEDIR/../../helper/vcf_to_oncotator_tsv.pl $TMPDIR/mutect.vcf $TMPDIR/tmp.tsv


	cp $TMPDIR/mutect.vcf $RESULTS_DIR/$SAMPLE.mutect.vcf
	chmod 660 $RESULTS_DIR/$SAMPLE.mutect.vcf

	cp $TMPDIR/tmp.tsv $RESULTS_DIR/$SAMPLE.mutect.tsv
	chmod 660 $RESULTS_DIR/$SAMPLE.mutect.tsv

fi

RAW_STATS_MERGED=$RESULTS_DIR/$SAMPLE.stats

cat $RAW_STATS_PATH|grep -P '^contig\tposition\t'|head -n 1 > $RAW_STATS_MERGED

cat $RAW_STATS_PATH|grep -vP '^#'|grep -vP '^contig\tposition\t' >> $RAW_STATS_MERGED
 
grep -v REJECT $RAW_STATS_MERGED > $RAW_STATS_MERGED.keep

perl $SUMMARY_SCRIPT_PATH



