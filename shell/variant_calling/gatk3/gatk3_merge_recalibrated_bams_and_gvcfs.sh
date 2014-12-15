#!/bin/bash

## script to merge realigned and recalibrated BAM files

#PBS -l walltime=72:00:00
#PBS -l select=1:ncpus=2:mem=5gb

#PBS -M cgi@imperial.ac.uk
#PBS -m ea
#PBS -j oe

#PBS -q pqcgi


# load required modules
module load picard/#picardVersion
module load samtools/#samtoolsVersion
module load java/#javaVersion
module load gatk/#gatkVersion

#TABIX module not installed properly...
#module load tabix/0.2.6
#...use own copy for the time being
TABIX_HOME=/groupvol/cgi/bin/tabix-#tabixVersion

NOW="date +%Y-%m-%d%t%T%t"
JAVA_XMX=4G
NT=2

SCRIPT_CODE="GATKMEBG"

# define variables
RECAL_IN_BAM="#recalibratedBam"
HC_IN_BAM="#HCbamFiles"
INPUT_DIR_RECALIBRATED=#inputDirRecalibrated
PATH_OUTPUT_DIR_RECALIBRATED=#pathOutputDirRecalibrated
PATH_STRIP_INDEL_QUALS_SCRIPT=#pathStripIndelQualsScript
INPUT_DIR_HC_BAM=#inputDirHCbam
PATH_OUTPUT_DIR_HC_BAM=#pathOutputDirHCbam
SAMPLE=#sample
RUN_LOG=#runLog
IN_GVCFS="#genomicVCFs"
INPUT_DIR_GVCF=#inputDirGVCF
PATH_OUTPUT_DIR_GVCF=#pathOutputDirGVCF
REFERENCE_FASTA=#referenceFasta
REFRENCE_SEQ_DICT=`echo $REFERENCE_FASTA | perl -pe 's/\.fa/\.dict/'`
PROJECT=#projectName
CAPTURE=#captureKit
GVCF_LIST=#gvcfList

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
	
IN_GVCF_COUNT=`echo $IN_GVCFS | perl -e '$in=<>; @tokens=split(/\s/,$in); $count=@tokens; print $count;'`
TMP_PATH_OUT_GVCF=$TMPDIR/tmp.genomic.vcf
OUT_GVCF=$PATH_OUTPUT_DIR_GVCF/$SAMPLE.genomic.vcf

# if there is more than one input GVCF file
if [ $IN_GVCF_COUNT -ge 2 ]; then

	# copy GVCF files to be merged to temp space and get total number of variants in all gvcf files before merging 
	echo "`${NOW}`INFO $SCRIPT_CODE copying GVCF files to $TMPDIR..."
	TMP_IN_GVCF=""
	VARIANT_COUNT_INPUT=0

	for GVCF in $IN_GVCFS; do		
	
		GVCF_BASENAME=`basename $GVCF`
		echo "`${NOW}`INFO $SCRIPT_CODE $GVCF_BASENAME"
		cp $INPUT_DIR_GVCF/$GVCF_BASENAME $TMPDIR
#		TMP_IN_GVCF="$TMP_IN_GVCF -V $TMPDIR/$GVCF_BASENAME"	#this line if GATK is used for merging
		TMP_IN_GVCF="$TMP_IN_GVCF $TMPDIR/$GVCF_BASENAME"

		# get number of variants in the input GVCF file
		VARIANT_COUNT=`grep -v '#' $TMPDIR/$GVCF_BASENAME | wc -l`
		VARIANT_COUNT_INPUT=$(($VARIANT_COUNT_INPUT + $VARIANT_COUNT))
	done

	# merge GVCF files

	echo "`${NOW}`INFO $SCRIPT_CODE merging GVCF files..."
	# use simple concatenation instead of GATK, as GATK is incredibly slow for larger files
	# we can use simple concatenation, because the positions never overlap
	
	# get header
	CHUNK_1_GVCF=$(echo $TMP_IN_GVCF | cut -f1 -d" " )
	awk '{if ( $1 ~ /^#/ ) print}' $CHUNK_1_GVCF > $TMPDIR/merged.genomic.vcf

	# now concatenate files without headers
	for FILE in $TMP_IN_GVCF; do
    	awk '{if ( $1 !~ /^#/ ) print}' $FILE >> $TMPDIR/merged.genomic.vcf
	done


	### below is old version where GATK is used for merging

#	java -Xmx$JAVA_XMX -XX:+UseSerialGC -Djava.io.tmpdir=$TMPDIR/tmp -jar $GATK_HOME/GenomeAnalysisTK.jar \
#		-nt $NT \
#		-R $TMPDIR/reference.fa \
#		-T CombineVariants \
#		$TMP_IN_GVCF \
#		--assumeIdenticalSamples \
#		--genotypemergeoption UNSORTED \
#		-o $TMPDIR/merged.genomic.vcf

	# get number of reads in the output GVCF file
	VARIANT_COUNT_OUTPUT=`grep -v '#' $TMPDIR/merged.genomic.vcf | wc -l`

	# copy merged GVCF to destination folder - only copy if correct variant number
#	echo "`${NOW}`INFO $SCRIPT_CODE copying merged GVCF to $OUT_GVCF..."
#	cp $TMPDIR/merged.genomic.vcf $OUT_GVCF
#	chmod 660 $OUT_GVCF
	
	# check if read counts are the same and if so, copy merged GVCF to destination folder:
 	echo "`${NOW}`INFO $SCRIPT_CODE input  GVCFs read count: $VARIANT_COUNT_INPUT"
	echo "`${NOW}`INFO $SCRIPT_CODE output GVCF read count: $VARIANT_COUNT_OUTPUT"

	echo "`$NOW`INFO $SCRIPT_CODE compressing GVCF file..."
	$TABIX_HOME/bgzip $TMPDIR/merged.genomic.vcf

	echo "`$NOW`INFO $SCRIPT_CODE indexing GVCF file..."
	$TABIX_HOME/tabix -p vcf $TMPDIR/merged.genomic.vcf.gz
	
	if [[ $VARIANT_COUNT_INPUT -eq $VARIANT_COUNT_OUTPUT ]]; then
		echo "`${NOW}`INFO $SCRIPT_CODE copying merged GVCF to $OUT_GVCF.gz ..."
		cp $TMPDIR/merged.genomic.vcf.gz $OUT_GVCF.gz
		cp $TMPDIR/merged.genomic.vcf.gz.tbi $OUT_GVCF.gz.tbi
		chmod 660 $OUT_GVCF*

	else
		echo "`${NOW}`ERROR $SCRIPT_CODE Output GVCF does not contain the same number of variants as the input GVCF files, merged GVCF file not copied!"
	fi

	## loging 
	if [[ -s $OUT_GVCF.gz ]]; then 
		STATUS=OK

		echo -e "$SAMPLE\t$OUT_GVCF.gz\t$PROJECT\t$CAPTURE" >> $GVCF_LIST
		chmod 660 $GVCF_LIST

		echo "`${NOW}`INFO $SCRIPT_CODE deleting intermediate gVCF files..."
		for GVCF in $IN_GVCFS; do 
			rm $INPUT_DIR_GVCF/$GVCF
		done
	else 
		STATUS=FAILED
	fi

	echo -e "`${NOW}`$SCRIPT_CODE\t$SAMPLE\tall\tgenomic_vcf\t$STATUS" >> $RUN_LOG
fi
	
# if there is only one input GVCF
# nothing to merge. Just copy GVCF to results directory
if [ $IN_GVCF_COUNT -eq 1 ]; then

	echo "`${NOW}`INFO $SCRIPT_CODE only one input GVCF file. Nothing to merge."
	GVCF_BASENAME=`basename $IN_GVCFS`
	mv $INPUT_DIR_GVCF/$GVCF_BASENAME $OUT_GVCF

	echo "`$NOW`INFO $SCRIPT_CODE compressing GVCF file..."
	$TABIX_HOME/bgzip $OUT_GVCF

	echo "`$NOW`INFO $SCRIPT_CODE indexing GVCF file..."
	$TABIX_HOME/tabix -p vcf $OUT_GVCF.gz

	#logging
	if [[ ! -s $OUT_GVCF.gz ]]; then
		STATUS=FAILED
	else 
		STATUS=OK
		echo -e "$SAMPLE\t$OUT_GVCF.gz\t$PROJECT\t$CAPTURE" >> $GVCF_LIST
		chmod 660 $GVCF_LIST
	fi
	echo -e "`${NOW}`$SCRIPT_CODE\t$SAMPLE\tall\tgenomic_vcf\t$STATUS" >> $RUN_LOG
fi


##########################
## now merge BAM files ##
##########################

#for each set of BAM files
#we will merge realigned and recalibrated 
#and haplotype BAM files

for IN_BAM in "$RECAL_IN_BAM" "$HC_IN_BAM"; do 

	if [[ $IN_BAM == *recalibrated* ]]; then
		echo "`${NOW}`INFO $SCRIPT_CODE processing recalibrated BAM files..."
		INPUT_DIR=$INPUT_DIR_RECALIBRATED
		PATH_OUTPUT_DIR=$PATH_OUTPUT_DIR_RECALIBRATED

		#as we only merge realigned and recalibrated 
		#BAM files we will save the merged file simply
		#under the sample name
		OUT_BAM=$PATH_OUTPUT_DIR/$SAMPLE.bam
	elif [[ $IN_BAM == *hc.bam* ]]; then
		echo "`${NOW}`INFO $SCRIPT_CODE processing haplotype BAM files..."
		INPUT_DIR=$INPUT_DIR_HC_BAM
		PATH_OUTPUT_DIR=$PATH_OUTPUT_DIR_HC_BAM
		OUT_BAM=$PATH_OUTPUT_DIR/$SAMPLE.hc.bam
	fi

	IN_BAM_COUNT=`echo $IN_BAM | perl -e '$in=<>; @tokens=split(/\s/,$in); $count=@tokens; print $count;'`
	TMP_PATH_OUT_BAM=$TMPDIR/tmp.bam

	READ_COUNT_INPUT=0
	READ_COUNT_OUTPUT=0

	# if there is more than
	# one input BAM file
	if [ $IN_BAM_COUNT -ge 2 ]; then
		# copy BAM files to be merged to temp space and get total number of reads
		# in all bam files before 	merging
		echo "`${NOW}`INFO $SCRIPT_CODE copying BAM files to $TMPDIR..."
		TMP_IN_BAM=""
		INPUT_READ_COUNT=0

		for BAM in $IN_BAM; do		
			BAM_BASENAME=`basename $BAM`
			echo "`${NOW}`INFO $SCRIPT_CODE $BAM_BASENAME"
			cp $INPUT_DIR/$BAM_BASENAME $TMPDIR

			#remove indel quality scores from recalibrated BAM chunk files
			#we do it before merging so less temporary space is needed for merging
			echo "`${NOW}`INFO $SCRIPT_CODE stripping BAM file $BAM_BASENAME of indel quality scores..."
			$PATH_STRIP_INDEL_QUALS_SCRIPT -i $TMPDIR/$BAM_BASENAME -o $TMPDIR/$BAM_BASENAME.stripped
			mv $TMPDIR/$BAM_BASENAME.stripped $TMPDIR/$BAM_BASENAME

			TMP_IN_BAM="$TMP_IN_BAM INPUT=$TMPDIR/$BAM_BASENAME"

			#get number of reads in the input BAM file
			READ_COUNT=`samtools flagstat $TMPDIR/$BAM_BASENAME | head -n 1 | perl -e 'while(<>){ if(/(\d*?)\s\+\s(\d*?)\s/) { $retval=$1+$2; print "$retval\n"; }  }'`	
			READ_COUNT_INPUT=$(($READ_COUNT_INPUT + $READ_COUNT))
		done
	    
		# merge BAM files with picard tools
		# picard allows to merge unsorted BAM files and
		# output a coordinate sorted BAM while samtools will output
		# an unsorted BAM file if the input files are unsorted.
		# (Increasing the number of reads in memory before spilling to disc
		#  has no impact on the runtime! This was tested with 10M reads which
		# requires ~8GB of RAM and 4 processors. The default is 500k reads.)	
			
		echo "`${NOW}`INFO $SCRIPT_CODE merging and coordinate sorted BAM files..."
		java -jar -Xmx$JAVA_XMX -XX:+UseSerialGC $PICARD_HOME/MergeSamFiles.jar $TMP_IN_BAM OUTPUT=$TMP_PATH_OUT_BAM SORT_ORDER=coordinate USE_THREADING=true  VALIDATION_STRINGENCY=SILENT TMP_DIR=$TMPDIR
    
    	#remove indel quality scores from recalibrated BAM files
		#commented out, because do it before merging to make merging more efficient
#		echo "`${NOW}`INFO $SCRIPT_CODE stripping BAM file of indel quality scores..."
#		$PATH_STRIP_INDEL_QUALS_SCRIPT -i $TMP_PATH_OUT_BAM -o $TMP_PATH_OUT_BAM.stripped
#		mv $TMP_PATH_OUT_BAM.stripped $TMP_PATH_OUT_BAM

    
		# get number of reads in the output bam file
		READ_COUNT_OUTPUT=`samtools flagstat $TMP_PATH_OUT_BAM | head -n 1 | perl -e 'while(<>){ if(/(\d*?)\s\+\s(\d*?)\s/) { $retval=$1+$2; print "$retval\n"; }  }'`	
    
		# index output BAM
		echo "`${NOW}`INFO $SCRIPT_CODE indexing merged BAM..."
		samtools index $TMP_PATH_OUT_BAM
	    
		# copy merged BAM and index
		# to destination folder

		if [[ $READ_COUNT_INPUT -eq $READ_COUNT_OUTPUT ]]; then
			echo "`${NOW}`INFO $SCRIPT_CODE copying merged BAM to $OUT_BAM..."
			cp $TMP_PATH_OUT_BAM $OUT_BAM
			chmod 660 $OUT_BAM
         
			echo "`${NOW}`INFO $SCRIPT_CODE copying BAM index to $OUT_BAM.bai"
			cp $TMP_PATH_OUT_BAM.bai $OUT_BAM.bai
			chmod 660 $OUT_BAM.bai

			echo "`${NOW}`INFO $SCRIPT_CODE deleting intermediate BAM files..."
			for BAM_FILE in $IN_BAM; do
				rm $INPUT_DIR/$BAM_FILE
				rm $INPUT_DIR/$BAM_FILE.bai
			done
			#...if no, keep input BAM files for re-run
		else
			echo "`${NOW}`WARN $SCRIPT_CODE Output BAM does not contain the same number of reads as the input BAM file(s)!"
			echo "`${NOW}`WARN $SCRIPT_CODE Keeping intermediate BAM files for re-run..."  		 	
		fi

  
    		#logging
		BAM=recalibrated_bam_merged
    	
		STATUS=OK
		if [[ ! -s $OUT_BAM ]]; then
			STATUS=FAILED
		fi
		echo -e "`${NOW}`$SCRIPT_CODE\t$SAMPLE\tall\t$BAM\t$STATUS" >> $RUN_LOG    
	fi

	# if there is only one input BAM
	# nothing to merge. Just copy BAM
	# file to temp space for sorting, indexing
	# flagstat and md5
	if [ $IN_BAM_COUNT -eq 1 ]; then

		echo "`${NOW}`INFO $SCRIPT_CODE only one input BAM file. Nothing to merge."
		BAM_BASENAME=`basename $IN_BAM`
		IN_BAM=$INPUT_DIR/$BAM_BASENAME
   
		#remove indel quality scores from recalibrated BAM files
	
		echo "`${NOW}`INFO $SCRIPT_CODE copying input BAM to $TMPDIR..."
		cp $IN_BAM $TMPDIR

		# get number of reads in the input bam file
		READ_COUNT_INPUT=`samtools flagstat $BAM_BASENAME | head -n 1 | perl -e 'while(<>){ if(/(\d*?)\s\+\s(\d*?)\s/) { $retval=$1+$2; print "$retval\n"; }  }'`	
	    			
		#strip indel qualities
		echo "`${NOW}`INFO $SCRIPT_CODE stripping BAM file of indel quality scores..."
		$PATH_STRIP_INDEL_QUALS_SCRIPT -i $BAM_BASENAME -o $TMP_PATH_OUT_BAM.stripped
		mv $TMP_PATH_OUT_BAM.stripped $TMP_PATH_OUT_BAM
  
		# get number of reads in the output bam file
		READ_COUNT_OUTPUT=`samtools flagstat $TMP_PATH_OUT_BAM | head -n 1 | perl -e 'while(<>){ if(/(\d*?)\s\+\s(\d*?)\s/) { $retval=$1+$2; print "$retval\n"; }  }'`	
	    	
		# index output BAM
		echo "`${NOW}`INFO $SCRIPT_CODE indexing stripped BAM..."
		samtools index $TMP_PATH_OUT_BAM

		if [[ $READ_COUNT_INPUT -eq $READ_COUNT_OUTPUT ]]; then
			echo "`${NOW}`INFO $SCRIPT_CODE copying merged BAM to $OUT_BAM..."
			cp $TMP_PATH_OUT_BAM $OUT_BAM
			chmod 660 $OUT_BAM
         
			echo "`${NOW}`INFO $SCRIPT_CODE copying BAM index to $OUT_BAM.bai"
			cp $TMP_PATH_OUT_BAM.bai $OUT_BAM.bai
			chmod 660 $OUT_BAM.bai

			echo "`${NOW}`INFO $SCRIPT_CODE deleting intermediate BAM files..."
			rm $IN_BAM
			rm $IN_BAM.bai

			#...if no, keep input BAM files for re-run
		else
			echo "`${NOW}`WARN $SCRIPT_CODE Output BAM does not contain the same number of reads as the input BAM file(s)!"
			echo "`${NOW}`WARN $SCRIPT_CODE Keeping intermediate BAM files for re-run..."  		 
	
		fi

    		#logging
		BAM=recalibrated_bam_merged
    	
		STATUS=OK
		if [[ ! -s $OUT_BAM ]]; then
		STATUS=FAILED
		fi
		echo -e "`${NOW}`$SCRIPT_CODE\t$SAMPLE\tall\t$BAM\t$STATUS" >> $RUN_LOG
    	   
	fi

	echo "`${NOW}`INFO $SCRIPT_CODE input  BAM(s) read count: $READ_COUNT_INPUT"
	echo "`${NOW}`INFO $SCRIPT_CODE output BAM    read count: $READ_COUNT_OUTPUT"

done	

