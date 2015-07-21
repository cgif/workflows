#!/bin/bash

## script to run GATK variant recalibration and filtering for whole genomes
## based on GATK best practise guidelines
## http://www.broadinstitute.org/gsa/wiki/index.php/Best_Practice_Variant_Detection_with_the_GATK_v3

#PBS -l walltime=72:00:00
#PBS -l select=1:ncpus=4:mem=32gb

#PBS -M cgi@imperial.ac.uk
#PBS -m ea
#PBS -j oe

#PBS -q pqcgi

# load modules

module load gatk/#gatkVersion
module load R/#rVersion
module load java/#javaVersion

#TABIX module not installed properly...
#module load tabix/0.2.6
#...use own copy for the time being
TABIX_HOME=/groupvol/cgi/bin/tabix-#tabixVersion

NOW="date +%Y-%m-%d%t%T%t"
JAVA_XMX=31G
NT=4

SCRIPT_CODE="GATKREVC"

LOG_INFO="`$NOW`INFO $SCRIPT_CODE"
LOG_ERR="`$NOW`ERROR $SCRIPT_CODE"
LOG_WARN="`$NOW`WARN $SCRIPT_CODE"
LOG_DEBUG="`$NOW`DEBUG $SCRIPT_CODE"

# define variables
REFERENCE_FASTA=#referenceFasta
REFRENCE_SEQ_DICT=`echo $REFERENCE_FASTA | perl -pe 's/\.fa/\.dict/'`
ANALYSIS_DIR=#analysisDir
RESULTS_DIR=#resultsDir
SAMPLE=#sampleName
RAW_VCF_FILES="#rawVcfFiles"
SEQUENCING_TYPE=#sequencingType
SAMPLE_LIST=#sampleList
SUMMARY_SCRIPT_PATH=#summaryScriptPath
RUN_LOG=#runLog
OUTPUT_SAMPLE_LIST=#outputSampleList

#GATK resources
OMNI_1000G=#Omni1000G
HIGH_CONF_SNP_1000G=#highConfSnp1000G
INDELS_GOLDSTD=#indelsGoldStd
DBSNP=#dbSnp
DBSNP_EX_POST_129=#dbSnpExPost129
HAPMAP_SITES=#hapmapSites

echo "`${NOW}`INFO $SCRIPT_CODE copying GATK resources to tmp directory..."
OMNI_1000G_FILENAME=`basename $OMNI_1000G`
HIGH_CONF_SNP_1000G_FILENAME=`basename $HIGH_CONF_SNP_1000G`
INDELS_GOLDSTD_FILENAME=`basename $INDELS_GOLDSTD`
DBSNP_FILENAME=`basename $DBSNP`
DBSNP_EX_POST_129_FILENAME=`basename $DBSNP_EX_POST_129`
HAPMAP_SITES_FILENAME=`basename $HAPMAP_SITES`


echo "`${NOW}`$OMNI_1000G"
cp $OMNI_1000G $TMPDIR/$OMNI_1000G_FILENAME
cp $OMNI_1000G.idx $TMPDIR/$OMNI_1000G_FILENAME.idx

echo "`${NOW}`$INDELS_GOLDSTD"
cp $INDELS_GOLDSTD $TMPDIR/$INDELS_GOLDSTD_FILENAME
cp $INDELS_GOLDSTD.idx $TMPDIR/$INDELS_GOLDSTD_FILENAME.idx

echo "`${NOW}`$HIGH_CONF_SNP_1000G"
cp $HIGH_CONF_SNP_1000G $TMPDIR/$HIGH_CONF_SNP_1000G_FILENAME
cp $HIGH_CONF_SNP_1000G.idx $TMPDIR/$HIGH_CONF_SNP_1000G_FILENAME.idx

echo "`${NOW}`$DBSNP"
cp $DBSNP $TMPDIR/$DBSNP_FILENAME
cp $DBSNP.idx $TMPDIR/$DBSNP_FILENAME.idx

echo "`${NOW}`$DBSNP_EX_POST_129"
cp $DBSNP_EX_POST_129 $TMPDIR/$DBSNP_EX_POST_129_FILENAME
cp $DBSNP_EX_POST_129.idx $TMPDIR/$DBSNP_EX_POST_129_FILENAME.idx

echo "`${NOW}`$HAPMAP_SITES"
cp $HAPMAP_SITES $TMPDIR/$HAPMAP_SITES_FILENAME
cp $HAPMAP_SITES.idx $TMPDIR/$HAPMAP_SITES_FILENAME.idx


echo "`${NOW}`INFO $SCRIPT_CODE copying chunk VCF files to tmp directory..."
TMP_RAW_VCF_FILES=""
VCF_COUNT=0
for FILE in $RAW_VCF_FILES; do
	
	VCF_COUNT=$(( $VCF_COUNT + 1 ))
	TMP_VCF=$TMPDIR/chunk_$VCF_COUNT.raw.vcf	
	cp $FILE $TMP_VCF
	TMP_RAW_VCF_FILES="$TMP_RAW_VCF_FILES $TMP_VCF"

done

echo "`${NOW}`INFO $SCRIPT_CODE copying reference fasta and indexto tmp directory..."
cp $REFERENCE_FASTA $TMPDIR/reference.fa
cp $REFERENCE_FASTA.fai $TMPDIR/reference.fa.fai
cp $REFRENCE_SEQ_DICT $TMPDIR/reference.dict


# make tmp folder for temporary java files
mkdir $TMPDIR/tmp

# step 1: merge raw VCF files of all chomosomes
# only need to concatenate files because no overlapping records are present
# write a header line for a combined VCFfile

echo "`$NOW`INFO $SCRIPT_CODE merging raw VCF files..."

awk '{if ( $1 ~ /^#/ ) print}' $TMPDIR/chunk_1.raw.vcf > $TMPDIR/raw_merged.vcf

# concatenate all files excluding the header lines
for FILE in $TMP_RAW_VCF_FILES; do
    awk '{if ( $1 !~ /^#/ ) print}' $FILE >> $TMPDIR/raw_merged.vcf
done

#remove additional samples from VCF file
#(the --excludeNonVariants option removes loci found 
#to be non-variant after the subsetting procedure)

cut -f 1 $OUTPUT_SAMPLE_LIST > $TMPDIR/output_samples.txt

echo "`${NOW}`INFO $SCRIPT_CODE removing variants from auxiliary samples..."
java -Xmx$JAVA_XMX -XX:+UseSerialGC -Djava.io.tmpdir=$TMPDIR/tmp -jar $GATK_HOME/GenomeAnalysisTK.jar \
  -T SelectVariants \
  -R $TMPDIR/reference.fa \
  --variant $TMPDIR/raw_merged.vcf \
  -o $TMPDIR/raw_merged.noaux.vcf \
  --sample_file $TMPDIR/output_samples.txt \
  --excludeNonVariants

echo "`$NOW`INFO $SCRIPT_CODE compressing merged, filtered VCF file..."
$TABIX_HOME/bgzip -c $TMPDIR/raw_merged.noaux.vcf > $TMPDIR/raw_merged.noaux.vcf.gz

echo "`$NOW`INFO $SCRIPT_CODE indexing merged VCF file..."
$TABIX_HOME/tabix -p vcf $TMPDIR/raw_merged.noaux.vcf.gz

echo  "`$NOW`INFO $SCRIPT_CODE copying merged VCF file and index to $RESULTS_DIR..."
cp $TMPDIR/raw_merged.noaux.vcf.gz $RESULTS_DIR/$SAMPLE.raw.vcf.gz
cp $TMPDIR/raw_merged.noaux.vcf.gz.tbi $RESULTS_DIR/$SAMPLE.raw.vcf.gz.tbi

FINAL_VCF="$TMPDIR/raw_merged.noaux.vcf.gz"

#logging
STATUS=OK
if [[ ! -s $RESULTS_DIR/$SAMPLE.raw.vcf.gz ]]
then
	STATUS=FAILED
fi

echo -e "`${NOW}`$SCRIPT_CODE\tmultisample\tall\traw_vcf\t$STATUS" >> $RUN_LOG

#skip variant recalibration for targeted (amplicon) sequencing
#the number of variants will be to small to build a recalibration
#model
if [[ "$SEQUENCING_TYPE" != "TARGETED" ]]
then

	# step 2: snp.model <- BuildErrorModelWithVQSR(raw.vcf, SNP)
	# Arguments for VariantRecalibrator command taken from GATK best practise guidelines v3 for whole genomes
	# http://gatkforums.broadinstitute.org/discussion/1259/what-vqsr-training-sets-arguments-should-i-use-for-my-specific-project
	# removed inbreeding coef -an InbreedingCoeff, because it is a population level statistic
        # parameters were updated according to The GATK Guide Book for Version 3.1-1: 
        # annotations include -an DP for whole genome sequencing projects
	# added --target_titv, which is used for plotting purposes only 

	echo "`$NOW`INFO $SCRIPT_CODE building SNP recalibration model..."
	ANNOTATIONS_SNP="-an QD -an MQ -an MQRankSum -an ReadPosRankSum -an SOR -an FS "
	if [[ "$SEQUENCING_TYPE" == "WGS" ]]; then
		ANNOTATIONS_SNP="$ANNOTATIONS_SNP -an DP"
		TARGET_TITV="2.15"
	elif  [[ "$SEQUENCING_TYPE" == "EXOME" ]]; then
		TARGET_TITV="3.2"
	else
		echo "ERROR: sequencing type contains wrong value: $SEQUENCING_TYPE"
		exit 1
	fi

	java -Xmx$JAVA_XMX -XX:+UseSerialGC -Djava.io.tmpdir=$TMPDIR/tmp -jar $GATK_HOME/GenomeAnalysisTK.jar \
		-T VariantRecalibrator \
		-nt $NT \
		-R $TMPDIR/reference.fa \
		-input $TMPDIR/raw_merged.vcf \
                -resource:hapmap,known=false,training=true,truth=true,prior=15.0 $HAPMAP_SITES_FILENAME \
                -resource:omni,known=false,training=true,truth=true,prior=12.0 $OMNI_1000G_FILENAME \
                -resource:1000G,known=false,training=true,truth=false,prior=10.0 $HIGH_CONF_SNP_1000G_FILENAME \
                -resource:dbsnp,known=true,training=false,truth=false,prior=2.0 $DBSNP_FILENAME \
		$ANNOTATIONS_SNP \
		-mode SNP \
		-recalFile $TMPDIR/$SAMPLE.SNP.recal \
		-tranchesFile $TMPDIR/$SAMPLE.SNP.tranches \
		-rscriptFile $TMPDIR/$SAMPLE.SNP.plots.R \
		--target_titv $TARGET_TITV

	echo "`$NOW`INFO $SCRIPT_CODE copying SNP recalibration reports to $RESULTS_DIR/recalibration..."
	cp $SAMPLE.SNP.recal $RESULTS_DIR/recalibration
	cp $SAMPLE.SNP.tranches $RESULTS_DIR/recalibration
	cp $SAMPLE.SNP.tranches.pdf $RESULTS_DIR/recalibration
	cp $SAMPLE.SNP.plots.R $RESULTS_DIR/recalibration
	cp $SAMPLE.SNP.plots.R.pdf $RESULTS_DIR/recalibration

	# step 3: indel.model <- BuildErrorModelWithVQSR(raw.vcf, INDEL)
	# Arguments for VariantRecalibrator command taken from GATK best practise guidelines for whole genomes
	# http://gatkforums.broadinstitute.org/discussion/1259/what-vqsr-training-sets-arguments-should-i-use-for-my-specific-project
	# removed inbreeding coef -an InbreedingCoeff , because it is a population level statistic
        # parameters were updated according to The GATK Guide Book for Version 3.1-1:
        # -an DP only taken into account for whole genome sequencing

	echo "`$NOW`INFO $SCRIPT_CODE building INDEL recalibration model"
        ANNOTATIONS_INDEL="-an QD -an FS -an SOR -an ReadPosRankSum -an MQRankSum"
        if [[ "$SEQUENCING_TYPE" = "WGS" ]]
        then
                ANNOTATIONS_INDEL="$ANNOTATIONS_INDEL -an DP"
        fi

	java -Xmx$JAVA_XMX -XX:+UseSerialGC -Djava.io.tmpdir=$TMPDIR/tmp -jar $GATK_HOME/GenomeAnalysisTK.jar \
		-T VariantRecalibrator \
		-nt $NT \
		-R $TMPDIR/reference.fa \
		-input $TMPDIR/raw_merged.vcf \
                --maxGaussians 4 \
		-resource:mills,VCF,known=false,training=true,truth=true,prior=12.0 $INDELS_GOLDSTD_FILENAME \
		-resource:dbsnp,VCF,known=true,training=false,truth=false,prior=2.0 $DBSNP_FILENAME \
		$ANNOTATIONS_INDEL \
		-mode INDEL \
		-recalFile $TMPDIR/$SAMPLE.INDEL.recal \
		-tranchesFile $TMPDIR/$SAMPLE.INDEL.tranches \
		-rscriptFile $TMPDIR/$SAMPLE.INDEL.plots.R

	echo "`$NOW`INFO $SCRIPT_CODE copying INDEL recalibration reports to $RESULTS_DIR/recalibration..."
	cp $SAMPLE.INDEL.recal $RESULTS_DIR/recalibration
	cp $SAMPLE.INDEL.tranches $RESULTS_DIR/recalibration
	cp $SAMPLE.INDEL.plots.R $RESULTS_DIR/recalibration
	cp $SAMPLE.INDEL.plots.R.pdf $RESULTS_DIR/recalibration


	# step 4: recalibratedSNPs.vcf <- ApplyRecalibration(raw.snps.vcf, snp.model)
	echo "`$NOW`INFO $SCRIPT_CODE applying SNP recalibration" ## will try for the whole file, if very slow, then will break down
	# parameters were updated according to The GATK Guide Book for Version 3.1-1: increased --ts_filter_level 99.5

	java -Xmx$JAVA_XMX -XX:+UseSerialGC -Djava.io.tmpdir=$TMPDIR/tmp -jar $GATK_HOME/GenomeAnalysisTK.jar \
		-T ApplyRecalibration \
		-nt $NT \
		-R $TMPDIR/reference.fa \
		-input $TMPDIR/raw_merged.vcf \
		--ts_filter_level 99.5 \
		-mode SNP \
		-tranchesFile $TMPDIR/$SAMPLE.SNP.tranches \
		-recalFile $TMPDIR/$SAMPLE.SNP.recal \
		-o $TMPDIR/$SAMPLE.recalibratedSNPs.rawINDELs.vcf

	# step 5: analysisReady.vcf <- ApplyRecalibration(recalibratedSNPs.rawIndels.vcf, indel.model, INDEL)
	echo "`$NOW`INFO $SCRIPT_CODE applying INDEL recalibration" ## will try for the whole file, if very slow, then will break down
	java -Xmx$JAVA_XMX -XX:+UseSerialGC -Djava.io.tmpdir=$TMPDIR/tmp -jar $GATK_HOME/GenomeAnalysisTK.jar \
		-T ApplyRecalibration \
		-nt $NT \
		-R $TMPDIR/reference.fa \
		-input $TMPDIR/$SAMPLE.recalibratedSNPs.rawINDELs.vcf \
		--ts_filter_level 99.0 \
		-mode INDEL \
		-tranchesFile $TMPDIR/$SAMPLE.INDEL.tranches \
		-recalFile $TMPDIR/$SAMPLE.INDEL.recal \
		-o $TMPDIR/$SAMPLE.recalibratedSNPs.recalibratedINDELs.vcf

	# step 6: extract PASS variants
	echo "`$NOW`INFO $SCRIPT_CODE extracting PASS variants..."
	awk '{if ( ($7 == "PASS") || ($1 ~ /^#/) ) print}' $TMPDIR/$SAMPLE.recalibratedSNPs.recalibratedINDELs.vcf > $TMPDIR/$SAMPLE.recalibratedSNPs.recalibratedINDELs.PASS.vcf

#	echo "`${NOW}`INFO $SCRIPT_CODE removing variants from auxiliary samples from SNP recalibrated VCF file..."
#	java -Xmx$JAVA_XMX -Djava.io.tmpdir=$TMPDIR/tmp -jar $GATK_HOME/GenomeAnalysisTK.jar \
#		-T SelectVariants \
#		-R $REFERENCE_FASTA \
#		--variant $TMPDIR/$SAMPLE.recalibratedSNPs.rawINDELs.vcf \
#		-o $SAMPLE.recalibratedSNPs.rawINDELs.noaux.vcf \
#		--sample_file $TMPDIR/output_samples.txt \
#		--excludeNonVariants

#	echo "`$NOW`INFO $SCRIPT_CODE compressing SNP recalibrated VCF file..."
#	$TABIX_HOME/bgzip $TMPDIR/$SAMPLE.recalibratedSNPs.rawINDELs.noaux.vcf

#	echo "`$NOW`INFO $SCRIPT_CODE copying SNP recalibrated VCF file to $ANALYSIS_DIR..."
#	cp $SAMPLE.recalibratedSNPs.rawINDELs.noaux.vcf.gz $ANALYSIS_DIR/$SAMPLE.recalibratedSNPs.rawINDELs.vcf.gz

	echo "`${NOW}`INFO $SCRIPT_CODE removing variants from auxiliary samples from SNP and INDEL recalibrated VCF file..."
	java -Xmx$JAVA_XMX -XX:+UseSerialGC -Djava.io.tmpdir=$TMPDIR/tmp -jar $GATK_HOME/GenomeAnalysisTK.jar \
		-T SelectVariants \
		-R $TMPDIR/reference.fa \
		--variant $TMPDIR/$SAMPLE.recalibratedSNPs.recalibratedINDELs.vcf \
		-o $TMPDIR/$SAMPLE.recalibratedSNPs.recalibratedINDELs.noaux.vcf \
		--sample_file $TMPDIR/output_samples.txt \
		--excludeNonVariants

	echo "`$NOW`INFO $SCRIPT_CODE compressing SNP and INDEL recalibrated VCF file..."
	$TABIX_HOME/bgzip $TMPDIR/$SAMPLE.recalibratedSNPs.recalibratedINDELs.noaux.vcf

	echo "`$NOW`INFO $SCRIPT_CODE indexing SNP and INDEL recalibrated VCF file..."
	$TABIX_HOME/tabix -p vcf $TMPDIR/$SAMPLE.recalibratedSNPs.recalibratedINDELs.noaux.vcf.gz

	echo "`$NOW`INFO $SCRIPT_CODE copying SNP and INDEL recalibrated VCF file and index to $RESULTS_DIR..."
	cp $SAMPLE.recalibratedSNPs.recalibratedINDELs.noaux.vcf.gz $RESULTS_DIR/$SAMPLE.recalibrated.vcf.gz
	cp $SAMPLE.recalibratedSNPs.recalibratedINDELs.noaux.vcf.gz.tbi $RESULTS_DIR/$SAMPLE.recalibrated.vcf.gz.tbi

	echo "`$NOW`INFO $SCRIPT_CODE compressing SNP and INDEL recalibrated VCF file for all samples..."
	$TABIX_HOME/bgzip $TMPDIR/$SAMPLE.recalibratedSNPs.recalibratedINDELs.vcf

	echo "`$NOW`INFO $SCRIPT_CODE indexing SNP and INDEL recalibrated VCF file for all samples..."
	$TABIX_HOME/tabix -p vcf $TMPDIR/$SAMPLE.recalibratedSNPs.recalibratedINDELs.vcf.gz

	echo "`$NOW`INFO $SCRIPT_CODE copying SNP and INDEL recalibrated VCF file and index to $RESULTS_DIR..."
	cp $SAMPLE.recalibratedSNPs.recalibratedINDELs.vcf.gz $RESULTS_DIR/$SAMPLE.recalibrated.all_samples.vcf.gz
	cp $SAMPLE.recalibratedSNPs.recalibratedINDELs.vcf.gz.tbi $RESULTS_DIR/$SAMPLE.recalibrated.all_samples.vcf.gz.tbi

	#logging
	STATUS=OK
	if [[ ! -s $RESULTS_DIR/$SAMPLE.recalibrated.vcf.gz ]]
	then
		STATUS=FAILED
	fi

	echo -e "`${NOW}`$SCRIPT_CODE\t$SAMPLE\tall\trecalibrated_vcf\t$STATUS" >> $RUN_LOG


	echo "`${NOW}`INFO $SCRIPT_CODE removing variants from auxiliary samples from PASS variants VCF file..."
	java -Xmx$JAVA_XMX -XX:+UseSerialGC -Djava.io.tmpdir=$TMPDIR/tmp -jar $GATK_HOME/GenomeAnalysisTK.jar \
		-T SelectVariants \
		-R $TMPDIR/reference.fa \
		--variant $TMPDIR/$SAMPLE.recalibratedSNPs.recalibratedINDELs.PASS.vcf \
		-o $TMPDIR/$SAMPLE.recalibratedSNPs.recalibratedINDELs.PASS.noaux.vcf \
		--sample_file $TMPDIR/output_samples.txt \
		--excludeNonVariants

	echo "`$NOW`INFO $SCRIPT_CODE compressing PASS variants VCF file..."
	$TABIX_HOME/bgzip $TMPDIR/$SAMPLE.recalibratedSNPs.recalibratedINDELs.PASS.noaux.vcf

	echo "`$NOW`INFO $SCRIPT_CODE indexing PASS variants VCF file..."
	$TABIX_HOME/tabix -p vcf $TMPDIR/$SAMPLE.recalibratedSNPs.recalibratedINDELs.PASS.noaux.vcf.gz

	echo "`$NOW`INFO $SCRIPT_CODE copying PASS variants VCF file to $RESULTS_DIR..."
	cp $TMPDIR/$SAMPLE.recalibratedSNPs.recalibratedINDELs.PASS.noaux.vcf.gz $RESULTS_DIR/$SAMPLE.recalibrated.PASS.vcf.gz
	cp $TMPDIR/$SAMPLE.recalibratedSNPs.recalibratedINDELs.PASS.noaux.vcf.gz.tbi $RESULTS_DIR/$SAMPLE.recalibrated.PASS.vcf.gz.tbi

	echo "`$NOW`INFO $SCRIPT_CODE compressing PASS variants VCF file for all samples..."
	$TABIX_HOME/bgzip $TMPDIR/$SAMPLE.recalibratedSNPs.recalibratedINDELs.PASS.vcf

	echo "`$NOW`INFO $SCRIPT_CODE indexing PASS variants VCF file..."
	$TABIX_HOME/tabix -p vcf $TMPDIR/$SAMPLE.recalibratedSNPs.recalibratedINDELs.PASS.vcf.gz

	echo "`$NOW`INFO $SCRIPT_CODE copying PASS variants VCF file to $RESULTS_DIR..."
	cp $TMPDIR/$SAMPLE.recalibratedSNPs.recalibratedINDELs.PASS.vcf.gz $RESULTS_DIR/$SAMPLE.recalibrated.PASS.all_samples.vcf.gz
	cp $TMPDIR/$SAMPLE.recalibratedSNPs.recalibratedINDELs.PASS.vcf.gz.tbi $RESULTS_DIR/$SAMPLE.recalibrated.PASS.all_samples.vcf.gz.tbi



	#logging
	STATUS=OK
	if [[ ! -s $RESULTS_DIR/$SAMPLE.recalibrated.PASS.vcf.gz ]]
	then
		STATUS=FAILED
	fi

	echo -e "`${NOW}`$SCRIPT_CODE\t$SAMPLE\tall\tfiltered_vcf\t$STATUS" >> $RUN_LOG


	FINAL_VCF="$TMPDIR/$SAMPLE.recalibratedSNPs.recalibratedINDELs.noaux.vcf.gz"
	
	
else

	echo "`$NOW`INFO $SCRIPT_CODE Targetet sequencing data -> Skipping variant recalibration as variant set size will be to small to build recalibration model. "
	
fi


# step 7: variant evaluation with VariantEval 
## or to do separately for SNPs and INDELs?
echo "`$NOW`running variant evalutation"
FINAL_VCF_UNCOMPRESSED=`basename $FINAL_VCF .gz`
gzip -d $FINAL_VCF
#against dbSNP 129

java -Xmx$JAVA_XMX -XX:+UseSerialGC -Djava.io.tmpdir=$TMPDIR/tmp -jar $GATK_HOME/GenomeAnalysisTK.jar \
  -T VariantEval \
  -nt $NT \
  -R $TMPDIR/reference.fa \
  --eval $FINAL_VCF_UNCOMPRESSED \
  --dbsnp $DBSNP_EX_POST_129_FILENAME \
  -ST Sample \
  -ST Filter \
  -noEV \
  -EV CompOverlap \
  -EV CountVariants \
  -EV TiTvVariantEvaluator \
  -o $TMPDIR/$SAMPLE.varianteval.dbSNP129.report


echo "`$NOW`INFO $SCRIPT_CODE copying variant evaluation report dbSNP 129 to $RESULTS_DIR/recalibration..."
cp $TMPDIR/$SAMPLE.varianteval.dbSNP129.report $RESULTS_DIR/recalibration

#against current version
 
java -Xmx$JAVA_XMX -XX:+UseSerialGC -Djava.io.tmpdir=$TMPDIR/tmp -jar $GATK_HOME/GenomeAnalysisTK.jar \
  -T VariantEval \
  -nt $NT \
  -R $TMPDIR/reference.fa \
  --eval $FINAL_VCF_UNCOMPRESSED \
  --dbsnp $DBSNP_FILENAME \
  -ST Sample \
  -ST Filter \
  -noEV \
  -EV CompOverlap \
  -EV CountVariants \
  -EV TiTvVariantEvaluator \
  -o $TMPDIR/$SAMPLE.varianteval.dbSNP_current.report


echo "`$NOW`INFO $SCRIPT_CODE copying variant evaluation report dbSNP 129 to $RESULTS_DIR/recalibration..."
cp $TMPDIR/$SAMPLE.varianteval.dbSNP_current.report $RESULTS_DIR/recalibration

find $ANALYSIS_DIR -type f | xargs chmod 750 
find $RESULTS_DIR -type f | xargs chmod 750

#logging
STATUS=OK
if [[ ! -s $RESULTS_DIR/recalibration/$SAMPLE.varianteval.report ]]
then
	STATUS=FAILED
fi

echo -e "`${NOW}`$SCRIPT_CODE\t$SAMPLE\tall\tvariant_eval_report\t$STATUS" >> $RUN_LOG


echo "`$NOW`done"

#run summary script
perl $SUMMARY_SCRIPT_PATH

ls -ahl

