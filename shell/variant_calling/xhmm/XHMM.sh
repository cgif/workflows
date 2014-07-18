#!/bin/bash

# script to run GATK DepthOfCoverage

#PBS -l walltime=24:00:00
#PBS -l select=1:ncpus=1:mem=5gb

#PBS -M cgi@imperial.ac.uk
#PBS -m bea
#PBS -j oe
#PBS -q pqcgi

# load modules
module load xhmm/#xhmmVersion
module load R/#rVersion
module load gatk/#gatkVersion
module load java/#javaVersion

NOW="date +%Y-%m-%d%t%T%t"
JAVA_XMX=4G

DEPTHS_LIST="#depthsList"
PARAMS=#params
REFERENCE_FASTA=#referenceFasta
RESULTS_DIR=#resultsDir
R_SCRIPT=#Rscript
TARGET_INTERVALS=#targetIntervals
REFERENCE_DICT=`echo $REFERENCE_FASTA | perl -pe 's/\.fa/\.dict/'`
CHROM=`basename $RESULTS_DIR`

cp $PARAMS $TMPDIR/tmp.params
cp $REFERENCE_FASTA $TMPDIR/tmp.fasta
cp $REFERENCE_FASTA.fai $TMPDIR/tmp.fasta.fai
cp $REFERENCE_DICT $TMPDIR/tmp.dict

if [ $CHROM == "AUT" ];then

    grep -vP '^X' $TARGET_INTERVALS > $TMPDIR/tmp.intervals

elif [ $CHROM == "CHROMX_M" ];then

    grep -P '^X' $TARGET_INTERVALS > $TMPDIR/tmp.intervals

elif [ $CHROM == "CHROMX_F" ];then

    grep -P '^X' $TARGET_INTERVALS > $TMPDIR/tmp.intervals

elif [ $CHROM == ""];then

    cp $TARGET_INTERVALS > $TMPDIR/tmp.intervals

fi

#combines GATK Depth-of-Coverage outputs for multiple samples
xhmm --mergeGATKdepths \
-o  $TMPDIR/DATA.RD.txt \
$DEPTHS_LIST

cp $TMPDIR/DATA.RD.txt  $RESULTS_DIR/DATA.RD.txt 


#filter outlier samples and targets and then mean-center the targets
xhmm --matrix \
-r $TMPDIR/DATA.RD.txt \
--minMeanSampleRD 25 --maxMeanSampleRD 150 \
--maxSdSampleRD 150 \
--minMeanTargetRD 5 --maxMeanTargetRD 300 \
--centerData \
--centerType target \
--outputExcludedTargets $TMPDIR/DATA.filtered_centered.RD.txt.filtered_targets.txt \
--outputExcludedSamples $TMPDIR/DATA.filtered_centered.RD.txt.filtered_samples.txt \
-o $TMPDIR/DATA.filtered_centered.RD.txt

cp $TMPDIR/DATA.filtered_centered.RD.txt.filtered_targets.txt $RESULTS_DIR/DATA.filtered_centered.RD.txt.filtered_targets.txt
cp $TMPDIR/DATA.filtered_centered.RD.txt.filtered_samples.txt $RESULTS_DIR/DATA.filtered_centered.RD.txt.filtered_samples.txt
cp $TMPDIR/DATA.filtered_centered.RD.txt $RESULTS_DIR/DATA.filtered_centered.RD.txt


#run PCA on filtered and centered data
xhmm --PCA \
-r $TMPDIR/DATA.filtered_centered.RD.txt \
--PCAfiles $TMPDIR/DATA.RD_PCA

cp $TMPDIR/DATA.RD_PCA* $RESULTS_DIR/


#normalize mean-centered data using PCA information
xhmm --normalize \
-r $TMPDIR/DATA.filtered_centered.RD.txt \
--PCAfiles $TMPDIR/DATA.RD_PCA \
--PCnormalizeMethod PVE_mean \
--PVE_mean_factor 0.7 \
--normalizeOutput $TMPDIR/DATA.PCA_normalized.txt

cp $TMPDIR/DATA.PCA_normalized.txt $RESULTS_DIR/DATA.PCA_normalized.txt \


#filter and z-score center (by sample) the PCA-normalized data
xhmm --matrix \
-r $TMPDIR/DATA.PCA_normalized.txt \ \
--maxSdTargetRD 30 \
--centerData \
--centerType sample \
--zScoreData \
--outputExcludedTargets $TMPDIR/DATA.PCA_normalized.filtered.sample_zscores.RD.txt.filtered_targets.txt \
--outputExcludedSamples $TMPDIR/DATA.PCA_normalized.filtered.sample_zscores.RD.txt.filtered_samples.txt \
-o $TMPDIR/DATA.PCA_normalized.filtered.sample_zscores.RD.txt

cp $TMPDIR/DATA.PCA_normalized.filtered.sample_zscores.RD.txt $RESULTS_DIR/DATA.PCA_normalized.filtered.sample_zscores.RD.txt 
cp $TMPDIR/DATA.PCA_normalized.filtered.sample_zscores.RD.txt.filtered_targets.txt $RESULTS_DIR/DATA.PCA_normalized.filtered.sample_zscores.RD.txt.filtered_targets.txt
cp $TMPDIR/DATA.PCA_normalized.filtered.sample_zscores.RD.txt.filtered_samples.txt $RESULTS_DIR/DATA.PCA_normalized.filtered.sample_zscores.RD.txt.filtered_samples.txt


#filter original counts data to be the same as filtered, normalized data
xhmm --matrix \
-r $TMPDIR/DATA.RD.txt \
--excludeTargets $TMPDIR/DATA.filtered_centered.RD.txt.filtered_targets.txt \
--excludeTargets $TMPDIR/DATA.PCA_normalized.filtered.sample_zscores.RD.txt.filtered_targets.txt \
--excludeSamples $TMPDIR/DATA.filtered_centered.RD.txt.filtered_samples.txt \
--excludeSamples $TMPDIR/DATA.PCA_normalized.filtered.sample_zscores.RD.txt.filtered_samples.txt \
-o $TMPDIR/DATA.same_filtered.RD.txt

cp $TMPDIR/DATA.same_filtered.RD.txt $RESULTS_DIR/DATA.same_filtered.RD.txt


#discover CNVs in normalized data
xhmm --discover \
-p $TMPDIR/tmp.params \
-r $TMPDIR/DATA.PCA_normalized.filtered.sample_zscores.RD.txt \
-R $TMPDIR/DATA.same_filtered.RD.txt \
-c $TMPDIR/DATA.xcnv \
-a $TMPDIR/DATA.aux_xcnv \
-s $TMPDIR/DATA \
-t 60

cp $TMPDIR/DATA* $RESULTS_DIR


#genotype discovered CNVs in all samples
xhmm --genotype \
-p $TMPDIR/tmp.params \
-r $TMPDIR/DATA.PCA_normalized.filtered.sample_zscores.RD.txt \
-R $TMPDIR/DATA.same_filtered.RD.txt \
-g $TMPDIR/DATA.xcnv \
-F $TMPDIR/tmp.fasta \
-v $TMPDIR/DATA.vcf

cp $TMPDIR/DATA* $RESULTS_DIR
chmod 660 $RESULTS_DIR/*

#generate plots
mkdir $RESULTS_DIR/plots
chmod 770 $RESULTS_DIR/plots

echo -n "" > $RESULTS_DIR/DATA.locus_GC.txt
echo -n "" > $RESULTS_DIR/DATA.locus_complexity.txt

R CMD BATCH --no-save --no-restore $R_SCRIPT $R_SCRIPT.log
chmod 660 $RESULTS_DIR/plots/*
