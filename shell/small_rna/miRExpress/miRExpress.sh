#!/bin/bash

## script run miRExpress to get known miRNA expression counts

#PBS -l walltime=10:00:00
#PBS -l select=1:ncpus=#ncpus:mem=2gb

#PBS -M cgi@imperial.ac.uk
#PBS -m ae
#PBS -j oe

#PBS -q pqcgi

NOW="date +%Y-%m-%d%t%T%t"
TODAY=`date +%Y-%m-%d`

module load mirexpress/2.1.4
module load R/3.1.0

INPUT_FASTQ=#inputFastq
OUTPUT_DIR=#outputDir
SAMPLE=#sample
ADAPTOR=#adaptor		
MIRBASE_DIR=#mirbaseDir
NCPUS=#ncpus
R_SCRIPT=#Rscript

echo "`${NOW}`copying fasq file $INPUT_FASTQ to temporary scratch space..."
cp $INPUT_FASTQ $TMPDIR/reads.fq.gz
gunzip $TMPDIR/reads.fq.gz



echo "`${NOW}`copying miRBase databases to temporary scratch space..."
cp $MIRBASE_DIR/hsa_precursor.txt $TMPDIR
cp $MIRBASE_DIR/hsa_miRNA.txt $TMPDIR
cp $MIRBASE_DIR/AllPrecursors_Structures.txt $TMPDIR

echo "`${NOW}`Parsing Raw data..."
Raw_data_parse -i $TMPDIR/reads.fq -o $SAMPLE.stp1

# we do adaptor trimming with cutadapt, which gives better distribution
# mirExpress removes all reads where adaptor is in the middle, 
# which is not appropriate for 50bp reads

#echo "`${NOW}`copying adaptor file $ADAPTOR to temporary scratch space..."
#cp $ADAPTOR $TMPDIR/adaptor.txt

#echo "`${NOW}`Trimming Adapter..."
#Trim_adapter -i $SAMPLE.stp1 -t adaptor.txt -o $SAMPLE.stp2

echo "`${NOW}`Generating statistics of reads..."
statistics_reads -i $SAMPLE.stp1 -o $SAMPLE.read_statistics.txt

echo "`${NOW}`Create output folder \"$TMPDIR/mirExpress_results\"..."
mkdir $TMPDIR/mirExpress_results

echo "`${NOW}`Aligning reads to precursor miRNA..."
alignmentSIMD -r $TMPDIR/hsa_precursor.txt -i $SAMPLE.stp1 -o $TMPDIR/mirExpress_results/ -u $NCPUS

echo "`${NOW}`Producing results..."
analysis -r $TMPDIR/hsa_precursor.txt -m $TMPDIR/hsa_miRNA.txt -d $TMPDIR/mirExpress_results/ -o $SAMPLE.read_align_premiRNA.txt -t $SAMPLE.miRNA_expression.txt -s $TMPDIR/AllPrecursors_Structures.txt

echo "`${NOW}`copying results to $OUTPUT_DIR..."
cp ${SAMPLE}* $OUTPUT_DIR
rm -R $TMPDIR/mirExpress_results/Temp
cp -R $TMPDIR/mirExpress_results $OUTPUT_DIR

echo "`${NOW}`running summary statistics script..."

R CMD BATCH --no-save --no-restore $R_SCRIPT ${R_SCRIPT}.log

chmod 660 $OUTPUT_DIR/*
chmod 770 $OUTPUT_DIR/mirExpress_results
chmod 660 $OUTPUT_DIR/mirExpress_results/*

echo "`${NOW}`listing $TMPDIR..."
ls -al $TMPDIR
ls -al $TMPDIR/mirExpress_results

echo "`${NOW}`Done"

