#!/bin/bash

## run miRDeep2 mapper.pl

#PBS -l walltime=24:00:00
#PBS -l select=1:ncpus=1:mem=5gb

#PBS -M igf@imperial.ac.uk
#PBS -m ae
#PBS -j oe

#PBS -q pqcgi

NOW="date +%Y-%m-%d%t%T%t"
TODAY=`date +%Y-%m-%d`

#load modules
module load mirdeep2/0.0.7 
module load bowtie/1.1.1
module load vienna-rna/2.1.8
module load squid/1.9g
module load randfold/2.0

#define variables
REFERENCE=#reference
INDEXED_PATH=#indexedPath
ADAPTER=#adapter
INPUT_FASTQ=#inputFastq
SAMPLE=#sample
OUTPUT_DIR=#outputDir
DATE=#today

#variables for quantifier

SPECIES=#species
MIRBASE_DIR=#mirbaseDir
PRECURSORS_ALL=$MIRBASE_DIR/hairpin.fa
MATURE_ALL=$MIRBASE_DIR/mature.fa

#deployment
DEPLOYMENT_SERVER=#deploymentServer
DEPLOYMENT_PATH=#deploymentPath

echo "`${NOW}`copying reference to temporary space"
cp $REFERENCE $TMPDIR
REFERENCE_PREFIX=`basename $REFERENCE .fa`

echo "`${NOW}`copying bowtie index reference to temporary space"
cp $INDEXED_PATH/* $TMPDIR

echo "`${NOW}`copying input to temporary space"
cp $INPUT_FASTQ $TMPDIR/fastq.fq.gz
gunzip $TMPDIR/fastq.fq.gz

echo "`${NOW}`running mapper.pl"
mapper.pl fastq.fq -e -h -k $ADAPTER -l 18 -m -p $REFERENCE_PREFIX \
			-s $SAMPLE.reads_collapsed.fa \
			-t $SAMPLE.reads_collapsed_vs_GRCh37.arf \
			-v

echo "`${NOW}`copying output to $OUTPUT_DIR/mapper"
cp $SAMPLE.reads_collapsed.fa $OUTPUT_DIR/mapper
cp $SAMPLE.reads_collapsed_vs_GRCh37.arf $OUTPUT_DIR/mapper
cp *log $OUTPUT_DIR/mapper
chmod 660 $OUTPUT_DIR/mapper/*

echo "`${NOW}`list temporary folder"
ls -al $TMPDIR

echo "`${NOW}`copying miRBase databases to temporary space"

cp $PRECURSORS_ALL $TMPDIR/hairpin.fa
cp $MATURE_ALL $TMPDIR/mature.fa

echo "`${NOW}`running quantifier.pl"
quantifier.pl -p hairpin.fa -m mature.fa -r $SAMPLE.reads_collapsed.fa -t $SPECIES -P -y $DATE

#make expression data file for edgeR
awk -v OFS='\t' '{ if ($2 != "0") print $3 "," $1, $2 }' $TMPDIR/expression_analyses/expression_analyses_$DATE/miRNA_expressed.csv | tail -n +2 | sort -k 2 -rn | uniq > $TMPDIR/expression_analyses/expression_analyses_$DATE/miRNA_expressed.tsv


echo "`${NOW}`copying output to $OUTPUT_DIR/quantifier"
cp $TMPDIR/expression_$DATE.html $OUTPUT_DIR/quantifier
cp $TMPDIR/miRNAs_expressed_all_samples_$DATE.csv $OUTPUT_DIR/quantifier
chmod 660 $OUTPUT_DIR/quantifier/*
cp -R $TMPDIR/expression_analyses/expression_analyses_$DATE $OUTPUT_DIR/quantifier
cp -R $TMPDIR/pdfs_$DATE $OUTPUT_DIR/quantifier
chmod 770 $OUTPUT_DIR/quantifier/expression_analyses_$DATE
chmod 660 $OUTPUT_DIR/quantifier/expression_analyses_$DATE/*
chmod 770 $OUTPUT_DIR/quantifier/pdfs_$DATE
chmod 660 $OUTPUT_DIR/quantifier/pdfs_$DATE/*

echo "`${NOW}`creating deployment directory $DEPLOYMENT_PATH on $DEPLOYMENT_SERVER..."
ssh $DEPLOYMENT_SERVER "mkdir -p -m 775 $DEPLOYMENT_PATH" > /dev/null 2>&1

echo "`${NOW}`deploying quantifier results to $DEPLOYMENT_SERVER:$DEPLOYMENT_PATH..."
scp $TMPDIR/expression_$DATE.html $DEPLOYMENT_SERVER:$DEPLOYMENT_PATH/$SAMPLE.expression.html > /dev/null 2>&1
scp $TMPDIR/expression_analyses/expression_analyses_$DATE/miRNA_expressed.tsv $DEPLOYMENT_SERVER:$DEPLOYMENT_PATH/$SAMPLE.miRNA_expressed.tsv  > /dev/null 2>&1
scp -r $TMPDIR/pdfs_$DATE $DEPLOYMENT_SERVER:$DEPLOYMENT_PATH/  > /dev/null 2>&1

ssh $DEPLOYMENT_SERVER "chmod 664 $DEPLOYMENT_PATH/$SAMPLE.expression.html" > /dev/null 2>&1
ssh $DEPLOYMENT_SERVER "chmod 664 $DEPLOYMENT_PATH/$SAMPLE.miRNA_expressed.tsv" > /dev/null 2>&1
ssh $DEPLOYMENT_SERVER "chmod 775 $DEPLOYMENT_PATH/pdfs_$DATE" > /dev/null 2>&1
ssh $DEPLOYMENT_SERVER "chmod 664 $DEPLOYMENT_PATH/pdfs_$DATE/*" > /dev/null 2>&1

echo "`${NOW}`list temporary folder"
ls -al $TMPDIR

echo "`${NOW}`Done"

