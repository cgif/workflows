#!/bin/bash

#PBS -l walltime=1:00:00
#PBS -l mem=50g
#PBS -M igf@imperial.ac.uk
#PBS -m bea
#PBS -j oe

BOWTIE_VERSION=#bowtieVersion
BLAST_VERSION=#blastVersion
BLAT_VERSION=#blatVersion

BASEDIR=#baseDir
ANALYSIS=#analysis
RESULTS_DIR=#resultPath

BOWTIE_HG_DATABASE=#bowtieHGdatabase
BLAT_HG_DATABASE=#blatHGdatabase
BLAT_HG_OOC_FILE=#blatHGOOCfile
BLAST_NT_DATABASE=#blastNTdatabase
BLAST_HG_DATABASE=#blastHGdatabase
TAXONOMY_DMP=#taxonomyDmp

module load bowtie/$BOWTIE_VERSION
BLAST_HOME=/apps/ncbi-blast/$BLAST_VERSION/bin
BLAT_HOME=/apps/blat/$BLAT_VERSION/bin

#copy imsa scripts to tmp space
mkdir $TMPDIR/$ANALYSIS
cp $BASEDIR/*.py $TMPDIR/$ANALYSIS

#send parameters to configuration file
sed  -i -e "s/#TMPDIR/${TMPDIR//\//\\/}/" $TMPDIR/$ANALYSIS/config.py
sed  -i -e "s/#BOWTIE_HG_DATABASE_TMP/${BOWTIE_HG_DATABASE//\//\\/}/" $TMPDIR/$ANALYSIS/config.py
sed  -i -e "s/#BLAT_HG_DATABASE_TMP/${BLAT_HG_DATABASE//\//\\/}/" $TMPDIR/$ANALYSIS/config.py
sed  -i -e "s/#BLAT_HG_OOC_FILE_TMP/${BLAT_HG_OOC_FILE//\//\\/}/" $TMPDIR/$ANALYSIS/config.py
sed  -i -e "s/#BLAST_NT_DATABASE_TMP/${BLAST_NT_DATABASE//\//\\/}/" $TMPDIR/$ANALYSIS/config.py
sed  -i -e "s/#BLAST_HG_DATABASE_TMP/${BLAST_HG_DATABASE//\//\\/}/" $TMPDIR/$ANALYSIS/config.py
sed  -i -e "s/#TAXONOMY_DMP_TMP/${TAXONOMY_DMP//\//\\/}/" $TMPDIR/$ANALYSIS/config.py
sed  -i -e "s/#BOWTIE_HOME/${BOWTIE_HOME//\//\\/}/" $TMPDIR/$ANALYSIS/config.py
sed  -i -e "s/#BLAST_HOME/${BLAST_HOME//\//\\/}/" $TMPDIR/$ANALYSIS/config.py
sed  -i -e "s/#BLAT_HOME/${BLAT_HOME//\//\\/}/" $TMPDIR/$ANALYSIS/config.py
sed  -i -e "s/#ANALYSIS/${ANALYSIS//\//\\/}/g" $TMPDIR/$ANALYSIS/config.py

#combine blast results from two reads
READ_PATH1=#readPath1
READ_PATH2=#readPath2
READ1=`basename $READ_PATH1`
READ2=`basename $READ_PATH2`
READ1=${READ1%%.*}
READ2=${READ2%%.*}
RESULTS_FILE1_BASENAME=`basename $(ls --color=never $RESULTS_DIR/*fa | grep $READ1 | awk '{ print length(), $0}' | sort -n | tail -n 1 | cut -f 2 -d " ") .fa`
RESULTS_FILE2_BASENAME=`basename $(ls --color=never $RESULTS_DIR/*fa | grep $READ2 | awk '{ print length(), $0}' | sort -n | tail -n 1 | cut -f 2 -d " ") .fa`
cat $RESULTS_DIR/$RESULTS_FILE1_BASENAME.bln $RESULTS_DIR/$RESULTS_FILE2_BASENAME.bln > $TMPDIR/combined.$RESULTS_FILE1_BASENAME.bln

#create tax files
python $TMPDIR/$ANALYSIS/postprocess.py -b $TMPDIR/combined.$RESULTS_FILE1_BASENAME.bln
cp $TMPDIR/combined* $RESULTS_DIR
ls -l 

chmod 770 $RESULTS_DIR
chmod 660 $RESULTS_DIR/*
