#!/bin/bash

#PBS -l walltime=#walltimeHours:00:00
#PBS -l ncpus=#threadsPerRun
#PBS -l mem=50g
#PBS -M igf@imperial.ac.uk
#PBS -m bea
#PBS -j oe

THREADS_PER_RUN=#threadsPerRun

BOWTIE_VERSION=#bowtieVersion
BLAST_VERSION=#blastVersion
BLAT_VERSION=#blatVersion

BASEDIR=#baseDir
ANALYSIS=#analysis
RESULTS_DIR=#resultPath
READ_PATH=#readPath
ACTIONS_FILE=#actionsFile

BOWTIE_HG_DATABASE=#bowtieHGdatabase
BLAT_HG_DATABASE=#blatHGdatabase
BLAT_HG_OOC_FILE=#blatHGOOCfile
BLAST_NT_DATABASE=#blastNTdatabase
BLAST_HG_DATABASE=#blastHGdatabase
TAXONOMY_DMP=#taxonomyDmp

#load modules
module load bowtie/$BOWTIE_VERSION
module load blast+/$BLAST_VERSION
module load blat/$BLAT_VERSION

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
sed  -i -e "s/#TMPDIR/${TMPDIR//\//\\/}/" $TMPDIR/$ANALYSIS/createLogTable.py
sed  -i -e "s/#ANALYSIS/${ANALYSIS//\//\\/}/g" $TMPDIR/$ANALYSIS/createLogTable.py

#replace empty and sequence lines (together with corresponding quality lines) created by cutadapt
#unzip and copy fastqs to tmp space
echo "`${NOW}`copying reads $READ_PATH to temporary scratch space..."

#generate string that will replace empty and short sequence lines in fastq files
READ_LENGTH=`gzip -d -c $READ_PATH | head -n 100 | awk '{if(NR%4==2) print length($1)}' | sort -n| uniq | tail -n 1`
STRING=$(for i in `eval echo {1..$READ_LENGTH}`;do printf "%s" "N";done;)

FASTQ_READ_NO_EXT=`basename $READ_PATH .gz`
echo "`${NOW}`$FASTQ_READ_NO_EXT"
gzip -c -d $READ_PATH | sed "s/^$/$STRING/g" | sed "s/\s/#0\//g" > $TMPDIR/$FASTQ_READ_NO_EXT

#create pipelineScript.py
python $TMPDIR/$ANALYSIS/master.py -i $TMPDIR/$FASTQ_READ_NO_EXT -d "#0/" -p $ACTIONS_FILE -a $THREADS_PER_RUN -o 33

#run pipelineScript.py to filter human sequences and blast non-human reads against nt database
python $TMPDIR/pipelineScript.py

#create log table with info for each action
python $TMPDIR/$ANALYSIS/createLogTable.py -i $TMPDIR/log* > $RESULTS_DIR/logTable.txt

FASTQ_READ_REMOVE_EXT=${FASTQ_READ_NO_EXT%.*}
cp $TMPDIR/log* $RESULTS_DIR
cp $TMPDIR/$FASTQ_READ_REMOVE_EXT* $RESULTS_DIR
chmod 660 $RESULTS_DIR/*

ls -l


