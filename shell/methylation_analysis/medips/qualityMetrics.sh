#!/bin/bash

# runs MEDIPS QC on set of BAM files

#PBS -l walltime=72:00:00
#PBS -l ncpus=1
#PBS -l mem=20gb

#PBS -M igf@imperial.ac.uk
#PBS -m bea
#PBS -j oe

module load R/#Rversion

R_SCRIPT="#Rscript"
BAM_FILE="#inputBam"

BAM_TMP=$TMPDIR/tmp.bam
cp $BAM_FILE $BAM_TMP
sed -i -e "s/#bamFile/${BAM_TMP//\//\\/}/" $R_SCRIPT

R --vanilla < $R_SCRIPT
