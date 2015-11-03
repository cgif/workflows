#!/bin/bash

#
# runs DESeq on filtered data
#

#PBS -l walltime=#walltimeHours:00:00
#PBS -l ncpus=1
#PBS -l mem=50gb

#PBS -M cgi@imperial.ac.uk
#PBS -m bea
#PBS -j oe

module load R/3.0.1

RESULT_DIR=#resultsDir
R_SCRIPT=#rScript 

R CMD BATCH --no-save --no-restore $R_SCRIPT ${R_SCRIPT}.log

chmod 0660 $RESULT_DIR/*
