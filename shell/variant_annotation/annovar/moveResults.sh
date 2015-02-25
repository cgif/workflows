#!/bin/bash

## script to move ANNOVAR results and scripts to /project/tgu

#PBS -l walltime=24:00:00
#PBS -l select=1:ncpus=1:mem=10gb

#PBS -M cgi@imperial.ac.uk
#PBS -m ea
#PBS -j oe

NOW="date +%Y-%m-%d%t%T%t"

# define variables

RUN_DIR=#runDir
RESULTS_DIR=#resultsDir
PROJECT=#project
TODAY=#today

RESULTS_TGU_DIR=/project/tgu/results/$PROJECT/annovar
RUN_TGU_DIR=/project/tgu/runs/$PROJECT/annovar/$TODAY

mkdir -p $RESULTS_TGU_DIR
mkdir -p $RUN_TGU_DIR

# annovar annotation
echo "`${NOW}`moving results..."
mv $RESULTS_DIR $RESULTS_TGU_DIR
chmod -R 770 /project/tgu/results/$PROJECT/annovar

echo "`${NOW}`moving scripts..."
mv $RUN_DIR $RUN_TGU_DIR
chmod -R 770 /project/tgu/runs/$PROJECT/annovar

echo "`${NOW}`done"

