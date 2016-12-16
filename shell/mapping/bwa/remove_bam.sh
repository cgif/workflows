#!/bin/bash
#
# script to run irods_deploy_fastq 
#

#PBS -l walltime=72:00:00
#PBS -l select=1:ncpus=1:mem=1024mb

#PBS -m ea
#PBS -M igf@imperial.ac.uk
#PBS -j oe

#PBS -q pqcgi


#now
NOW="date +%Y-%m-%d%t%T%t"

#today
TODAY=`date +%Y-%m-%d`

#set up script
SEQ_RUN_DATE=#seqRunDate
SEQ_RUN_NAME=#seqRunName
PATH_TO_DESTINATION=#pathToDestination


echo "`$NOW` removing bam files of $SEQ_RUN_DATE ..."
ssh login.cx1.hpc.ic.ac.uk "cd $PATH_TO_DESTINATION; rm -Rf $SEQ_RUN_DATE"	

echo "`$NOW` removing bam files of $SEQ_RUN_DATE completed"

