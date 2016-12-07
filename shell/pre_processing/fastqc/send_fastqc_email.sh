#!/bin/bash
#
# script to run send_fastq_email 
#

#PBS -l walltime=72:00:00
#PBS -l select=1:ncpus=1:mem=1mb

#PBS -m ea
#PBS -M igf@imperial.ac.uk
#PBS -j oe

#PBS -q pqcgi


#now
NOW="date +%Y-%m-%d%t%T%t"

#today
TODAY=`date +%Y-%m-%d`

#set up script
MAIL_PATH=#mailPath

sendmail -t < $MAIL_PATH
