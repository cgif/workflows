#!/bin/bash

#
# script run cutadapt algorithm for removing adapters and low quality ends
#

#PBS -l walltime=walltimeHours:00:00
#PBS -l select=1:ncpus=1:mem=4gb

#PBS -m bea
#PBS -M cgi@imperial.ac.uk
#PBS -j oe

#PBS -q pqcgi

module load cutadapt/1.0

#now
NOW="date +%Y-%m-%d%t%T%t"

#today
TODAY=`date +%Y-%m-%d`

#path to fastq file 
PATH_READS_FASTQ=pathReadsFastq

#path to trimmed fastq file
PATH_TRIMMED_DIR=pathTrimmedDir

#cutadapt parameters
ADAPTER5=adapter5
ADAPTER3=adapter3
OVERLAP=overlap
CUTOFF=cutoff

FASTQ_NAME=`basename $PATH_READS_FASTQ`
FASTQ_BASENAME=`echo $FASTQ_NAME | perl -pe 's/(.*)\.f.*q.*/$1/g'`

#copy files to temp directory
echo "`${NOW}`copying read $FASTQ_NAME to temporary scratch space..."
cp $PATH_READS_FASTQ $TMPDIR/$FASTQ_NAME

#make ADAPTERS variable
ADAPTER5_REV=`echo $ADAPTER5 | perl -pe 'tr/ATCGatcg/TAGCtagc/'|rev`
ADAPTER3_REV=`echo $ADAPTER3 | perl -pe 'tr/ATCGatcg/TAGCtagc/'|rev`
ADAPTERS="-g $ADAPTER5 -a $ADAPTER3 -g $ADAPTER3_REV -a $ADAPTER5_REV"

if [[ $ADAPTER5 == "none" ]]; then
	ADAPTERS="-a $ADAPTER3"
fi

echo "`${NOW}`running cutadapt"
cutadapt $ADAPTERS -O $OVERLAP -o $TMPDIR/${FASTQ_BASENAME}_trimmed.fq.gz -q $CUTOFF $TMPDIR/$FASTQ_NAME

echo "`${NOW}`copying trimmed fastqc to $PATH_TRIMMED_DIR"
cp $TMPDIR/${FASTQ_BASENAME}_trimmed.fq.gz $PATH_TRIMMED_DIR
chmod 660 $PATH_TRIMMED_DIR/*

