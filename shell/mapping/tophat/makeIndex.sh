#!/bin/bash

#
#script to create transcriptome sequence file and Bowtie index
#

#PBS -l walltime=72:00:00
#PBS -l select=1:ncpus=16:mem=25gb:tmpspace=50gb
#PBS -M igf@imperial.ac.uk
#PBS -m bea
#PBS -j oe
#PBS -q pqcgi

module load tophat/2.0.10
module load bowtie/2.1.0
module load samtools/0.1.19

USAGE='qsub -q pqcgi -v PATH_TRANSCRIPTOME_GFF="path_to_transcriptome_gff_file",PATH_REFERENCE_FASTA="path_to_bowtie_indexed_reference_fasta_file" ./index.sh'

if [ -z $PATH_TRANSCRIPTOME_GFF ] || \
    [ -z $PATH_REFERENCE_FASTA ] 
    then
    echo $USAGE
    exit 1
fi

PATH_TRANSCRIPTOME_INDEX=${PATH_TRANSCRIPTOME_GFF%.*}
PATH_REFERENCE_INDEX=${PATH_REFERENCE_FASTA%.*}

tophat -G $PATH_TRANSCRIPTOME_GFF --transcriptome-index=$PATH_TRANSCRIPTOME_INDEX $PATH_REFERENCE_INDEX
