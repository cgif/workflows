#!/bin/bash

#PBS -l walltime=10:00:00
#PBS -l select=1:ncpus=1:mem=20gb:tmpspace=50gb

#PBS -M cgi@imperial.ac.uk
#PBS -m bea
#PBS -j oe

#PBS -q pqcgi

module load samtools/0.1.18

BAM_PATH=#bamPath
RESULTS_DIR=#resultsDir

BAM_NAME=`basename $BAM_PATH .bam`

cp $BAM_PATH $TMPDIR/tmp.bam

#keep only 'read paired' and 'read mapped in proper pair' (-f 3)
#remove 'read unmapped' (-F 4), 'mate unmapped' (-F 8), 'read is PCR or optical duplicate' (-F 1024), MAPQ < 20 (removes reads that map to multiple locations)

samtools view -b -h -f 3 -F 4 -F 8 -F 1024 -q 20 -o $TMPDIR/$BAM_NAME.filtered.bam $TMPDIR/tmp.bam 
samtools index $TMPDIR/$BAM_NAME.filtered.bam
samtools flagstat $TMPDIR/$BAM_NAME.filtered.bam > $TMPDIR/$BAM_NAME.filtered.stats 

cp $TMPDIR/$BAM_NAME.filtered* $RESULTS_DIR
chmod 660 $RESULTS_DIR/$BAM_NAME.filtered*

