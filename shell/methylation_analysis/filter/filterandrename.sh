#!/bin/bash

# runs filtering of BAM files

#PBS -l walltime=72:00:00
#PBS -l ncpus=2
#PBS -l mem=100gb

#PBS -M cgi@imperial.ac.uk
#PBS -m bea
#PBS -j oe

module load samtools/#samtoolsVersion

MERGETAG_DIR="#mergetagDir"
FILTERED_DIR="#filteredDir"
SAMPLE="#sample"

cp $MERGETAG_DIR/$SAMPLE/$SAMPLE.bam $TMPDIR/$SAMPLE.bam

echo "removing duplicates and unmapped reads and renaming chromosome names to hg19 from sample $SAMPLE..."

samtools view -hF 1028 $TMPDIR/$SAMPLE.bam \
| perl -e 'open(CHROM,"/groupvol/cgi/resources/reference/eukaryote/human/hs37d5_to_hg19.chrom_names.txt"); %rename = (); while (<CHROM>){ /(.*)\t(.*)/; $rename{$1} = $2;} while(<>) {if (/^@/ && /SN:/) { s/(.*\tSN:)(.*)(\t.*)/$1$rename{$2}$3/; print } elsif (/^@/) { print } else { @cols = split(/\t/); $cols[2] =~ s/(.*)/$rename{$1}/; $cols[6] =~ s/(.*)/$rename{$1}/ unless $cols[6] eq "="; $line = ""; foreach $col (@cols) { $line .= $col."\t"; } chop($line); print "$line" } }' \
| samtools view -bS - > $TMPDIR/$SAMPLE.nondup.rename.filt.bam

samtools index $TMPDIR/$SAMPLE.nondup.rename.filt.bam

cp $TMPDIR/$SAMPLE.nondup.rename.filt.bam $FILTERED_DIR/$SAMPLE/$SAMPLE.nondup.rename.filt.bam
cp $TMPDIR/$SAMPLE.nondup.rename.filt.bam.bai $FILTERED_DIR/$SAMPLE/$SAMPLE.nondup.rename.filt.bam.bai





