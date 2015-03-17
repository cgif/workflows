#!/bin/bash

## script to run BICseq

#PBS -l walltime=72:00:00
#PBS -l select=1:ncpus=1:mem=10gb:tmpspace=#tmpSpacegb

#PBS -M cgi@imperial.ac.uk
#PBS -m ea
#PBS -j oe

#PBS -q pqcgi

module load samtools/#samtoolsVersion
module load R/#Rversion

BIC_SCRIPT=#BICscript
LAMBDA=#lambda
BIN_SIZE=#binSize
MULTIPLICITY=#multiplicity
WINDOW=#window
SINGLE_READ=#singleRead

RESULTS_DIR=#resultsDir
PREFIX=#prefix

TUMOR_BAM=#tumorBam
cp $TUMOR_BAM $TMPDIR/tumor.bam

NORMAL_BAM=#normalBam
cp $NORMAL_BAM $TMPDIR/normal.bam

#filter reads with non-uniq mapping, extract alignment locations, create file with locations for each chromosome
samtools view $TMPDIR/tumor.bam|awk '!($5 == 0) {print $3"\t"$4;}' > $TMPDIR/tumor.seq

#create configuration file with chromosome name, tumor and normal files with alignment locations
printf "chrom\tcase\tcontrol\n" > $TMPDIR/ConfigFile

for CHROM in `samtools view -H $TMPDIR/tumor.bam|grep -P '^@SQ'|grep -vP 'MT|GL|NC_007605|hs37d5'|cut -f 2|cut -f 2 -d :`; do

	cat $TMPDIR/tumor.seq|grep -P "^"$CHROM"\t"|cut -f 2 > $TMPDIR/tumor.$CHROM.seq

	if [[ -f $TMPDIR/tumor.$CHROM.seq ]] && [[ -s $TMPDIR/tumor.$CHROM.seq ]]; then
		printf "$CHROM\t$TMPDIR/tumor.$CHROM.seq\t$TMPDIR/normal.$CHROM.seq\n" >> $TMPDIR/ConfigFile
	fi

done

more $TMPDIR/ConfigFile
rm $TMPDIR/tumor.bam
rm $TMPDIR/tumor.seq

samtools view $TMPDIR/normal.bam|awk '!($5 == 0) {print $3"\t"$4;}' > $TMPDIR/normal.seq

for CHROM in `cut -f 1 $TMPDIR/ConfigFile | grep -v chrom`; do

	cat $TMPDIR/normal.seq|grep -P "^"$CHROM"\t"|cut -f 2 > $TMPDIR/normal.$CHROM.seq

done

rm $TMPDIR/normal.bam
rm $TMPDIR/normal.seq

if [ $SINGLE_READ == "F" ]; then

	perl $BIC_SCRIPT --lambda=$LAMBDA --bin_size=$BIN_SIZE --multiplicity=$MULTIPLICITY --window=$WINDOW --paired $TMPDIR/ConfigFile $TMPDIR/BICseq $PREFIX

else

	perl $BIC_SCRIPT --lambda=$LAMBDA --bin_size=$BIN_SIZE --multiplicity=$MULTIPLICITY --window=$WINDOW $TMPDIR/ConfigFile $TMPDIR/BICseq $PREFIX

fi 

cp $TMPDIR/BICseq/* $RESULTS_DIR
chmod 0660 $RESULTS_DIR/*



