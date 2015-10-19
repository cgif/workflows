#!/bin/bash

#PBS -l walltime=24:00:00
#PBS -l select=1:ncpus=1:mem=10gb

#PBS -M cgi@imperial.ac.uk
#PBS -m bea
#PBS -j oe

#PBS -q pqcgi

module load verifybamid/#verifyBamIDVersion
module load samtools/#samtoolsVersion

BAM=#bamPath
RESULTS_DIR=#resultsDir
VCF=#vcfFile
TARGETED=#targeted
SAMPLE=`basename $RESULTS_DIR`

NOW="date +%Y-%m-%d%t%T%t"

echo "`${NOW}`copying VCF file to tmp directory..."
cp $VCF $TMPDIR/input.vcf

echo "`${NOW}`copying BAM and index file to tmp directory..."
cp $BAM $TMPDIR/input.bam
cp $BAM.bai $TMPDIR/input.bam.bai

echo "`${NOW}`sample 10% of the bam file..."
samtools view -s 0.1 -b $TMPDIR/input.bam > $TMPDIR/input.reduced.bam

echo "`${NOW}`sort and index sampled reads..."
samtools index $TMPDIR/input.reduced.bam

if [ $TARGETED == TRUE ]; then
	TARGETED_OPTIONS="--maxDepth 1000 --precise"
fi

echo "`${NOW}`running verifyBamID"
verifyBamID --vcf $TMPDIR/input.vcf \
	--bam $TMPDIR/input.reduced.bam \
	--out $SAMPLE \
	--ignoreRG \
	--chip-none \
	$TARGETED_OPTIONS

echo "`${NOW}`copying output to results directory"
cp $TMPDIR/${SAMPLE}.selfSM $RESULTS_DIR/
chmod 660 $RESULTS_DIR/${SAMPLE}.selfSM

ls -lh 

