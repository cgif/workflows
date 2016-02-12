#!/bin/bash

# script to run GATK DepthOfCoverage

#PBS -l walltime=24:00:00
#PBS -l select=1:ncpus=1:mem=5gb

#PBS -M igf@imperial.ac.uk
#PBS -m bea
#PBS -j oe
#PBS -q pqcgi

# load modules
module load gatk/#gatkVersion
module load java/#javaVersion

NOW="date +%Y-%m-%d%t%T%t"
JAVA_XMX=4G

BAM_PATH=#bamPath
TARGET_INTERVALS=#targetIntervals
REFERENCE_FASTA=#referenceFasta
RESULTS_DIR=#resultsDir
REFERENCE_DICT=`echo $REFERENCE_FASTA | perl -pe 's/\.fa/\.dict/'`
CHROM_X=#chromX

BAM_NAME=`basename ${BAM_PATH%%.*}`

cp $BAM_PATH $TMPDIR/tmp.bam
cp $BAM_PATH.bai $TMPDIR/tmp.bam.bai
cp $TARGET_INTERVALS $TMPDIR/tmp.intervals
cp $REFERENCE_FASTA $TMPDIR/tmp.fasta
cp $REFERENCE_FASTA.fai $TMPDIR/tmp.fasta.fai
cp $REFERENCE_DICT $TMPDIR/tmp.dict

#run GATK depth of coverage
java -Xmx$JAVA_XMX -Djava.io.tmpdir=$TMPDIR/tmp -jar $GATK_HOME/GenomeAnalysisTK.jar \
    -T DepthOfCoverage \
    -I $TMPDIR/tmp.bam \
    -L $TMPDIR/tmp.intervals \
    -R $TMPDIR/tmp.fasta \
    -rf BadCigar \
    -dt BY_SAMPLE \
    -dcov 5000 \
    --omitDepthOutputAtEachBase \
    --omitLocusTable \
    --minBaseQuality 0 \
    --minMappingQuality 20 \
    --start 1 \
    --stop 5000 \
    --nBins 200 \
    --includeRefNSites \
    --countType COUNT_FRAGMENTS \
    -o $TMPDIR/${BAM_NAME}.counts

cp $TMPDIR/${BAM_NAME}.counts* $RESULTS_DIR

#divide counts summary for intervals on autosomes and chromX into separate files
if [ "$CHROM_X" == "T" ];then

    grep -vP '^X' $TMPDIR/${BAM_NAME}.counts.sample_interval_summary > $RESULTS_DIR/${BAM_NAME}.autosomes_counts.sample_interval_summary
    head -n 1 $TMPDIR/${BAM_NAME}.counts.sample_interval_summary > $RESULTS_DIR/${BAM_NAME}.chromX_counts.sample_interval_summary
    grep -P '^X' $TMPDIR/${BAM_NAME}.counts.sample_interval_summary >> $RESULTS_DIR/${BAM_NAME}.chromX_counts.sample_interval_summary

fi

chmod 0660 $RESULTS_DIR/${BAM_NAME}*