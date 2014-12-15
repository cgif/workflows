#!/bin/bash

## script to annotate with ANNOVAR

#PBS -l walltime=24:00:00
#PBS -l select=1:ncpus=1:mem=20gb

#PBS -M cgi@imperial.ac.uk
#PBS -m ea
#PBS -j oe

# load modules
ANNOVAR_PATH=#annovarPath
ANNOVAR_DB=#annovarDb
BUILD=#build

NOW="date +%Y-%m-%d%t%T%t"

# define variables

AVINPUT=#avInput
ANALYSIS_DIR=#analysisDir
RESULTS_DIR=#resultsDir
OUTPUT_PREFIX=#outputPrefix
SAMPLE=#sample
HGMD_VCF=#hgmdVcf
EXTRA_PROTOCOL="#extraProtocol"
EXTRA_OPERATION="#exraOperation"
GI_BED=#giBed
PROTOCOL="#protocol"
OPERATION="#operation"


OUTPUT=$AVINPUT.${BUILD}_multianno.txt

# annovar annotation
echo "`${NOW}`annotating calls with annovar"

#$ANNOVAR_PATH/table_annovar.pl $AVINPUT $ANNOVAR_DB/ -buildver $BUILD -protocol $PROTOCOL$EXTRA_PROTOCOL -operation $OPERATION$EXTRA_OPERATION --bedfile $GI_BED --vcfdbfile $HGMD_VCF --argument annovarArgumentextraArg -nastring NA --otherinfo --remove

$ANNOVAR_PATH/table_annovar.pl $AVINPUT $ANNOVAR_DB/ -buildver $BUILD -protocol $PROTOCOL$EXTRA_PROTOCOL -operation $OPERATION$EXTRA_OPERATION --bedfile $GI_BED --vcfdbfile $HGMD_VCF --argument annovarArgumentextraArg -nastring NA --otherinfo --remove


echo "`${NOW}`compressing output..."
gzip -f $OUTPUT

echo "`${NOW}`copying annovar annotations to $RESULTS_DIR/$SAMPLE"

#copy summary annotation to results directory
mkdir -m 770 -p $RESULTS_DIR/$SAMPLE
cp $OUTPUT.gz $RESULTS_DIR/$SAMPLE/$OUTPUT_PREFIX.$BUILD.$SAMPLE.multianno.txt.gz
chmod 660 $RESULTS_DIR/$SAMPLE/$OUTPUT_PREFIX.$BUILD.$SAMPLE.multianno.txt.gz

#extract exonic variants
zcat $RESULTS_DIR/$SAMPLE/$OUTPUT_PREFIX.$BUILD.$SAMPLE.multianno.txt.gz | head -n 1 > $RESULTS_DIR/$SAMPLE/$OUTPUT_PREFIX.$BUILD.$SAMPLE.multianno.exonic.txt
zcat $RESULTS_DIR/$SAMPLE/$OUTPUT_PREFIX.$BUILD.$SAMPLE.multianno.txt.gz | grep -E "exonic|splicing" >> $RESULTS_DIR/$SAMPLE/$OUTPUT_PREFIX.$BUILD.$SAMPLE.multianno.exonic.txt

gzip -f $RESULTS_DIR/$SAMPLE/$OUTPUT_PREFIX.$BUILD.$SAMPLE.multianno.exonic.txt

#move annovar log file
mv $AVINPUT.log $ANALYSIS_DIR/run 

#remove intermediate files
#rm $AVINPUT*

echo "`${NOW}`done"

