#!/bin/bash

## script to run CREST

#PBS -l walltime=72:00:00
#PBS -l select=1:ncpus=1:mem=2gb

#PBS -m ea
#PBS -M cgi@imperial.ac.uk
#PBS -j oe

# load modules
module load crest/#crestVersion
module load samtools/#samtoolsVersion
module load blat/#blatVersion

#now
NOW="date +%Y-%m-%d%t%T%t"

# define variables
REFERENCE_FASTA=#referenceFasta
REFERENCE_2BIT=#reference2Bit

BAM_NAME=#sample
TUMOR_BAM=#tumorBam
GERMLINE_BAM=#germlineBam
ANALYSIS_PATH=#analysisPath
PREF=`basename $TUMOR_BAM`

READ_LENGTH=#readLength
MIN_SC_READS=#minScReads
SENSITIVE=#sensitive
PORT=#serverPort

cp $REFERENCE_FASTA $TMPDIR/tmp.fasta
cp $REFERENCE_FASTA.fai $TMPDIR/tmp.fasta.fai
cp $REFERENCE_2BIT $TMPDIR/tmp.2bit

cp $TUMOR_BAM $TMPDIR/tumor.bam
cp $TUMOR_BAM.bai $TMPDIR/tumor.bam.bai

cp $GERMLINE_BAM $TMPDIR/germline.bam
cp $GERMLINE_BAM.bai $TMPDIR/germline.bam.bai

######################################
#step 1. Get soft-clipping positions.
echo "`${NOW}`getting soft-clipping positions"
echo "`${NOW}`file $ANALYSIS_PATH/crestSM/$PREF.cover will be generated"

extractSClip.pl -i $TMPDIR/tumor.bam --ref_genome $TMPDIR/tmp.fasta 

# change permissions of the results files
cp $TMPDIR/tumor.bam.cover $ANALYSIS_PATH/crestSM/$PREF.cover
chmod 0660 $ANALYSIS_PATH/crestSM/$PREF.cover

##################################
# step 2: Run SV detection script.
echo "`${NOW}`running SV detection script"
echo "`${NOW}`file $ANALYSIS_PATH/crestSM/$PREF.predSV.txt will be generated"

gfServer start localhost $PORT $TMPDIR/tmp.2bit &

STAT=0
while [ $STAT -lt 2 ]
do
    sleep 20
    STAT=`ps -elf|grep "gfServer start localhost $PORT"|grep " S "|wc -l`
done 

if [[ "$SENSITIVE" == "TRUE" ]]; then
    CREST.pl -f $TMPDIR/tumor.bam.cover -d $TMPDIR/tumor.bam -g $TMPDIR/germline.bam --ref_genome $TMPDIR/tmp.fasta -t $TMPDIR/tmp.2bit --blatserver localhost --blatport $PORT -l $READ_LENGTH --min_sclip_reads $MIN_SC_READS --sensitive 
else
    CREST.pl -f $TMPDIR/tumor.bam.cover -d $TMPDIR/tumor.bam -g $TMPDIR/germline.bam --ref_genome $TMPDIR/tmp.fasta -t $TMPDIR/tmp.2bit --blatserver localhost --blatport $PORT -l $READ_LENGTH --min_sclip_reads $MIN_SC_READS
fi

# change permissions of the results files
cp $TMPDIR/tumor.bam.predSV.txt $ANALYSIS_PATH/crestSM/$PREF.predSV.txt
chmod 0660 $ANALYSIS_PATH/crestSM/$PREF.predSV.txt

# kill gfserver
gfServer stop localhost $PORT

##################################
# step 3: Create html file to display alignment at breakpoint.
echo "`${NOW}`creating breakpoint html file"
echo "`${NOW}`file $ANALYSIS_PATH/crestSM/$PREF.predSV.html will be generated"

bam2html.pl -d $TMPDIR/tumor.bam -g $TMPDIR/germline.bam -i $TMPDIR/tumor.bam.predSV.txt --ref_genome $TMPDIR/tmp.fasta -o $TMPDIR/tumor.bam.predSV.html
ls -l

# change permissions of the results files
cp $TMPDIR/tumor.bam.predSV.html $ANALYSIS_PATH/crestSM/$PREF.predSV.html
chmod 0660 $ANALYSIS_PATH/crestSM/$PREF.predSV.html

