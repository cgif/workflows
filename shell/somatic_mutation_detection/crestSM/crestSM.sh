#!/bin/bash

## script to run CREST

#PBS -l walltime=72:00:00
#PBS -l select=1:ncpus=1:mem=10gb

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
CHROM=#chrom
TUMOR_BAM=#tumorBam
GERMLINE_BAM=#germlineBam
ANALYSIS_PATH=#analysisPath
CHROM=#chrom
TUMOR_NAME=`basename $TUMOR_BAM .bam`
TUMOR_PREF=$TUMOR_NAME.$CHROM
GERMLINE_NAME=`basename $GERMLINE_BAM .bam`
GERMLINE_PREF=$GERMLINE_NAME.$CHROM

READ_LENGTH=#readLength
MIN_SC_READS=#minScReads
SENSITIVE=#sensitive
PORT=#serverPort

######################################
#step 1. Get soft-clipping positions.
echo "`${NOW}`getting soft-clipping positions"
echo "`${NOW}`files $ANALYSIS_PATH/crestSM/$TUMOR_PREF.cover and $ANALYSIS_PATH/crestSM/$GERMLINE_PREF.cover will be generated"

extractSClip.pl -i $TUMOR_BAM --ref_genome $REFERENCE_FASTA -r $CHROM -o $ANALYSIS_PATH/crestSM -p $TUMOR_NAME
extractSClip.pl -i $GERMLINE_BAM --ref_genome $REFERENCE_FASTA -r $CHROM -o $ANALYSIS_PATH/crestSM -p $GERMLINE_NAME

######################################
#step 2. Remove germline events.
echo "`${NOW}`removing germline events"
echo "`${NOW}`file $ANALYSIS_PATH/crestSM/$TUMOR_PREF.cover.somatic.cover will be generated"

$ANALYSIS_PATH/run/countDiff.pl -d $ANALYSIS_PATH/crestSM/$TUMOR_PREF.cover -g $ANALYSIS_PATH/crestSM/$GERMLINE_PREF.cover > $ANALYSIS_PATH/crestSM/$TUMOR_PREF.soft_clip.dist.txt

##################################
# step 3: Run SV detection script.
echo "`${NOW}`running SV detection script"
echo "`${NOW}`file $ANALYSIS_PATH/crestSM/$TUMOR_PREF.predSV.txt will be generated"

gfServer start localhost $PORT $REFERENCE_2BIT &

STAT=0
while [ $STAT -lt 2 ]
do
    sleep 20
    STAT=`ps -elf|grep "gfServer start localhost $PORT"|grep " S "|wc -l`
done 

if [[ "$SENSITIVE" == "TRUE" ]]; then
    CREST.pl -f $ANALYSIS_PATH/crestSM/$TUMOR_PREF.cover.somatic.cover -d $TUMOR_BAM -g $GERMLINE_BAM --ref_genome $REFERENCE_FASTA -r $CHROM -t $REFERENCE_2BIT --blatserver localhost --blatport $PORT -l $READ_LENGTH --min_sclip_reads $MIN_SC_READS --sensitive -o $ANALYSIS_PATH/crestSM -p $TUMOR_PREF
else
    CREST.pl -f $ANALYSIS_PATH/crestSM/$TUMOR_PREF.cover.somatic.cover -d $TUMOR_BAM -g $GERMLINE_BAM --ref_genome $REFERENCE_FASTA -r $CHROM -t $REFERENCE_2BIT --blatserver localhost --blatport $PORT -l $READ_LENGTH --min_sclip_reads $MIN_SC_READS -o $ANALYSIS_PATH/crestSM -p $TUMOR_PREF
fi

# kill gfserver
gfServer stop localhost $PORT

##################################
# step 4: Create html file to display alignment at breakpoint.
echo "`${NOW}`creating breakpoint html file"
echo "`${NOW}`file $ANALYSIS_PATH/crestSM/$TUMOR_PREF.predSV.html will be generated"

bam2html.pl -d $TUMOR_BAM -g $GERMLINE_BAM -i $ANALYSIS_PATH/crestSM/$TUMOR_PREF.predSV.txt --ref_genome $REFERENCE_FASTA -o $ANALYSIS_PATH/crestSM/$TUMOR_PREF.predSV.html
chmod 0660 $ANALYSIS_PATH/crestSM/$TUMOR_PREF*

