#!/bin/bash

## script to add read group tag to bam file and clip reads if required

#PBS -l walltime=24:00:00
#PBS -l select=1:ncpus=1:mem=5gb:tmpspace=#tmpSpacemb

#PBS -M cgi@imperial.ac.uk
#PBS -m ae
#PBS -j oe

#PBS -q pqcgi

# load modules
module load gatk/#gatkVersion
module load java/#javaVersion
module load samtools/#samtoolsVersion
module load picard/#picardVersion

NOW="date +%Y-%m-%d%t%T%t"
JAVA_XMX=4800M

# define variables
INPUT_BAM=#inputBam
INPUT_BAM_NAME=`basename $INPUT_BAM .bam`
READ_GROUP=#readGroup
CLIP_FILE=#clipFile
RG_FILE=#RGFile
REFERENCE_FASTA=#referenceFasta
REFERENCE_FASTA_DICT=`echo $REFERENCE_FASTA | perl -pe 's/\.fa/\.dict/g'`
RESULTS_DIR=#resultsDir

cp $INPUT_BAM $TMPDIR/original.bam
cp $INPUT_BAM.bai $TMPDIR/original.bam.bai
cp $REFERENCE_FASTA $TMPDIR/reference.fa
cp $REFERENCE_FASTA.fai $TMPDIR/reference.fa.fai
cp $REFERENCE_FASTA_DICT $TMPDIR/reference.dict
cp $CLIP_FILE $TMPDIR/clip.txt
cp $RG_FILE $TMPDIR/RG.txt

READS_BEFORE_CLIP=`samtools view -c $TMPDIR/original.bam`

CLIP_CYCLE=`grep -P "^$INPUT_BAM\t" $TMPDIR/clip.txt| perl -e '$_=<>; @cols=split(/\t/, $_); print $cols[2];'`

CLIP_READ=`grep -P "^$INPUT_BAM\t" $TMPDIR/clip.txt| perl -e '$_=<>; @cols=split(/\t/, $_); print $cols[3];'`
echo "Clipping cycles $CLIP_CYCLE reads $CLIP_READ from file $INPUT_BAM"
CLIP_READ1=`echo $CLIP_READ|grep 1`
CLIP_READ2=`echo $CLIP_READ|grep 2`

RG_INFO=`grep -P "ID:$READ_GROUP\t" $TMPDIR/RG.txt |perl -e '$_=<>; @cols=split(/\t/,$_); foreach $col(@cols) {next if $col =~ /^PI/; next unless $col=~/(\w+):(.*)/; $tag=$1; $val=$2; $val =~ s/\s/_/g; print "$tag=$val\t"}'`

#add read group info to the bam file
echo "`${NOW}`Adding read group info to the bam file"
java -XX:+UseSerialGC -Xmx$JAVA_XMX -jar $PICARD_HOME/AddOrReplaceReadGroups.jar INPUT=$TMPDIR/original.bam OUTPUT=$TMPDIR/rg.bam SORT_ORDER=coordinate $RG_INFO VALIDATION_STRINGENCY=SILENT TMP_DIR=$TMPDIR VERBOSITY=WARNING

if [ "$CLIP_CYCLE" != "" ]
then

    #make separate bam files for read1, read2 and unmapped reads (ClipReads will filter unmapped reads)
    samtools view -b -f 64 -F 4 $TMPDIR/rg.bam > $TMPDIR/first.bam
    samtools view -b -f 128 -F 4 $TMPDIR/rg.bam > $TMPDIR/second.bam
    samtools view -b -f 4 $TMPDIR/rg.bam > $TMPDIR/unmapped.bam
    samtools view -H $TMPDIR/rg.bam > $TMPDIR/header.sam

    #clip read1 if required
    if [ "$CLIP_READ1" != "" ] 
    then

	echo "`${NOW}`Clipping read1"
	samtools index $TMPDIR/first.bam
	java -Xmx$JAVA_XMX -XX:+UseSerialGC -Djava.io.tmpdir=$TMPDIR -jar $GATK_HOME/GenomeAnalysisTK.jar \
	     -T ClipReads \
	     -R $TMPDIR/reference.fa \
	     -I $TMPDIR/first.bam \
	     -o $TMPDIR/first.clipped.bam \
	     -os $TMPDIR/first.clipped.stats \
	     -CR SOFTCLIP_BASES \
	     -CT "$CLIP_CYCLE" \
	     -rf UnmappedRead

    else

	cp $TMPDIR/first.bam $TMPDIR/first.clipped.bam

    fi

    #clip read2 if required
    if [ "$CLIP_READ2" != "" ] 
    then

	echo "`${NOW}`Clipping read2"
	samtools index $TMPDIR/second.bam
	java -Xmx$JAVA_XMX -XX:+UseSerialGC -Djava.io.tmpdir=$TMPDIR -jar $GATK_HOME/GenomeAnalysisTK.jar \
	     -T ClipReads \
	     -R $TMPDIR/reference.fa \
	     -I $TMPDIR/second.bam \
	     -o $TMPDIR/second.clipped.bam \
	     -os $TMPDIR/second.clipped.stats \
	     -CR SOFTCLIP_BASES \
	     -CT "$CLIP_CYCLE" \
	     -rf UnmappedRead

    else

	cp $TMPDIR/second.bam $TMPDIR/second.clipped.bam

    fi

    #merge back read1, read2, unmapped and add header
    echo "`${NOW}`Merging back read1, read2, unmapped and add header"
    samtools merge -h $TMPDIR/header.sam $TMPDIR/complete.bam $TMPDIR/first.clipped.bam $TMPDIR/second.clipped.bam $TMPDIR/unmapped.bam 

else

    cp $TMPDIR/rg.bam $TMPDIR/complete.bam

fi

READS_AFTER_CLIP=`samtools view -c $TMPDIR/complete.bam`

#compare number of reads befor and after clipping
if [ $READS_BEFORE_CLIP != $READS_AFTER_CLIP ]
then
    echo "`${NOW}`Number of reads before and after clipping is not the same!!!!!!!!!!!!!!!"
    exit 1
fi

#copy clipped files
echo "`${NOW}`Copy bam files"

if [[ ! -d $RESULTS_DIR ]]; then
	mkdir -p $RESULTS_DIR
fi

cp $TMPDIR/complete.bam $RESULTS_DIR/$INPUT_BAM_NAME.bam
echo "`${NOW}`Done"

