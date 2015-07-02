#!/bin/bash

## script returns list of genes that overlap with CNVs in at least 2 patients

#PBS -l walltime=72:00:00
#PBS -l select=1:ncpus=1:mem=10gb

#PBS -m ea
#PBS -M cgi@imperial.ac.uk
#PBS -j oe

#PBS -q pqcgi

#module load bedtools/#bedtoolsVersion

#RESULTS_DIR=#resultsDir
#EXON_BED=#exonBed
#DEPLOYMENT_SERVER=#deploymentServer
#SUMMARY_DEPLOYMENT=#summaryDeployment
#PROJECT=#project

module load bedtools/2.13.3

RESULTS_DIR=/groupvol/cgi/results/johnson_glioma/FREEC/2015-04-22
EXON_BED=/groupvol/cgi/resources/annotations/eukaryote/human/GRCh37/GRCh37.74.exon.bed
DEPLOYMENT_SERVER=eliot.med.ic.ac.uk
SUMMARY_DEPLOYMENT=/www/html/report/project/johnson_glioma/FREEC/2015-02-22
PROJECT=johnson_glioma

CNV_BED=$TMPDIR/FREEC_for_ann.bed
CNV_TRACK=$TMPDIR/FREEC.bed
CNV_GEN=$TMPDIR/FREEC.gen
CNV_ANN=$TMPDIR/FREEC.ann

cp $EXON_BED $TMPDIR/exon.bed

#create index HTML file 
INDEX=$TMPDIR/index.html
printf "<HTML><BODY>\n" > $INDEX
printf "<H2>Graphical presentation of the discovered CNVs</H2>\n" >> $INDEX
printf "<TABLE CELLPADDING = 5>\n" >> $INDEX

#create bed files for annotation and for opening in Genome Browser
echo -n "" > $CNV_BED
echo "browser position chr1\n" > $CNV_TRACK
echo "browser hide all\n" >> $CNV_TRACK

SAMPLE_COUNT=0
for SAMPLE in `ls $RESULTS_DIR|grep -v multisample`; do 

	echo "Analysing at CNVs for sample $SAMPLE"
	FREEC=$RESULTS_DIR/$SAMPLE/tumor.bam_CNVs.p.value.txt
	#FREEC=$RESULTS_DIR/$SAMPLE/${SAMPLE}.CNVs.txt
	FREEC_TMP=$TMPDIR/$SAMPLE.CNVs.txt
	cp $FREEC $FREEC_TMP

	#make bed file
	export SAMPLE_NAME=`echo ${SAMPLE##*vs.}`
	cat $FREEC_TMP | grep -vE '^chr|X|Y' | perl -MPOSIX -e 'while(<>) {chomp(); @data=split(/\t/,$_); print "$data[0]\t$data[1]\t$data[2]\t$ENV{SAMPLE_NAME}:$data[4]\n";}' >> $CNV_BED
	echo "track name=$SAMPLE_NAME type=bed visibility=3 db=hg19 itemRgb=On\n" >> $CNV_TRACK
	cat $FREEC_TMP | grep -vE '^chr' | perl -MPOSIX -e 'while(<>) {chomp(); @data=split(/\t/,$_); print "$data[0]\t$data[1]\t$data[2]\t$ENV{SAMPLE_NAME}:$data[4]\n";}' >> $CNV_TRACK

	#copy FREEC images
	PNG=$RESULTS_DIR/$SAMPLE/tumor.bam_ratio.txt.png
	#PNG=$RESULTS_DIR/$SAMPLE/${SAMPLE}.png
	scp -r $PNG $DEPLOYMENT_SERVER:$SUMMARY_DEPLOYMENT/$PROJECT.$SAMPLE_NAME.png > /dev/null 2>&1
	if [[ $((SAMPLE_COUNT % 3)) -eq 0 ]]; then
		printf "<TR>" >> $INDEX
	fi
	printf "<TD><FONT SIZE = +1><B>$SAMPLE_NAME</B><P><FONT COLOR = 'blue'><IMG BORDER=1 WIGTH=300 HEIGHT=300 HSPACE=5 VSPACE=5 SRC=$PROJECT.$SAMPLE_NAME.png>\n" >> $INDEX
	SAMPLE_COUNT=$(( $SAMPLE_COUNT + 1 ))

done

printf "</TABLE>\n" >> $INDEX

#run bestools to select genes overlapping CNVs
intersectBed -a $CNV_BED -b $TMPDIR/exon.bed -wa -wb > $CNV_ANN

#print list of genes that were affected in > 2 patients
cat $CNV_ANN | cut -f 4,8 | sort | uniq | perl -e 'while(<>) {chomp(); @data=split(/\t/,$_); $data[0]=~s/(.*):(.*)/$2/; print "$data[0]\t$data[1]\n"}' | sort | uniq -c > $CNV_GEN


printf "<P><H2>Number of genes that overlap with CNVs in multiple samples (autosomes only)</H2>\n" >> $INDEX
printf "<TABLE CELLPADDING = 5>\n" >> $INDEX
printf "<TR><TH><FONT SIZE = +1>#Samples<TH><FONT SIZE = +1>#Genes\n" >> $INDEX

#print counts table
for I in `seq 1 9`; do

	COUNT=`cat $CNV_GEN|grep -E "^\s+$I\s"|wc -l`
	cat $CNV_GEN|grep -E "^\s+$I\s" > $TMPDIR/$PROJECT.$I.genes.txt
	scp -r $TMPDIR/$PROJECT.$I.genes.txt $DEPLOYMENT_SERVER:$SUMMARY_DEPLOYMENT/$PROJECT.$I.genes.txt > /dev/null 2>&1
	printf "<TR><TD>$I<TD><A HREF = $PROJECT.$I.genes.txt>$COUNT</A>\n" >> $INDEX

done

printf "</TABLE>\n" >> $INDEX

printf "<P><FONT SIZE = +1>List of CNVs in <A HREF = $PROJECT.CNVs.bed>BED</A> format\n" >> $INDEX
printf "</BODY></HTML>\n" >> $INDEX

cp $INDEX $RESULTS_DIR/multisample/index.html
cp $CNV_TRACK $RESULTS_DIR/multisample/${PROJECT}.bed
cp $CNV_ANN $RESULTS_DIR/multisample/${PROJECT}.annotation
chmod -R 0664 $RESULTS_DIR/multisample/*

scp -r $INDEX $DEPLOYMENT_SERVER:$SUMMARY_DEPLOYMENT/index.html
scp -r $CNV_TRACK $DEPLOYMENT_SERVER:$SUMMARY_DEPLOYMENT/$PROJECT.CNVs.bed
ssh $DEPLOYMENT_SERVER chmod -R 664 $SUMMARY_DEPLOYMENT/*

