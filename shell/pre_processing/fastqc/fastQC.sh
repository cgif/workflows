#!/bin/bash

#
# script to run FastQC on a fastq file
# on cx1
#

#PBS -l walltime=24:00:00
#PBS -l select=1:ncpus=1:mem=500mb

#PBS -m ea
#PBS -M cgi@imperial.ac.uk
#PBS -j oe

#PBS -q pqcgi

module load fastqc/0.11.2

#now
NOW="date +%Y-%m-%d%t%T%t"

#today
TODAY=`date +%Y-%m-%d`

#CONFIGURATION
##############

#path to reads fastq file
PATH_READS_DIRECTORY=#pathReadsFastq
#PATTERN_READ1=#patternRead1
#PATTERN_READ2=#patternRead2
FASTQ_READ1=#fastqRead1
FASTQ_READ2=#fastqRead2
SINGLE_READS=#singleReads


#QC output directory
PATH_QC_REPORT_DIR=#pathQcReportDir

#deployment
DEPLOYMENT_SERVER=#deploymentServer
DEPLOYMENT_PATH=#deploymentPath
SUMMARY_PATH=#summaryPath

#create temporary QC report output directory
mkdir $TMPDIR/qc

#copy fastqs to tmp space
echo "`${NOW}`copying fasq files to temporary scratch space..."
echo "`${NOW}`$FASTQ_READ1"
cp $PATH_READS_DIRECTORY/$FASTQ_READ1 $TMPDIR/$FASTQ_READ1

# checks if FASTQ_READ2 exists. If doesn't exists assume RUN sigle read
if [ ! -f "$PATH_READ_DIRECTORY/$FASTQ_READ2" ]; then
	SINGLE_READS = "T"
fi

if [[ "$SINGLE_READS" == "F" ]]; then
	echo "`${NOW}`$FASTQ_READ2"
	cp $PATH_READS_DIRECTORY/$FASTQ_READ2 $TMPDIR/$FASTQ_READ2
fi
           
#check if mate file found and the number of lines in mate files is the same
gzip -t $TMPDIR/$FASTQ_READ1
if [ $? -ne "0" ]; then
	echo "`${NOW}`ERROR:File $FASTQ_READ1 is corrupted. Skipped." 
elif [[ "$SINGLE_READS" == "F" ]]; then

	gzip -t $TMPDIR/$FASTQ_READ2
	if [ $? -ne "0" ]; then
		echo "`${NOW}`ERROR:File $FASTQ_READ2 is corrupted. Skipped." 
    else 

		#compare number of reads
		COUNT_LINES_READ1=`gzip -d -c $TMPDIR/$FASTQ_READ1 | wc -l | cut -f 1 -d ' '`
		COUNT_LINES_READ2=`gzip -d -c $TMPDIR/$FASTQ_READ2 | wc -l | cut -f 1 -d ' '`

		if [ $COUNT_LINES_READ1 -ne $COUNT_LINES_READ2 ]; then
			echo "ERROR:Unequal number of lines in the mate files. Skipped." 
		else

			#run FastQC 
			#--noextract   creating a zip archived report
			#--nogroup     disable grouping of bases for reads >50bp 

			echo "`${NOW}`running FastQC..."
			$FASTQC_HOME/bin/fastqc -o $TMPDIR/qc --noextract --nogroup  $TMPDIR/$FASTQ_READ1
			$FASTQC_HOME/bin/fastqc -o $TMPDIR/qc --noextract --nogroup  $TMPDIR/$FASTQ_READ2

		fi	
	fi
else
#run FastQC for single reads
	echo "`${NOW}`running FastQC..."
	$FASTQC_HOME/bin/fastqc -o $TMPDIR/qc --noextract --nogroup  $TMPDIR/$FASTQ_READ1

fi


#copy results to output folder
echo "`${NOW}`copying zipped QC report to $PATH_QC_REPORT_DIR..."
cp $TMPDIR/qc/*zip $PATH_QC_REPORT_DIR
chmod 660 $PATH_QC_REPORT_DIR/*zip

echo "`${NOW}`creating deployment directory $DEPLOYMENT_PATH on $DEPLOYMENT_SERVER..."
ssh $DEPLOYMENT_SERVER "mkdir -p -m 775 $DEPLOYMENT_PATH" < /dev/null

for ZIP in `ls $TMPDIR/qc/*.zip`
do

        echo "`${NOW}`uncompressing QC report zip archive $ZIP..."
	unzip $ZIP
	REPORT_DIR=`basename $ZIP .zip`	
		
	echo "`${NOW}`deploying $REPORT_DIR report directory to $DEPLOYMENT_SERVER:$DEPLOYMENT_PATH..."
	scp -r $TMPDIR/$REPORT_DIR $DEPLOYMENT_SERVER:$DEPLOYMENT_PATH/  < /dev/null 
	ssh $DEPLOYMENT_SERVER "chmod 775 $DEPLOYMENT_PATH/$REPORT_DIR" < /dev/null 
	ssh $DEPLOYMENT_SERVER "chmod 775 $DEPLOYMENT_PATH/$REPORT_DIR/*"  < /dev/null
	ssh $DEPLOYMENT_SERVER "chmod 664 $DEPLOYMENT_PATH/$REPORT_DIR/*/*"  < /dev/null

    echo "`${NOW}`copying unzipped QC report to $PATH_QC_REPORT_DIR/$REPORT_DIR"
	mkdir -p -m 770 $PATH_QC_REPORT_DIR/$REPORT_DIR
	cp $TMPDIR/$REPORT_DIR/*.txt  $PATH_QC_REPORT_DIR/$REPORT_DIR
	chmod 660 $PATH_QC_REPORT_DIR/$REPORT_DIR/*.txt

done



  

