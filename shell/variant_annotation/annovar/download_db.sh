#!/bin/bash

# script to download annovar databases
# specified in a tab separated text file
# in the format <db_name>\t<webfrom>\t<description>

BUILD=hg19
TODAY=`date +%Y-%m-%d`
ANNOVAR_PATH=/ax3-cgi/resources/annovar/$TODAY

DATABASE_LIST=$1


NOW="date +%Y-%m-%d%t%T%t"

LOG_OUT=$ANNOVAR_PATH/download_db.$TODAY.log

#redirect stdout and stderr to terminal and log file
exec > >(tee $LOG_OUT)
exec 2>&1

ANNOVAR_DB_PATH=$ANNOVAR_PATH/db/$BUILD
echo "`${NOW}`creating output directory $ANNOVAR_DB_PATH..."
mkdir -p $ANNOVAR_DB_PATH


echo "`${NOW}`downloading databases specified in $DATABASE_LIST..."

cut -f1,2 $DATABASE_LIST | while read DB WEBFROM
do

	echo "`${NOW}`$DB..."
	$ANNOVAR_PATH/annotate_variation.pl -downdb -buildver $BUILD -webfrom $WEBFROM $DB $ANNOVAR_DB_PATH
	
done

echo "`${NOW}`creating sequence files..."
#create sequence files for Encode/Genecode gene annotations
GENCODE_COUNT=`grep 'EncodeGencode' $DATABASE_LIST | wc -l`

if [[ $GENCODE_COUNT -ge 1 ]]
then

	mkdir -p $ANNOVAR_DB_PATH/tmp_seq
	cd $ANNOVAR_DB_PATH/tmp_seq
	wget --timestamping 'ftp://hgdownload.cse.ucsc.edu/goldenPath/hg19/chromosomes/*'
	gzip -d *

	cd $ANNOVAR_PATH

	cut -f1 $DATABASE_LIST | grep 'EncodeGencode' | while read DB
	do

		$ANNOVAR_PATH/retrieve_seq_from_fasta.pl -format genericGene -seqdir $ANNOVAR_DB_PATH/tmp_seq -outfile $ANNOVAR_DB_PATH/${BUILD}_${DB}Mrna.fa $ANNOVAR_DB_PATH/${BUILD}_${DB}.txt

	done

	rm -r $ANNOVAR_DB_PATH/tmp_seq

fi

chmod -R 770 $ANNOVAR_PATH/db
