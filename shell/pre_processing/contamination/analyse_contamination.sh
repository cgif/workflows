#!/bin/bash

#
# script to create a summary file
#

#PBS -l walltime=1:00:00
#PBS -l select=1:ncpus=1:mem=10gb
#PBS -M cgi@imperial.ac.uk
#PBS -m ea
#PBS -j oe

INPUT_DIR="#input_dir"
OUTPUT_FILE="#output_file"
PATTERN_READ_1="#pattern_read1"
PATTERN_READ_2="#pattern_read2"

printf "sample\treads_R1\tnon-human_pct\tnon-nr_pct\treads_R2\tnon-human_pct\tnon-nr_pct\tcontaminant_species\tpct\n" > $OUTPUT_FILE

for SAMPLE in `ls $INPUT_DIR|grep -v multisample`; do

	SPECIES=""
	TOTAL_READS=0
	TOTAL_READS_1=0
	TOTAL_READS_2=0
	TOTAL_NON_HUMAN_COUNT_1=0
	TOTAL_NON_HUMAN_COUNT_2=0
	TOTAL_UNMAPPED_COUNT_1=0
	TOTAL_UNMAPPED_COUNT_2=0

	for INPUT_FASTQ_1 in `ls $INPUT_DIR/$SAMPLE/imsa/*fastq|grep $PATTERN_READ_1`; do

		LINE_COUNT_FASTQ_1=`cat $INPUT_FASTQ_1|wc -l`
		INPUT_COUNT_FASTQ_1=$((LINE_COUNT_FASTQ_1/4))

		READ_GROUP_1=`basename $INPUT_FASTQ_1 .fastq`

		NON_HUMAN_1="$INPUT_DIR/$SAMPLE/imsa/$READ_GROUP_1"_bhg_lhg_lhg_shg.fa
		NON_HUMAN_COUNT_1=`grep -P "^>" $NON_HUMAN_1|wc -l`

		UNMAPPED_1="$INPUT_DIR/$SAMPLE/imsa/$READ_GROUP_1"_bhg_lhg_lhg_shg_snt.fa
		UNMAPPED_COUNT_1=`grep -P "^>" $UNMAPPED_1|wc -l`

		INPUT_FASTQ_2=`echo $INPUT_FASTQ_1| perl -pe "s/$PATTERN_READ_1/$PATTERN_READ_2/"`
		LINE_COUNT_FASTQ_2=`cat $INPUT_FASTQ_2|wc -l`
		INPUT_COUNT_FASTQ_2=$((LINE_COUNT_FASTQ_2/4))

		READ_GROUP_2=`basename $INPUT_FASTQ_2 .fastq`

		NON_HUMAN_2="$INPUT_DIR/$SAMPLE/imsa/$READ_GROUP_2"_bhg_lhg_lhg_shg.fa
		NON_HUMAN_COUNT_2=`grep -P "^>" $NON_HUMAN_2|wc -l`

		UNMAPPED_2="$INPUT_DIR/$SAMPLE/imsa/$READ_GROUP_2"_bhg_lhg_lhg_shg_snt.fa
		UNMAPPED_COUNT_2=`grep -P "^>" $UNMAPPED_2|wc -l`

		SPECIES_COUNT="$INPUT_DIR/$SAMPLE/imsa/"combined."$READ_GROUP_1"_bhg_lhg_lhg_shg_snt.tax_species.txt
		
		NEW_SPECIES=`cut -f2 $SPECIES_COUNT|grep -v 'Scientific Name'`
		SPECIES=`printf "$SPECIES\n$NEW_SPECIES"`

		TOTAL_READS=`echo "scale=1;$TOTAL_READS + $INPUT_COUNT_FASTQ_1 + $INPUT_COUNT_FASTQ_2"|bc`   
		TOTAL_READS_1=`echo "scale=1;$TOTAL_READS_1 + $INPUT_COUNT_FASTQ_1"|bc` 
		TOTAL_READS_2=`echo "scale=1;$TOTAL_READS_2 + $INPUT_COUNT_FASTQ_2"|bc` 

		TOTAL_NON_HUMAN_COUNT_1=`echo "scale=1;$TOTAL_NON_HUMAN_COUNT_1 + $NON_HUMAN_COUNT_1"|bc` 
		TOTAL_NON_HUMAN_COUNT_2=`echo "scale=1;$TOTAL_NON_HUMAN_COUNT_2 + $NON_HUMAN_COUNT_2"|bc` 

		TOTAL_UNMAPPED_COUNT_1=`echo "scale=1;$TOTAL_UNMAPPED_COUNT_1 + $UNMAPPED_COUNT_1"|bc` 
		TOTAL_UNMAPPED_COUNT_2=`echo "scale=1;$TOTAL_UNMAPPED_COUNT_2 + $UNMAPPED_COUNT_2"|bc` 

	done

	TOTAL_READS_1_PCT=`echo "scale=1;$TOTAL_READS_1/100"|bc`
	TOTAL_NON_HUMAN_COUNT_1_PCT=`echo "scale=1;$TOTAL_NON_HUMAN_COUNT_1/$TOTAL_READS_1_PCT"|bc`
	TOTAL_UNMAPPED_COUNT_1_PCT=`echo "scale=1;$TOTAL_UNMAPPED_COUNT_1/$TOTAL_READS_1_PCT"|bc`

	TOTAL_READS_2_PCT=`echo "scale=1;$TOTAL_READS_2/100"|bc`
	TOTAL_NON_HUMAN_COUNT_2_PCT=`echo "scale=1;$TOTAL_NON_HUMAN_COUNT_2/$TOTAL_READS_2_PCT"|bc`
	TOTAL_UNMAPPED_COUNT_2_PCT=`echo "scale=1;$TOTAL_UNMAPPED_COUNT_2/$TOTAL_READS_2_PCT"|bc`

	printf "$SAMPLE\t$TOTAL_READS_1\t$TOTAL_NON_HUMAN_COUNT_1_PCT\t$TOTAL_UNMAPPED_COUNT_1_PCT\t$TOTAL_READS_2\t$TOTAL_NON_HUMAN_COUNT_2_PCT\t$TOTAL_UNMAPPED_COUNT_2_PCT" >> $OUTPUT_FILE

	TOTAL_READS_PCT=`echo "scale=1;$TOTAL_READS/100"|bc`
	
	for ID in `echo "$SPECIES"|sort|uniq|grep -vP '^$'|sed 's/ /_/g'|grep -v 'Homo_sapiens'`; do

		ID=`echo $ID|sed 's/_/ /g'`
		TOTAL_COUNT=0
		
		for COUNT in `grep "$ID" $INPUT_DIR/$SAMPLE/imsa/*tax_species.txt|cut -f 3`; do
	
			TOTAL_COUNT=`echo "scale=1;$TOTAL_COUNT+$COUNT"|bc`

		done

		ID_PCT=`echo "scale=1;$TOTAL_COUNT/$TOTAL_READS_PCT"|bc`
		if [ $ID_PCT \> 1 ]; then
			printf "\t$ID\t$ID_PCT" >> $OUTPUT_FILE
		fi

	done

	printf "\n" >> $OUTPUT_FILE

done

