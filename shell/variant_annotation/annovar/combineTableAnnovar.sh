
#!/bin/bash

#script to combine summary variant annotation 
# generated with ANNOVAR table_annovar.pl

#PBS -l walltime=5:00:00
#PBS -l select=1:ncpus=1:mem=20gb

#PBS -M cgi@imperial.ac.uk
#PBS -m ea
#PBS -j oe

ANNOVAR_ANNOTATIONS="#annovarAnnotations"
OUTPUT_PATH=#outputPath
OUTPUT_PATH_EXONIC=#outputPathExonic
GENE_TYPE=#geneType
BASEDIR=#baseDir
ANNOVAR_DB=#annovarDb

for FILE in `cat $ANNOVAR_ANNOTATIONS`
do
	ANNOVAR_ANNOTATION_FILES="$ANNOVAR_ANNOTATION_FILES $FILE"
done

$BASEDIR/combineTableAnnovar.pl $ANNOVAR_ANNOTATION_FILES > $OUTPUT_PATH.temp

#add associated gene names if Ensembl based annotation
#if [[ GENE_TYPE == "ensgene" ]]
#then
#	$BASEDIR/addAssociatedGeneName.pl $OUTPUT_PATH.gz $ANNOVAR_DB/hg19_ensembl_associated_gene_names.tsv | gzip > $OUTPUT_PATH.tmp.gz
#	mv $OUTPUT_PATH.tmp.gz $OUTPUT_PATH.gz
#fi

# add sample names columns
$BASEDIR/addSampleColumns.pl $OUTPUT_PATH.temp $OUTPUT_PATH

head -n 1 $OUTPUT_PATH > $OUTPUT_PATH_EXONIC
grep -E "exonic|splicing" $OUTPUT_PATH >> $OUTPUT_PATH_EXONIC

gzip -f $OUTPUT_PATH
gzip -f $OUTPUT_PATH_EXONIC

chmod 660 $OUTPUT_PATH.gz
chmod 660 $OUTPUT_PATH_EXONIC.gz

rm $OUTPUT_PATH.temp



