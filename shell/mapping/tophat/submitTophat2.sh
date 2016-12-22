#!/bin/bash

#
# sets up directory structure and configures scripts
# for tophat2 mapping run and submits jobs to pqcgi queue
#

#CONFIGURATION
##############

#now
NOW="date +%Y-%m-%d%t%T%t"

#today
TODAY=`date +%Y-%m-%d`

WALLTIME_HOURS_PER_RUN=72
MEMORY_PER_RUN=20gb
TMPSPACE_PER_RUN=25gb

#trhreads per run
THREADS_PER_RUN=16
#CRAM
THREADS_CRAM2BAM=4


ANALYSIS_NAME=tophat

TOPHAT_VERSION=2.0.10
SAMTOOLS_VERSION=0.1.19
BOWTIE_VERSION=2.1.0
JAVA_VERSION=jdk-7u25
PICARD_VERSION=1.85
PATTERN_READ1='_R1_'
PATTERN_READ2='_R2_'
MULT_READS=F
EDIT_DIST0=F
LIBRARY_TYPE=NONE

FASTQ_GEN_DATE=#fastqGenDate
DATA_VOL_IGF=#dataVolIgf
PATH_SAMPLE_SHEET_REFORMATTED=#pathSampleSheetReformatted
PATH_PROJECT_TAG_DIR=#pathProjectTagDir
SEQ_RUN_DATE=#mseqRunDate
SEQ_RUN_NAME=#mseqRunName
PROJECT_TAG=#mprojectTag
TOPHAT_SCRIPTS_DIR=#TophatScriptDir
DEPLOYMENT_SERVER=#sdeploymentServer
DEPLOYMENT_BASE_DIR=#deploymentBaseDir
QUEUE=#queue
INPUT_SEQRUN_DIR=#inputSeqrunDir

template_path=$TOPHAT_SCRIPTS_DIR/../../../../data-management/templates

#bam_dependencies=afterok
bam_dependencies=""
cram_dependencies=afterok
for d in $PATH_PROJECT_TAG_DIR/fastq/$SEQ_RUN_DATE/*/
do
        sample_name=$(basename $d)
        [[ $sample_name =~ ^(SampleSheet)$ ]] && continue
	sample_sheet_row=`grep  ",$sample_name," $PATH_SAMPLE_SHEET_REFORMATTED`
        ## checks if there was an error 
        retval=$?
        if [ $retval -ne 0 ]; then
                echo "`$NOW` sample $sample_name not found in sample sheet"
                continue
        fi
        lane=`echo $sample_sheet_row | cut -f2 -d ',' | perl -pe 's/\s//g'`
        #..for samples in current lane...
        index=`echo $sample_sheet_row | cut -f5 -d ',' | perl -pe 's/\s//g'`
	path_reads_dir=$DATA_VOL_IGF/rawdata/$PROJECT_TAG/fastq/$SEQ_RUN_DATE/$sample_name
        fastq_read1=${SEQ_RUN_NAME}_${index}_L00${lane}_R1_001.fastq.gz
        fastq_read2=${SEQ_RUN_NAME}_${index}_L00${lane}_R2_001.fastq.gz

        fastq_read1_no_ext=`basename $fastq_read1 .gz`
        fastq_read2_no_ext=`basename $fastq_read2 .gz`

        #get reference FASTA file names
        species=`echo $sample_sheet_row | cut -f4 -d ',' | cut -f1 -d ':' | perl -pe 's/\s//g'`
        assembly=`echo $sample_sheet_row | cut -f4 -d ',' | cut -f2 -d ':' | perl -pe 's/\s//g'`

	PATH_REFERENCE_FASTA=/project/tgu/resources/reference/$species/$assembly/fasta/$assembly.fa
	PATH_BOWTIE2_INDEX=/project/tgu/resources/reference/$species/$assembly/index/bowtie2/$assembly
	PATH_ANNOTATION_GFF=/project/tgu/resources/reference/$species/$assembly/annotation/$assembly.transcripts.gff
	PATH_DICTIONARY=/project/tgu/resources/reference/$species/$assembly/dict/$assembly.dict

        #create directory structure for bam files;
        #each file will contain the same reads as in the original fastq file
        path_results_dir=$DATA_VOL_IGF/rawdata/$PROJECT_TAG/bam/$SEQ_RUN_DATE/$sample_name
        path_results_dir_cram=$DATA_VOL_IGF/rawdata/$PROJECT_TAG/cram/$SEQ_RUN_DATE/$sample_name

        if [ ! -e path_results_dir ]; then
                echo "`$NOW`creating BAM results directory $path_results_dir"
                mkdir -m 770 -p $path_results_dir
                echo "`$NOW`creating CRAM directory $path_results_dir_cram"
                mkdir -m 770 -p $path_results_dir_cram
		chmod -R 770 $DATA_VOL_IGF/rawdata/$PROJECT_TAG/bam
		chmod -R 770 $DATA_VOL_IGF/rawdata/$PROJECT_TAG/cram
        fi

        path_to_tophat_dir=$DATA_VOL_IGF/runs/$PROJECT_TAG/tophat/$FASTQ_GEN_DATE
        path_run_dir=$DATA_VOL_IGF/runs/$PROJECT_TAG/tophat/$FASTQ_GEN_DATE/$sample_name
        path_scripts_dir=$path_run_dir/run
        path_mapping_dir=$path_run_dir/mapping
        path_tmp_dir=$path_run_dir/tmp

        if [ ! -e path_run_dir ]; then
                echo "`$NOW`creating BAM analysis directories in $path_results_dir"
                mkdir -m 770 -p $path_run_dir
                mkdir -m 770 -p $path_scripts_dir
                mkdir -m 770 -p $path_mapping_dir
                mkdir -m 770 -p $path_tmp_dir
		chmod -R 770 $DATA_VOL_IGF/runs/$PROJECT_TAG
        fi

        setup_log=$path_scripts_dir/setup.log
        echo -n "" > $setup_log
        job_id_list=$path_scripts_dir/job_id.list
        echo -n "" > $job_id_list
        chmod 660 $job_id_list

	for fastq_read in `ls --color=never $path_reads_dir/*.f*q*.gz | grep $PATTERN_READ1`
	do
 
	        fastq_read1=`basename $fastq_read`
    		fastq_read2=`echo $fastq_read1 | perl -pe "s/$PATTERN_READ1/$PATTERN_READ2/"`

                #right filter the shortest match (chop end extension)
		path_reference_fasta_no_ext=${PATH_REFERENCE_FASTA%.*}
		reference_fasta_name=`basename $path_reference_fasta_no_ext`

                #right filter the longest match (chop end extension)
		read_group_name=${fastq_read1%%.*}
     
	        #output prefix
		output_prefix=$read_group_name.vs.$reference_fasta_name

		echo "`$NOW`setting up $ANALYSIS_NAME mapping run for $read_group_name... "

		echo "`$NOW`setting up $ANALYSIS_NAME run" >> $setup_log
		echo "`$NOW`read directory: $path_reads_dir" >> $setup_log
		echo "`$NOW`fastq file 1: $fastq_read1" >> $setup_log
		echo "`$NOW`fastq file 2: $fastq_read2" >> $setup_log
		echo "`$NOW`reference file: $PATH_REFERENCE_FASTA" >> $setup_log
		echo "`$NOW`annotation file: $PATH_ANNOTATION_GFF" >> $setup_log
		echo "`$NOW`bowtie2 index directory: $PATH_BOWTIE2_INDEX" >> $setup_log
		echo "`$NOW`dictionary file: $PATH_DICTIONARY" >> $setup_log
		echo "`$NOW`script directory: $path_scripts_dir" >> $setup_log
		echo "`$NOW`result directory: $path_results_dir" >> $setup_log
		echo "`$NOW`creating and submitting job scripts:" >> $setup_log

		 script_path=$path_scripts_dir/tophat2.$output_prefix.sh
		cp $TOPHAT_SCRIPTS_DIR/tophat2.sh $script_path
		chmod 770 $script_path

	        #set variables 
		sed -i -e "s/#tophatVersion/$TOPHAT_VERSION/" $script_path
		sed -i -e "s/#bowtieVersion/$BOWTIE_VERSION/" $script_path
		sed -i -e "s/#samtoolsVersion/$SAMTOOLS_VERSION/" $script_path
		sed -i -e "s/#javaVersion/$JAVA_VERSION/" $script_path
		sed -i -e "s/#picardVersion/$PICARD_VERSION/" $script_path
		sed -i -e "s/#walltimeHours/$WALLTIME_HOURS_PER_RUN/" $script_path
		sed -i -e "s/#threads/$THREADS_PER_RUN/" $script_path
                sed -i -e "s/#memory/$MEMORY_PER_RUN/" $script_path
                sed -i -e "s/#tmpspc/$TMPSPACE_PER_RUN/" $script_path
		sed -i -e "s/#outputPrefix/$output_prefix/" $script_path
		sed -i -e "s/#multReads/$mult_reads/" $script_path
		sed -i -e "s/#editDist0/$edit_dist0/" $script_path
		sed -i -e "s/#libraryType/$library_type/" $script_path
		sed -i -e "s/#pathOutputDir/${path_results_dir//\//\\/}/" $script_path
		sed -i -e "s/#pathReferenceFasta/${PATH_REFERENCE_FASTA//\//\\/}/" $script_path
		sed -i -e "s/#pathAnnotation/${PATH_ANNOTATION_GFF//\//\\/}/" $script_path
		sed -i -e "s/#pathBowtie2Index/${PATH_BOWTIE2_INDEX//\//\\/}/" $script_path
		sed -i -e "s/#pathDictionary/${PATH_DICTIONARY//\//\\/}/" $script_path
		sed -i -e "s/#pathReadsDirectory/${path_reads_dir//\//\\/}/" $script_path
		sed -i -e "s/#read1/$fastq_read1/" $script_path
		sed -i -e "s/#read2/$fastq_read2/" $script_path

	        #submit job and save job ID to dependency variable 
		log_output_path=`echo $script_path | perl -pe 's/\.sh/\.log/g'`
		echo "`$NOW`$script_path" >> $setup_log
		echo -n "`$NOW`" >> $setup_log
		job_id=`qsub -q $QUEUE -o $log_output_path $script_path`
		bam_dependencies=$bam_dependencies:$job_id
		echo $job_id >> $setup_log
		echo -e "$job_id" >> $job_id_list
		echo "#############################################################################" >> $setup_log
	done

	#submit cram conversion
        bam2cram_script_path=$path_scripts_dir/bam2cram.$output_prefix.sh
        cp $TOPHAT_SCRIPTS_DIR/../samtools/bam2cram.sh $bam2cram_script_path
        chmod 770 $bam2cram_script_path

        path_input_bam=$path_results_dir/$output_prefix.sorted.bam

        sed -i -e "s/threads/$THREADS_CRAM2BAM/" $bam2cram_script_path
        sed -i -e "s/pathInputBam/${path_input_bam//\//\\/}/" $bam2cram_script_path
        sed -i -e "s/outputPrefix/$output_prefix/" $bam2cram_script_path
        sed -i -e "s/pathOutputCramDir/${path_results_dir_cram//\//\\/}/" $bam2cram_script_path
        sed -i -e "s/pathReferenceFastaNoExt/${PATH_REFERENCE_FASTA//\//\\/}/" $bam2cram_script_path

        LOG_PATH=`echo $bam2cram_script_path | perl -pe 's/\.sh/\.log/g'`

        echo "`$NOW`submitting bam2cram job:" 
        echo "`$NOW`bam2cram.$output_prefix.sh"
        echo -n "`$NOW`"

        cram_job_id=`qsub -q $QUEUE -W depend=afterok${bam_dependencies} -o $LOG_PATH $bam2cram_script_path`
        echo "qsub -q $QUEUE -W depend=$bam_dependencies -o $LOG_PATH $bam2cram_script_path"
        echo $cram_job_id

done

#summary script
tophat_summary_results=$DATA_VOL_IGF/runs/$PROJECT_TAG/tophat/$FASTQ_GEN_DATE
tophat_summary_deployment=$DEPLOYMENT_BASE_DIR/project/$PROJECT_TAG/tophat/$FASTQ_GEN_DATE
SUMMARY_SCRIPT=$path_to_tophat_dir/summary_tophat2.$output_prefix.pl
cp $TOPHAT_SCRIPTS_DIR/summary_tophat2.pl $SUMMARY_SCRIPT
chmod 770 $SUMMARY_SCRIPT

sed -i -e "s/#projectDirAnalysis/${path_run_dir//\//\\/}/" $SUMMARY_SCRIPT
sed -i -e "s/#projectDirResults/${path_results_dir//\//\\/}/" $SUMMARY_SCRIPT
sed -i -e "s/#deploymentServer/$DEPLOYMENT_SERVER/" $SUMMARY_SCRIPT
sed -i -e "s/#summaryDeployment/${tophat_summary_deployment//\//\\/}/" $SUMMARY_SCRIPT
sed -i -e "s/#summaryResults/${tophat_summary_results//\//\\/}/" $SUMMARY_SCRIPT

SUM_DEPENDENCIES=afterany${bam_dependencies}
SUMMARY_LOG=`echo $SUMMARY_SCRIPT | perl -pe 's/\.pl/\.log/g'`
echo "`$NOW`submitting summary script:"
echo "`$NOW`$SUMMARY_SCRIPT"
echo -n "`$NOW`"

SUM_JOB_ID=`qsub -q $QUEUE -o $SUMMARY_LOG  -j oe -W depend=$SUM_DEPENDENCIES -M igf@imperial.ac.uk $SUMMARY_SCRIPT`
echo "qsub -q $QUEUE -o $SUMMARY_LOG  -j oe -W depend=$SUM_DEPENDENCIES -M igf@imperial.ac.uk $SUMMARY_SCRIPTr"
echo $SUM_JOB_ID

#############################################################################
###       configure script that tar bam files and  deploys in irods      ####
#############################################################################
deploy_irods_bam_script=$path_to_tophat_dir/irods_deploy_bam.${PROJECT_TAG}.sh

cp $TOPHAT_SCRIPTS_DIR/../bwa/irods_deploy_bam.sh $deploy_irods_bam_script
chmod 770 $deploy_irods_bam_script
path_project_tag_dir=$DATA_VOL_IGF/rawdata/$PROJECT_TAG/bam
sed -i -e "s/#seqRunDate/$SEQ_RUN_DATE/" $deploy_irods_bam_script
sed -i -e "s/#seqRunName/$SEQ_RUN_NAME/" $deploy_irods_bam_script
sed -i -e "s/#runDirBcl2Fastq/${path_to_tophat_dir//\//\\/}/" $deploy_irods_bam_script
sed -i -e "s/#customerFilePath/${INPUT_SEQRUN_DIR//\//\\/}/" $deploy_irods_bam_script
sed -i -e "s/#projectTag/$PROJECT_TAG/" $deploy_irods_bam_script
sed -i -e "s/#mailTemplatePath/${template_path//\//\\/}/" $deploy_irods_bam_script
sed -i -e "s/#pathToDestination/${path_project_tag_dir//\//\\/}/" $deploy_irods_bam_script

#submit job 
log_output_path=`echo $deploy_irods_bam_script | perl -pe 's/\.sh/\.log/g'`
echo -n "" > $log_output_path
echo -n "`$NOW`submitting tarresult job: " 
echo "$deploy_irods_bam_script"

## deploy BAM files
#job_id=`qsub -q $QUEUE -W depend=$bam_dependencies -o $log_output_path -j oe $deploy_irods_bam_script`
#echo "qsub -q $QUEUE -W depend=$bam_dependencies -o $log_output_path -j oe $deploy_irods_bam_script"
#echo "`$NOW`Job ID:$job_id"
#chmod 660 $log_output_path

