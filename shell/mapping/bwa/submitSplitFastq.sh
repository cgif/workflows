#!/bin/bash
#
# script to run splitFastq for each sample by project_tag 
#

#PBS -l walltime=72:00:00
#PBS -l select=1:ncpus=1:mem=1024mb

#PBS -m ea
#PBS -M igf@imperial.ac.uk
#PBS -j oe

#PBS -q pqcgi


#now
NOW="date +%Y-%m-%d%t%T%t"

#today
TODAY=`date +%Y-%m-%d`

#DATA_VOL_IGF=/project/tgu
#PATH_SAMPLE_SHEET_REFORMATTED=/project/tgu/runs/seqrun/150724_SN172_0481_BC7D1MACXX-TEST4/bcl2fastq/2015-11-30/BC7D1MACXX-TEST4.csv
#PATH_PROJECT_TAG_DIR=/project/tgu/rawdata/bam_test
#SEQ_RUN_DATE=2015-07-24
#SEQ_RUN_NAME=150724_SN172_0481_BC7D1MACXX-TEST4
#PROJECT_TAG=bam_test
#BWA_SCRIPTS_DIR=/home/mcosso/git.mdev/workflows/shell/mapping/bwa


#set up script
DATA_VOL_IGF=#dataVolIgf
PATH_SAMPLE_SHEET_REFORMATTED=#pathSampleSheetReformatted
PATH_PROJECT_TAG_DIR=#pathProjectTagDir
SEQ_RUN_DATE=#seqRunDate
SEQ_RUN_NAME=#seqRunName
PROJECT_TAG=#projectTag
BWA_SCRIPTS_DIR=#bwaScriptDir
DEPLOYMENT_SERVER=#deploymentServer
DEPLOYMENT_BASE_DIR=#deploymentBaseDir
QUEUE=#queue
INPUT_SEQRUN_DIR=#inputSeqrunDir
REMOVE_BAMS=#removeBams



#SPLIT FASTQ
FASTQ_FILE_SIZE_KB=20000000
THREADS_PER_RUN_SPLITFASTQ=4
WALLTIME_HOURS_SPLITFASTQ=72

#BWA
READS_PER_RUN_BWA=10000000
THREADS_PER_RUN_BWA=2
WALLTIME_HOURS_PER_RUN_BWA=30

splitfastq_dependencies="afterok"
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


        echo "`$NOW`creating and submitting job scripts:"

	 path_to_bwa_dir=$DATA_VOL_IGF/runs/$PROJECT_TAG/bwa/$TODAY 
         path_run_dir=$DATA_VOL_IGF/runs/$PROJECT_TAG/bwa/$TODAY/$sample_name
         path_scripts_dir=$path_run_dir/run
         path_tmp_dir=$path_run_dir/tmp

        if [ ! -e path_run_dir ]; then
        	echo "`$NOW`creating BWA analysis directories in $path_results_dir"
        	mkdir -m 770 -p $path_run_dir
                mkdir -m 770 -p $path_scripts_dir
                mkdir -m 770 -p $path_tmp_dir
        fi

        SETUP_LOG=$path_scripts_dir/setup.log
        echo -n "" > $SETUP_LOG

        #temporary directory for splitfastq results
        splitfastq_output_dir=$path_tmp_dir/${fastq_read1_no_ext}_split
        mkdir -m 770 -p $splitfastq_output_dir

        ########################################
        #submit splitting jobs

        echo "`$NOW`submitting jobs to split fastq files into $READS_PER_RUN chunks... " >> $SETUP_LOG
        for fastq in $path_reads_dir/$fastq_read1 $path_reads_dir/$fastq_read2
        do
	        fastq_name=`basename $fastq .gz`

                #calculate required temp space
                 file_size_mb=$(( $FASTQ_FILE_SIZE_KB / 1024 ))
                 tmp_space_mb=$(( $file_size_mb * 2 ))

                 script_path=$path_scripts_dir/splitFastq.$fastq_name.sh
                cp $BWA_SCRIPTS_DIR/splitFastq.sh $script_path
                chmod 770 $script_path

                sed -i -e "s/#walltimeHours/$WALLTIME_HOURS_SPLITFASTQ/" $script_path
                sed -i -e "s/#threads/$THREADS_PER_RUN_SPLITFASTQ/" $script_path
                sed -i -e "s/#tmpSpace/$tmp_space_mb/" $script_path
                sed -i -e "s/#inputFastq/${fastq//\//\\/}/" $script_path
                sed -i -e "s/#outputDir/${splitfastq_output_dir//\//\\/}/" $script_path
                sed -i -e "s/#readsPerChunk/$READS_PER_RUN_BWA/" $script_path

                 log_path=`echo $script_path | perl -pe 's/\.sh/\.log/g'`

                echo "`$NOW`submitting fastq splitting job:" >> $SETUP_LOG
                echo "`$NOW`splitFastq.$fastq_name.sh" >> $SETUP_LOG
                echo  -n "`$NOW`" >> $SETUP_LOG

                echo "qsub -o $log_path $script_path"
                split_job_id=`qsub -o $log_path $script_path`
                echo $split_job_id >> $SETUP_LOG

                splitfastq_dependencies="$splitfastq_dependencies:$split_job_id"
       done;
       echo "splitfastq dependencies $splitfastq_dependencies"

done
mapping_script_path=$path_to_bwa_dir/submitMappingFastq.$PROJECT_TAG.sh
cp $BWA_SCRIPTS_DIR/submitMappingFastq.sh $mapping_script_path
chmod 770 $mapping_script_path

#set variables 
sed -i -e "s/#mfastqGenDate/$TODAY/" $mapping_script_path
sed -i -e "s/#mdataVolIgf/${DATA_VOL_IGF//\//\\/}/" $mapping_script_path
sed -i -e "s/#mpathSampleSheetReformatted/${PATH_SAMPLE_SHEET_REFORMATTED//\//\\/}/" $mapping_script_path
sed -i -e "s/#mpathProjectTagDir/${PATH_PROJECT_TAG_DIR//\//\\/}/" $mapping_script_path
sed -i -e "s/#mseqRunDate/$SEQ_RUN_DATE/" $mapping_script_path
sed -i -e "s/#mseqRunName/$SEQ_RUN_NAME/" $mapping_script_path
sed -i -e "s/#mprojectTag/$PROJECT_TAG/" $mapping_script_path
sed -i -e "s/#mbwaScriptDir/${BWA_SCRIPTS_DIR//\//\\/}/" $mapping_script_path
sed -i -e "s/#mdeploymentServer/$DEPLOYMENT_SERVER/" $mapping_script_path
sed -i -e "s/#mdeploymentBaseDir/${DEPLOYMENT_BASE_DIR//\//\\/}/" $mapping_script_path
sed -i -e "s/#mqueue/$QUEUE/" $mapping_script_path
sed -i -e "s/#minputSeqrunDir/${INPUT_SEQRUN_DIR//\//\\/}/" $mapping_script_path
sed -i -e "s/#mremoveBams/$REMOVE_BAMS/" $mapping_script_path


#submit job and save job ID to dependency variable 
log_output_path=`echo $mapping_script_path | perl -pe 's/\.sh/\.log/g'`
echo "`$NOW` mapping_project_fastq.$PROJECT_TAG.sh"
echo -n "`$NOW`"
align_job_id=`qsub -W depend=$splitfastq_dependencies -o $log_output_path $mapping_script_path`
echo $align_job_id
