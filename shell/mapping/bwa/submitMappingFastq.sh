#!/bin/bash
#
# script to run tarBcl2FastqResults 
#

#PBS -l walltime=72:00:00
#PBS -l select=1:ncpus=1:mem=1024mb

#PBS -m ea
#PBS -M cgi@imperial.ac.uk
#PBS -j oe

#PBS -q pqcgi


#now
NOW="date +%Y-%m-%d%t%T%t"

#today
TODAY=`date +%Y-%m-%d`

#set up script
DATA_VOL_IGF=#mdataVolIgf
PATH_SAMPLE_SHEET_REFORMATTED=#mpathSampleSheetReformatted
PATH_PROJECT_TAG_DIR=#mpathProjectTagDir
SEQ_RUN_DATE=#mseqRunDate
SEQ_RUN_NAME=#mseqRunName
PROJECT_TAG=#mprojectTag
BWA_SCRIPTS_DIR=#mbwaScriptDir
DEPLOYMENT_SERVER=#mdeploymentServer
DEPLOYMENT_BASE_DIR=#mdeploymentBaseDir
QUEUE=#mqueue
INPUT_SEQRUN_DIR=#minputSeqrunDir


#SPLIT FASTQ
FASTQ_FILE_SIZE_KB=20000000
THREADS_PER_RUN_SPLITFASTQ=4
WALLTIME_HOURS_SPLITFASTQ=72

#BWA
READS_PER_RUN_BWA=10000000
THREADS_PER_RUN_BWA=2
PATTERN_READ_1=_R1_
PATTERN_READ_2=_R2_
WALLTIME_HOURS_PER_RUN_BWA=30

SAMTOOLS_VERSION=1.1
BWA_VERSION=0.7.5a

template_path=$BWA_SCRIPTS_DIR/../../../../data-management/templates

bam_dependencies=afterok
for sample_name in `ls --color=never $PATH_PROJECT_TAG_DIR/fastq/$SEQ_RUN_DATE/`
do
	sample_sheet_row=`grep $sample_name $PATH_SAMPLE_SHEET_REFORMATTED`
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
        echo " XXXXX species $species"
         assembly=`echo $sample_sheet_row | cut -f4 -d ',' | cut -f2 -d ':' | perl -pe 's/\s//g'`
        echo " XXXXX assembly $assembly"
         sequence_type=`echo $sample_sheet_row | cut -f4 -d ',' | cut -f3 -d ':' | perl -pe 's/\s//g'`
        echo " XXXXX sequence_type $sequence_type"

         reference_fasta_index=$DATA_VOL_IGF/resources/reference/eukaryote/$species/$assembly/index/bwa/$assembly.fa
        if [[ $sequence_type = "rna" ]]
        then
        	reference_fasta_index=$DATA_VOL_IGF/resources/reference/eukaryote/$species/$assembly/index/bowtie2/$assembly.fa
        fi
        echo " XXXXX reference_fasta_index $reference_fasta_index"

         reference_fasta_name=`basename $reference_fasta_index .gz`

        #get reference FASTA directory path
         path_reference_fasta_dir=`dirname $reference_fasta_index`

         path_reference_fasta_no_ext=$path_reference_fasta_dir/$reference_fasta_name
         path_reference_dict_dir=`echo $path_reference_fasta_dir | perl -pe 's/\/index\/bwa\//\/dict\//'`
         reference_dict_name=`basename $reference_fasta_name .fa`
         reference_dict_name=${reference_dict_name}.dict
         path_reference_dict=$path_reference_dict_dir/$reference_dict_name

        echo "`$NOW`creating and submitting job scripts:"

        #create directory structure for merged bam files;
        #each file will contain the same reads as in the original fastq file
         path_results_dir=$DATA_VOL_IGF/rawdata/$PROJECT_TAG/bam/$SEQ_RUN_DATE/$sample_name
        # XXXX non lo sto usando
         path_results_dir_cram=$DATA_VOL_IGF/rawdata/$PROJECT_TAG/cram/$SEQ_RUN_DATE/$sample_name

        if [ ! -e path_results_dir ]; then
        	echo "`$NOW`creating BWA results directory $path_results_dir"
        	mkdir -m 770 -p $path_results_dir
                # XXXX  non los to usando
        	echo "`$NOW`creating CRAM directory $path_results_dir_cram"
                mkdir -m 770 -p $path_results_dir_cram
        fi
	

	 path_to_bwa_dir=$DATA_VOL_IGF/analysis/$PROJECT_TAG/bwa/$TODAY
         path_run_dir=$DATA_VOL_IGF/analysis/$PROJECT_TAG/bwa/$TODAY/$sample_name
         path_scripts_dir=$path_run_dir/run
         path_mapping_dir=$path_run_dir/mapping
         path_tmp_dir=$path_run_dir/tmp

        if [ ! -e path_run_dir ]; then
        	echo "`$NOW`creating BWA analysis directories in $path_results_dir"
        	mkdir -m 770 -p $path_run_dir
                mkdir -m 770 -p $path_scripts_dir
                mkdir -m 770 -p $path_mapping_dir
                mkdir -m 770 -p $path_tmp_dir
        fi

        SETUP_LOG=$path_scripts_dir/setup.log
        echo -n "" > $SETUP_LOG


	merge_files=""
	merge_dependencies=afterok
	for fastq_read1_split in `ls --color=never $path_tmp_dir/${fastq_read1_no_ext}_split/*.f*q* | grep -v $PATTERN_READ_2`
	do

       		#remove path information
       	 	fastq_read1_split=`basename $fastq_read1_split`

        	#get name of read2 fastq by preplacing read pair tag
        	fastq_read2_split=`echo $fastq_read1_split | perl -pe "s/$PATTERN_READ_1/$PATTERN_READ_2/"`

        	path_reads_fastq_read1_split=$path_tmp_dir/${fastq_read1_no_ext}_split/$fastq_read1_split
        	path_reads_fastq_read2_split=$path_tmp_dir/${fastq_read1_no_ext}_split/$fastq_read2_split


        	#output prefix
        	output_prefix=$fastq_read1_split.vs.$reference_fasta_name

        	mapping_script_path=$path_scripts_dir/bwaAlignPe.$output_prefix.sh
        	cp $BWA_SCRIPTS_DIR/bwaAlignPe.sh $mapping_script_path
        	chmod 770 $mapping_script_path

        	#set variables 
        	sed -i -e "s/#walltimeHours/$WALLTIME_HOURS_PER_RUN_BWA/" $mapping_script_path
        	sed -i -e "s/#threads/$THREADS_PER_RUN_BWA/" $mapping_script_path
        	sed -i -e "s/#outputPrefix/$output_prefix/" $mapping_script_path
        	sed -i -e "s/#multReads/$MULT_READS/" $mapping_script_path
        	sed -i -e "s/#pathOutputDir/${path_mapping_dir//\//\\/}/" $mapping_script_path
        	sed -i -e "s/#pathReferenceFastaNoExt/${path_reference_fasta_no_ext//\//\\/}/" $mapping_script_path
        	sed -i -e "s/#pathReadsFastqRead1NoExt/${path_reads_fastq_read1_split//\//\\/}/" $mapping_script_path
        	sed -i -e "s/#pathReadsFastqRead2NoExt/${path_reads_fastq_read2_split//\//\\/}/" $mapping_script_path
        	sed -i -e "s/#pathReferenceDict/${path_reference_dict//\//\\/}/" $mapping_script_path
        	sed -i -e "s/#pathReferenceIdxDir/${path_reference_dict_dir//\//\\/}/" $mapping_script_path

        	#submit job and save job ID to dependency variable 
        	log_output_path=`echo $mapping_script_path | perl -pe 's/\.sh/\.log/g'`
        	echo "`$NOW`bwaAlignPe.$output_prefix.sh"
        	echo -n "`$NOW`"
        	#align_job_id=`qsub -W depend=$splitfastq_dependencies -o $log_output_path $mapping_script_path`
        	align_job_id=`qsub -o $log_output_path $mapping_script_path`
        	echo $align_job_id
        	merge_dependencies=$merge_dependencies:$align_job_id
        	merge_files="$merge_files $output_prefix.unsorted.bam"
	done;

	#submit merging jobs
	output_prefix=$fastq_read1_no_ext.vs.$reference_fasta_name

	merge_script_path=$path_scripts_dir/samtoolsMerge.$output_prefix.sh
	cp $BWA_SCRIPTS_DIR/../samtools/samtoolsMergeAndDelete.sh $merge_script_path
	chmod 770 $merge_script_path

	sed -i -e "s/outputPrefix/$output_prefix/" $merge_script_path
	sed -i -e "s/inputDir/${path_mapping_dir//\//\\/}/" $merge_script_path
	sed -i -e "s/pathOutputDir/${path_results_dir//\//\\/}/" $merge_script_path
	sed -i -e "s/inBam/\"${merge_files//\//\\/}\"/" $merge_script_path

	log_output_path=`echo $merge_script_path | perl -pe 's/\.sh/\.log/g'`
	echo "`$NOW`submitting merging script:"
	echo "`$NOW`samtoolsMerge.$output_prefix.sh"
	echo -n "`$NOW`"

	merge_job_id=`qsub  -o $log_output_path -W depend=$merge_dependencies $merge_script_path`
	echo "merge_job_id: $merge_job_id"
	bam_dependencies=$bam_dependencies:$merge_job_id

	#summary script
	bwa_summary_results=$DATA_VOL_IGF/analysis/$PROJECT_TAG/bwa/$TODAY
        bwa_summary_deployment=$DEPLOYMENT_BASE_DIR/project/$PROJECT_TAG/bwa/$TODAY
	SUMMARY_SCRIPT=$path_scripts_dir/summary_bwa.$output_prefix.pl
	cp $BWA_SCRIPTS_DIR/summary_bwa.pl $SUMMARY_SCRIPT
	chmod 770 $SUMMARY_SCRIPT

	sed -i -e "s/projectDirAnalysis/${path_run_dir//\//\\/}/" $SUMMARY_SCRIPT
	sed -i -e "s/projectDirResults/${path_results_dir//\//\\/}/" $SUMMARY_SCRIPT
	sed -i -e "s/deploymentServer/$DEPLOYMENT_SERVER/" $SUMMARY_SCRIPT
	sed -i -e "s/summaryDeployment/${bwa_summary_deployment//\//\\/}/" $SUMMARY_SCRIPT
	sed -i -e "s/summaryResults/${bwa_summary_results//\//\\/}/" $SUMMARY_SCRIPT

	SUM_DEPENDENCIES=afterany:$merge_job_id									
	SUMMARY_LOG=`echo $SUMMARY_SCRIPT | perl -pe 's/\.pl/\.log/g'`
	echo "`$NOW`submitting summary script:"
	echo "`$NOW`$SUMMARY_SCRIPT"
	echo -n "`$NOW`"

	SUM_JOB_ID=`qsub -q $QUEUE -o $SUMMARY_LOG  -j oe -W depend=$SUM_DEPENDENCIES -M cgi@imperial.ac.uk $SUMMARY_SCRIPT`	
	echo $SUM_JOB_ID

done

#########################################
### create e configure tarJobScript ####
deploy_irods_bam_script=$path_to_bwa_dir/irods_deploy_bam.${PROJECT_TAG}.sh

cp $BWA_SCRIPTS_DIR/irods_deploy_bam.sh $deploy_irods_bam_script
chmod 770 $deploy_irods_bam_script
path_project_tag_dir=$DATA_VOL_IGF/rawdata/$PROJECT_TAG/bam
sed -i -e "s/#seqRunDate/$SEQ_RUN_DATE/" $deploy_irods_bam_script
sed -i -e "s/#seqRunName/$SEQ_RUN_NAME/" $deploy_irods_bam_script
sed -i -e "s/#runDirBcl2Fastq/${path_to_bwa_dir//\//\\/}/" $deploy_irods_bam_script
sed -i -e "s/#customerFilePath/${INPUT_SEQRUN_DIR//\//\\/}/" $deploy_irods_bam_script
sed -i -e "s/#projectTag/$PROJECT_TAG/" $deploy_irods_bam_script
sed -i -e "s/#mailTemplatePath/${template_path//\//\\/}/" $deploy_irods_bam_script
sed -i -e "s/#pathToDestination/${path_project_tag_dir//\//\\/}/" $deploy_irods_bam_script

#submit job 
log_output_path=`echo $deploy_irods_bam_script | perl -pe 's/\.sh/\.log/g'`
echo -n "" > $log_output_path
echo -n "`$NOW`submitting tarresult job: " 
echo "$deploy_irods_bam_script"

job_id=`qsub -q $QUEUE -W depend=$bam_dependencies -o $log_output_path -j oe $deploy_irods_bam_script`
echo "qsub -q $QUEUE -W depend=$bam_dependencies -o $log_output_path -j oe $deploy_irods_bam_script"
echo "`$NOW`Job ID:$job_id"
chmod 660 $log_output_path
