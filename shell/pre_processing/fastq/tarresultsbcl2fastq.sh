#!/bin/bash
#
# script to run tarBcl2FastqResults 
#

#PBS -l walltime=#walltimeHours:00:00
#PBS -l select=1:ncpus=#threads:mem=1024mb:tmpspace=#tmpSpacegb

#PBS -m ea
#PBS -M cgi@imperial.ac.uk
#PBS -j oe

#PBS -q #queu


#now
NOW="date +%Y-%m-%d%t%T%t"

#today
TODAY=`date +%Y-%m-%d`

#set up script
PATH_PROJECT_TAG_DIR=#pathProjectTagDir
SEQ_RUN_DATE=#seqRunDate
PATH_TO_DESTINATION=#pathToDestination
DEPLOYMENT_SERVER=#deploymentServer
DEPLOYMENT_TAR_BASE_DIR=#deploymentTarPath
DEPLOYMENT_SYMBOLIC_LINK=#deploymentSymbolicLink

echo "`$NOW` tarring the archive of $SEQ_RUN_DATE ..."
ssh login.cx1.hpc.ic.ac.uk "tar cfz $PATH_TO_DESTINATION/$SEQ_RUN_DATE.tar.gz  $PATH_PROJECT_TAG_DIR/$SEQ_RUN_DATE"
echo "`$NOW` tar of $SEQ_RUN_DATE completed"

#generate an md5 checksum for the tarball
#need to change to location of archive to generate md5
echo "`$NOW` Generating md5 checksum for TAR archive..."
ssh login.cx1.hpc.ic.ac.uk "cd $PATH_TO_DESTINATION; md5sum $SEQ_RUN_DATE.tar.gz > $SEQ_RUN_DATE.tar.gz.md5; chmod 664 $SEQ_RUN_DATE.tar.gz $SEQ_RUN_DATE.tar.gz.md5"
echo "`$NOW` md5 checksum Generated"

#change to location where the tar and the md5 file are & check 
MD5_STATUS=`ssh login.cx1.hpc.ic.ac.uk "cd $PATH_TO_DESTINATION; md5sum -c $SEQ_RUN_DATE.tar.gz.md5 2>&1 | head -n 1 | cut -f 2 -d ' '"`
echo  $MD5_STATUS

#abort if md5 check fails
if [[ $MD5_STATUS == 'FAILED' ]]
then
        #send email alert...
        echo -e "subject:Sequencing Run $SEQ_RUN_DATE TAR Processing Error - MD5 check failed\nThe MD5 check for the file transfer of sequencing run $SEQ_RUN_DATE failed. Processing aborted." | sendmail -f igf -F "Imperial BRC Genomics Facility" "mmuelle1@ic.ac.uk"

        #...and exit
        exit 1
fi

# creates rnd name for result directory
rnddir_results=`cat /dev/urandom | tr -dc 'a-zA-Z0-9' | fold -w ${1:-15} | head -n 1` 
PATH_TO_RNDDIR=$DEPLOYMENT_TAR_BASE_DIR/$rnddir_results
ssh $DEPLOYMENT_SERVER "mkdir -m 770 -p $PATH_TO_RNDDIR" 

echo "`$NOW` coping TAR archive on eliot server ..."
scp -r $PATH_TO_DESTINATION/$SEQ_RUN_DATE.tar.gz $DEPLOYMENT_SERVER:$PATH_TO_RNDDIR 
#create project_tag dir & symbolic link
ssh $DEPLOYMENT_SERVER "mkdir -m 770 -p $DEPLOYMENT_SYMBOLIC_LINK"
ssh $DEPLOYMENT_SERVER "ln -s  $PATH_TO_RNDDIR $DEPLOYMENT_SYMBOLIC_LINK/fastq"

#change to location where the tar and the md5 file are & check
MD5_STATUS=`ssh $DEPLOYMENT_SERVER "cd $PATH_TO_RNDDIR; md5sum -c $SEQ_RUN_DATE.tar.gz.md5 2>&1 | head -n 1 | cut -f 2 -d ' '"`
echo  $MD5_STATUS

#abort if md5 check fails
if [[ $MD5_STATUS == 'FAILED' ]]
then
        #send email alert...
        echo -e "subject:Sequencing Run $SEQ_RUN_DATE Deploying Error - MD5 check failed\nThe MD5 check for the file transfer of sequencing run $SEQ_RUN_DATE failed. Processing aborted." | sendmail -f igf -F "Imperial BRC Genomics Facility" "mmuelle1@ic.ac.uk"

        #...and exit
        exit 1
fi

#now remove the tar file
echo "`$NOW` remove tar from eliot server ..."
ssh login.cx1.hpc.ic.ac.uk "rm $PATH_TO_DESTINATION/$SEQ_RUN_DATE.tar.gz*" 
echo "`$NOW` Files have been deployed, Well done!"

