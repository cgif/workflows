#!/bin/bash
#
# script to run irods_deploy_fastq 
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
PATH_PROJECT_TAG_DIR=#pathProjectTagDir
SEQ_RUN_DATE=#seqRunDate
SEQ_RUN_NAME=#seqRunName
RUN_DIR_BCL2FASTQ=#runDirBcl2Fastq
CUSTOMER_FILE_PATH=#customerFilePath
PROJECT_TAG=#projectTag
MAIL_TEMPLATE_PATH=#mailTemplatePath
PATH_TO_DESTINATION=#pathToDestination
#DEPLOYMENT_SERVER=#deploymentServer
#DEPLOYMENT_TAR_BASE_DIR=#deploymentTarPath
#DEPLOYMENT_SYMBOLIC_LINK=#deploymentSymbolicLink
HIGHTLIGHT="iRODSUserTagging:Star"

IRODS_USER=igf
IRODS_PWD=igf

#ADDING FISTAQ FILES TO ELIOT(eliotResc)
module load irods/4.2.0
iinit igf

echo "`$NOW` tarring the archive of $SEQ_RUN_DATE ..."
ssh login.cx1.hpc.ic.ac.uk "cd $PATH_TO_DESTINATION; tar hcfz $SEQ_RUN_DATE.tar.gz  $SEQ_RUN_DATE"	

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
        echo -e "subject:Sequencing Run $SEQ_RUN_NAME TAR Processing Error - MD5 check failed\nThe MD5 check for the file transfer of sequencing run $SEQ_RUN_NAME failed. Processing aborted." | sendmail -f igf -F "Imperial BRC Genomics Facility" "igf@ic.ac.uk"

        #...and exit
        exit 1
fi

#now send mail to the customer
customers_info=`grep -w $PROJECT_TAG $CUSTOMER_FILE_PATH/customerInfo.csv`
customer_name=`echo $customers_info|cut -d ',' -f2`
customer_username=`echo $customers_info|cut -d ',' -f3`
customer_passwd=`echo $customers_info|cut -d ',' -f4`
customer_email=`echo $customers_info|cut -d ',' -f5`

echo "UTENTE $customer_username"
# check if is internal customer
ldapUser=`ldapsearch -x -h unixldap.cc.ic.ac.uk | grep "uid: $customer_username"`
retval=$?
if [ $retval -ne 0 ]; then
    echo "External customer"
    externalUser="Y"
fi

echo "$NOW checking if user already exists ..."
irods_user=`iadmin lu | grep $customer_username | cut -d "#" -f1`
echo "$NOW irods_user $irods_user"
# if the user has not yet been created, then we create him
if [ "$irods_user" = "" ]
then
	echo "$NOW creating user ..."
	# make user
	iadmin mkuser $customer_username#igfZone rodsuser
	#external user set a password
	if [ $externalUser eq "Y" ]; then
		iadmin moduser $customer_username#igfZone password $customer_passwd
	fi
fi
#ichmod -rM own igf /igfZone/home/$customer_username
#ichmod -rM inherit /igf/Zone/home/$customer_username

# creates the deploy structure
imkdir -p /igfZone/home/$customer_username/$PROJECT_TAG/bam/$SEQ_RUN_DATE

#ichmod -rM own igf /igfZone/home/$customer_username/$PROJECT_TAG/bam/$SEQ_RUN_DATE
#ichmod -rM inherit /igf/Zone/home/$customer_username/$PROJECT_TAG/bam/$SEQ_RUN_DATE
echo "$NOW attaching meta-data run_name to run_date collection ..."
imeta add -C /igfZone/home/$customer_username/$PROJECT_TAG/bam/$SEQ_RUN_DATE run_name $SEQ_RUN_NAME

echo "$NOW storing file in irods .... checksum"
iput -k -fP -R eliotResc $PATH_TO_DESTINATION/$SEQ_RUN_DATE.tar.gz  /igfZone/home/$customer_username/$PROJECT_TAG/bam/$SEQ_RUN_DATE

#set expire date
isysmeta mod /igfZone/home/$customer_username/$PROJECT_TAG/bam/$SEQ_RUN_DATE/$SEQ_RUN_DATE.tar.gz '+30d'
imeta add -d /igfZone/home/$customer_username/$PROJECT_TAG/bam/$SEQ_RUN_DATE/$SEQ_RUN_DATE.tar.gz "$TODAY - bam - $PROJECT_TAG" $customer_username $HIGHTLIGHT
imeta add -d /igfZone/home/$customer_username/$PROJECT_TAG/bam/$SEQ_RUN_DATE/$SEQ_RUN_DATE.tar.gz retention "30" "days"

ichmod -rM read $customer_username /igfZone/home/$customer_username/

#now remove the tar file
echo "`$NOW` remove tar from eliot server ..."
ssh login.cx1.hpc.ic.ac.uk "rm $PATH_TO_DESTINATION/$SEQ_RUN_DATE.tar.gz*" 
echo "`$NOW` Files have been deployed, Well done!"

if [[ $customer_email != *"@"* ]]; then
	#send email alert...
	echo -e "subject:Sequencing Run $SEQ_RUN_NAME Deploying Warning - the email address for $customer_username is unknown." | sendmail -f igf -F "Imperial BRC Genomics Facility" "igf@ic.ac.uk"
fi
customer_mail=customer_mail.$PROJECT_TAG
if [ $externalUser = "Y" ]; then
	cp $MAIL_TEMPLATE_PATH/ecustomer_mail.tml $RUN_DIR_BCL2FASTQ/$customer_mail
else
	cp $MAIL_TEMPLATE_PATH/icustomer_mail.tml $RUN_DIR_BCL2FASTQ/$customer_mail
fi
chmod 770 $RUN_DIR_BCL2FASTQ/$customer_mail
sed -i -e "s/#customerEmail/$customer_email/" $RUN_DIR_BCL2FASTQ/$customer_mail
sed -i -e "s/#customerName/$customer_name/" $RUN_DIR_BCL2FASTQ/$customer_mail
sed -i -e "s/#customerUsername/$customer_username/" $RUN_DIR_BCL2FASTQ/$customer_mail
sed -i -e "s/#passwd/$customer_passwd/" $RUN_DIR_BCL2FASTQ/$customer_mail
sed -i -e "s/#projectName/$PROJECT_TAG/" $RUN_DIR_BCL2FASTQ/$customer_mail
sed -i -e "s/#projectRunDate/$SEQ_RUN_DATE/g" $RUN_DIR_BCL2FASTQ/$customer_mail
sendmail -t < $RUN_DIR_BCL2FASTQ/$customer_mail 
#now remove 
rm $RUN_DIR_BCL2FASTQ/$customer_mail
sed -i /$PROJECT_TAG/d $CUSTOMER_FILE_PATH/customerInfo.csv
