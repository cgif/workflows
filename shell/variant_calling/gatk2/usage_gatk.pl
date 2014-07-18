#!/usr/bin/perl
 
use DBI;
use strict;

my $setup_log = "/groupvol/cgi/analysis/gatk2/project/cgi_vcopt/setup.2013-07-05.log";
my $pwd = "/groupvol/cgi/analysis/gatk2/project/cgi_vcopt/multisample/2013-07-05/run";

my $workflow="gatk2";
my $system="cx1";

#set up database connection
my $db_host="eliot.med.ic.ac.uk";
my $db="wordpress";
my $db_usr="wordpress";
my $db_pwd="wordpress";

my $dbh = DBI->connect("dbi:mysql:$db;host=$db_host",$db_usr,$db_pwd) or die "Connection Error: $DBI::errstr\n";

#get run information
my @tokens = split(/\//, $pwd);
my $run_date = $tokens[@tokens-2];
$run_date=$run_date." 00:00:00";
print "$run_date\n";

my @tokens=split(/\//, $pwd);
my $project_tag = $tokens[@tokens-4];


#update cgi_project table and get project_id
my $project_id = fetch_project_id($project_tag);
print "project_id = $project_id\n";

my $run_id = fetch_run_id($project_id, $run_date);
print "run_id = $run_id\n";






#my $sql = "INSERT INTO cgi_run (project_id, date_time) VALUES ($project_id,$run_date);";
#my $sth = $dbh->prepare($sql);
#$sth->execute or die "SQL Error: $DBI::errstr\n";
#
#while (my @row = $sth->fetchrow_array()) {
#	print "@row\n";
#}


#functions
##########

sub fetch_project_id {

	my $project_tag = $_[0];

	my $sql_get_project_id = "SELECT id FROM cgi_project WHERE project_tag = '$project_tag';";
	my $sth = $dbh->prepare($sql_get_project_id);
	$sth->execute or die "SQL Error: $DBI::errstr\n";

	my $project_id = 0;

	if($sth->rows == 0){

		my $sql_insert_project = "INSERT INTO cgi_project (project_tag) VALUES ('$project_tag');";
		$sth = $dbh->prepare($sql_insert_project);
		$sth->execute or die "SQL Error: $DBI::errstr\n";

		$sth = $dbh->prepare($sql_get_project_id);
		$sth->execute or die "SQL Error: $DBI::errstr\n";

		my @row = $sth->fetchrow_array();
		$project_id = $row[0];

	} else {

		my @row = $sth->fetchrow_array();
		$project_id = $row[0];

	}
	
	return $project_id;
	
}

sub fetch_run_id {

	my $project_id = $_[0];
	my $run_date = $_[1];

	#update cgi_run table and get run_id...
	my $sql_insert_run = "INSERT INTO cgi_run (project_id, date_time) VALUES ($project_id, '$run_date');";
	my $sth = $dbh->prepare($sql_insert_run);
	$sth->execute or die "SQL Error: $DBI::errstr\n";

	my $sql_get_run_id = "SELECT id FROM cgi_run WHERE project_id=project_id AND date_time='$run_date' ORDER BY id DESC;";
	$sth = $dbh->prepare($sql_get_run_id);
	$sth->execute or die "SQL Error: $DBI::errstr\n";
	my @row = $sth->fetchrow_array();
	my $run_id = $row[0];

	return $run_id;
}
