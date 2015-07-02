#!/usr/bin/perl -w

use strict;
use warnings;

my $project_dir_analysis = "projectDirAnalysis";
my $project_dir_results = "projectDirResults";
my $summary_results = "summaryResults";
my $deployment_server = "deploymentServer";
my $summary_deployment = "summaryDeployment";

$project_dir_analysis =~ s/^(\S+)\/\S+$/$1/;
$project_dir_results =~ s/^(\S+)\/\S+$/$1/;

#collect all sample names into an array
my @samples = ();
my ($sample_dir, $log, $stat, $name, $split, $full_name);

opendir (PROJECT, "$project_dir_analysis"); 
while (defined(my $sample = readdir(PROJECT))){
    $sample_dir = "$project_dir_analysis/$sample/run";
    next unless (-d "$sample_dir");
    push (@samples,$sample);
}

#collect data from bwa log files
my %sum = (); 
my %fastq = (); 
my %mate = ();

#for each sample in alphabetical order
foreach my $sample (sort {$a cmp $b} @samples){

    $sample_dir = "$project_dir_analysis/$sample/run";

    #check if fastq files are there
    if (`grep "does not contain any fastq files" $sample_dir/setup.log`){
        $fastq{$sample} = "FAIL";
    }else{
        $fastq{$sample} = "PASS";
    }

    #check if mate fastq files are there
    if (`grep "No mate file found" $sample_dir/setup.log`){
        $mate{$sample} = "FAIL";
    }else{
        $mate{$sample} = "PASS";
    }

    #check if fastq split and bam merge are completed successfully; collect statistics from flagstat file
    opendir (SAMPLE, "$sample_dir"); 
    while (defined($log = readdir(SAMPLE))){
        next unless $log =~ /\.log/;
    
        if ($log =~ /^bwaAlignPe\.(\S+)\.f.*q\.(\w+)\.vs\..*\.log/){
	    $name = $1;
	    $split = $2;

            if (`grep "deleting temporary fastq files" $sample_dir/$log`){
                $sum{$sample}{$name}{$split}{'bam'} = "PASS";
            }else{
                $sum{$sample}{$name}{$split}{'bam'} = "FAIL";
            }

	}elsif ($log =~ /^samtoolsMerge\.((\S+)\.f.*q\.vs\..*)\.log/){
	    $full_name = $1;
	    $name = $2;

            if (`grep "deleting intermediate BAM files" $sample_dir/$log`){
                $sum{$sample}{$name}{'00'}{'merge'} = "PASS";
            }else{
                $sum{$sample}{$name}{'00'}{'merge'} = "FAIL";
            }

	    $stat = "$project_dir_results/$sample/$full_name".".sorted.bam.flagstat";
	    open (STAT,"$stat")||print "Can't open $stat\n";
	    while(<STAT>){
	        $sum{$sample}{$name}{'00'}{'total'} = "$1<BR>" if /^(\d+)\s+\+\s+\d+\s+in total/;
                $sum{$sample}{$name}{'00'}{'mapped'} = "$1($2)<BR>" if  /^(\d+)\s+\+\s+\d+\s+mapped \((\S+):.*/;
	    }
	}
    }
}

#print extracted data into summary file
open (OUT, ">$summary_results/index.html");
print OUT "<HTML>";
print OUT "<TABLE CELLPADDING=5><TR>";
print OUT "<TH><CENTER>Sample";
print OUT "<TH><CENTER>Fastq files";
print OUT "<TH><CENTER>Mate files";
print OUT "<TH><CENTER>Fastq input";
print OUT "<TH><CENTER>Split";
print OUT "<TH><CENTER>BAM alignment";
print OUT "<TH><CENTER>Merging BAM";
print OUT "<TH><CENTER>Total reads";
print OUT "<TH><CENTER>Mapped reads\n";

#f1 is flagging for the first read group for each sample
my $f1 = 0;

#for each sample in alphabetical order
foreach my $sample (sort {$a cmp $b} @samples){

    #print sample name, results of fastq check and mate check 
    print OUT "<TR><TD>$sample";
    $f1 = 1;

    if ($fastq{$sample} eq "PASS"){
	print OUT "<TD><CENTER><IMG SRC=tick.png ALT=PASS>";
    }elsif ($fastq{$sample} eq "FAIL"){
	print OUT "<TD><CENTER><IMG SRC=error.png ALT=FAIL>";
    }

    if ($mate{$sample} eq "PASS"){
	print OUT "<TD><CENTER><IMG SRC=tick.png ALT=PASS>";
    }elsif ($mate{$sample} eq "FAIL"){
	print OUT "<TD><CENTER><IMG SRC=error.png ALT=FAIL>";
    }

    #foreach read group
    foreach my $group (sort {$a cmp $b} keys %{$sum{$sample}}){

    	#print read group name, results of fastq split check and bam merge check, mapping stats
	if ($f1){
	    print OUT "<TD>$group";
	    $f1 = 0;
	}else{
	    print OUT "<TR><TD> <TD> <TD> <TD>$group";
	}

        foreach my $split (sort {$a cmp $b} keys %{$sum{$sample}{$group}}){

	    if ($split eq '00'){
                print OUT "<TD>$split\n";
	    }else{
                print OUT "<TR><TD> <TD> <TD> <TD> <TD>$split\n";
	    }

            if ($sum{$sample}{$group}{$split}{'bam'} eq "PASS"){
	        print OUT "<TD><CENTER><IMG SRC=tick.png ALT=PASS>";
            }elsif ($sum{$sample}{$group}{$split}{'bam'} eq "FAIL"){
	        print OUT "<TD><CENTER><IMG SRC=error.png ALT=FAIL>";
	    }
	 
            if ($split eq '00'){

                if ($sum{$sample}{$group}{$split}{'merge'} eq "PASS"){
	            print OUT "<TD><CENTER><IMG SRC=tick.png ALT=PASS>";
                }elsif ($sum{$sample}{$group}{$split}{'merge'} eq "FAIL"){
	            print OUT "<TD><CENTER><IMG SRC=error.png ALT=FAIL>";
	        }

                print OUT "<TD><CENTER>$sum{$sample}{$group}{$split}{'total'}" if defined($sum{$sample}{$group}{$split}{'total'});
                print OUT "<TD><CENTER>$sum{$sample}{$group}{$split}{'mapped'}" if defined($sum{$sample}{$group}{$split}{'mapped'});

            }

	    print OUT "\n";
	}

    }

}
print OUT "</TABLE></HTML>";

#deploy HTML summary to the Web server
system("scp -r $summary_results/index.html $deployment_server:$summary_deployment/index.html");
system("ssh $deployment_server chmod -R 664 $summary_deployment/index.html");
