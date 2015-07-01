#!/usr/bin/perl -w

use strict;
use warnings;

#This script generates on-line summary that reflects jobs progress. 
#Once job is finished and non-empty log file is present in the run directory, 
#its content is tested for the presence of error messages from gatk. 
#If (i) there are no error messages in the log file and 
#(ii) expected output file is not empty then green tick appears in the summary table.
#If at least one of these conditions is not satisfied then red cross appears in the summary table.

my $project_dir_analysis = "#projectDirAnalysis";
my $project_dir_results = "#projectDirResults";
my $date = "#today";
my $summary_results = "#summaryResults";
my $deployment_server = "#deploymentServer";
my $summary_deployment = "#summaryDeployment";
my $sample_list = "#sampleList";
my $ref_chunks = "#refChunks";
my $ref_intervals = "#refIntervals";
my $type = "#type";
my $summary_link = $1 if  $summary_deployment =~ /\/www\/html\/(.*)/;

my %sum = ();

#count number of chunks
open (CHUNKS, "$ref_chunks");
my $total_chunks = 0;
my $chunk = 0;
while (<CHUNKS>){
    chomp();
    $chunk = $1 if /^\S+\t\S+\t\S+\t\S+\t(\S+)/;
    $total_chunks = $chunk if $chunk > $total_chunks;
}

#count number of intervals
if ($type eq "TARGETED" && $ref_intervals ne "null"){
    open (INTERVALS, "$ref_intervals");
    my $total_intervals = 0;
    my $interval = 0;
    while (<INTERVALS>){
        chomp();
        $interval = $1 if /^\S+\t\S+\t\S+\t\S+\t(\S+)/;
        $total_intervals = $interval if $interval > $total_intervals;
    }
}

#collect list of sample names
my @sample = (); my @list = ();
open (LIST, "$sample_list");
while (<LIST>){    
    next unless /\S/; 
    @list = split(/\t/, $_);
    push(@sample, $list[0]);
    push(@sample, $list[1]);
}

#collect status of SID jobs per sample per chunk
my $i = 0;
foreach my $sample_norm (@sample){
    $i ++;
    next if $i % 2 != 1;
    my $sample_tumor = "$sample[$i]";
    my $sample_pair = "$sample_norm.vs.$sample_tumor";

    foreach my $chunk (1..$total_chunks){
	my $chunk_formatted=`printf "%.3d\n" $chunk`;
	chomp($chunk_formatted);
	my $indels_log = "$project_dir_analysis/$date/$sample_pair/run/SI_$sample_pair"."_$chunk_formatted.log";

	$chunk = "chunk_$chunk";
	if (-s $indels_log){
	    my $indels_grep = "";
	    $indels_grep = `grep -vP '^WARN' $indels_log|grep ERROR`;
	    my $indels_file = "$project_dir_analysis/$date/$sample_pair/SomaticIndelDetector/$sample_pair.$chunk.vcf";
	    if (-s $indels_file && $indels_grep eq ""){
		$sum{$sample_pair}{$chunk}{'indels_chunk'} = "PASS";
	    }else{
		$sum{$sample_pair}{$chunk}{'indels_chunk'} = "FAIL";
	    }
	}
    }
    my $merge_log = "$project_dir_analysis/$date/$sample_pair/run/MI_$sample_pair"."_000.log";
    if (-s $merge_log){
	my $merge_grep = "";
	$merge_grep = `grep ERROR $merge_log`;
	my $merge_file = "$project_dir_results/$date/$sample_pair/$sample_pair.vcf";
	if (-s $merge_file && $merge_grep eq ""){
	    $sum{$sample_pair}{'merge_indels'} = "PASS";
	}else{
	    $sum{$sample_pair}{'merge_indels'} = "FAIL";
	}
    }
}

#print IndelGenotyper status and variant stats
#############
open (OUT, ">$summary_results/index.html");
print OUT "<HTML>";
print OUT "<HEAD><META HTTP-EQUIV='refresh' CONTENT='60'></HEAD>";
print OUT "<TABLE CELLPADDING=5><TR>";
print OUT "<TH><CENTER>Sample";
print OUT "<TH><CENTER>Chunk";
print OUT "<TH><CENTER>SomaticIndelDetector";

$i = 0;
foreach my $sample_norm (@sample){
    $i++;
    next if $i % 2 != 1;
    my $sample_tumor = "$sample[$i]";
    my $sample_pair = "$sample_norm.vs.$sample_tumor";
    foreach my $chunk (1..$total_chunks){
	if ($chunk == 1){
	    print OUT "<TR><TD>$sample_pair<TD><CENTER>$chunk\n";
	}else{
	    print OUT "<TR><TD> <TD><CENTER>$chunk\n";
	}

	my $chunk_formatted=`printf "%.3d\n" $chunk`;
	chomp($chunk_formatted);
	my $indels_log = "$project_dir_analysis/$date/$sample_pair/run/SI_$sample_pair"."_$chunk_formatted.log";

        $chunk = "chunk_$chunk";
	if (-s $indels_log){
	    print OUT "<TD><CENTER><IMG SRC=tick.png ALT=PASS>\n" if $sum{$sample_pair}{$chunk}{'indels_chunk'} eq "PASS";
	    print OUT "<TD><CENTER><IMG SRC=error.png ALT=FAIL>\n" if $sum{$sample_pair}{$chunk}{'indels_chunk'} eq "FAIL";
	}
    }
    print OUT "<TR><TD> <TD><CENTER>merged";

    my $merge_log = "$project_dir_analysis/$date/$sample_pair/run/MI_$sample_pair"."_000.log";
    if (-s $merge_log){
	if ($sum{$sample_pair}{'merge_indels'} eq "PASS"){
	    my $merge_file = "$project_dir_results/$date/$sample_pair/$sample_pair.stats.somatic";
	    my $merge_file_easy = "$project_dir_results/$date/$sample_pair/$sample_pair.stats.somatic.easy";

	    open (EASY,">$merge_file_easy");
	    print EASY "\tCONTIG\tSTART_POS\tEND_POS\tALT_ALLELE\t";
	    my $title = `head -n 1 $merge_file`;
	    chomp($title);
	    my @title = split(/\t/, $title);
	    foreach my $column (@title){
	        if ($column =~ /(.*):(.*)/){
		    print EASY "$1\t";
		}
	    }
	    print EASY "JUDGEMENT\n";

	    open (STATS,"$merge_file");
	    while(<STATS>){
		chomp;
		print EASY "\t";
		my @data = split(/\t/, $_);
		foreach my $data (@data){
		    if ($data =~ /(.*):(.*)/){
			print EASY "$2\t";
		    }else {
			print EASY "$data\t";
		    }
		}
		print EASY "\n";
	    }

	    my $php_file = "$project_dir_results/$date/$sample_pair/$sample_pair.stats.somatic.php";
            my $indels_deployment = "$summary_deployment/stats";
            system("ssh $deployment_server mkdir -p -m 775 $indels_deployment > /dev/null 2>&1");
            system("scp -r $merge_file.easy $deployment_server:$indels_deployment/$sample_pair.stats.somatic > /dev/null 2>&1");
            system("scp -r $php_file $deployment_server:$indels_deployment/$sample_pair.stats.somatic.php > /dev/null 2>&1");
            print OUT "<TD><CENTER><A HREF = http://$deployment_server/$summary_link/stats/$sample_pair.stats.somatic.php><IMG BORDER='5' SRC=tick.png ALT=PASS></A>\n";
	}elsif ($sum{$sample_pair}{'merge_stats'} eq "FAIL"){
	    print OUT "<TD><CENTER><IMG SRC=error.png ALT=FAIL>\n"
	}
    }
} 

print OUT "</TABLE>\n";

#deploy summary
system("scp -r $summary_results/index.html $deployment_server:$summary_deployment/index.html > /dev/null 2>&1");
system("ssh $deployment_server chmod -R 775 $summary_deployment/* > /dev/null 2>&1");






