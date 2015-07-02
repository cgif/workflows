#!/usr/bin/perl -w

use strict;
use warnings;

#This script generates on-line summary that reflects jobs progress. 
#Once job is finished and non-empty log file is present in the run directory, 
#its content is tested for the presence of error messages from gatk. 
#If (i) there are no error messages in the log file and 
#(ii) expected output file is not empty then green tick appears in the summary table.
#If at least one of these conditions is not satisfied then red cross appears in the summary table.

my $project_dir_analysis = "projectDirAnalysis";
my $project_dir_results = "projectDirResults";
my $project = "#project";
my $date = "#today";
my $summary_results = "summaryResults";
my $deployment_server = "deploymentServer";
my $summary_deployment = "summaryDeployment";
my $sample_list = "sampleList";
my $ref_chunks = "refChunks";
my $ref_intervals = "refIntervals";
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
    @list = split(/\t/, $_);
    push(@sample, $list[0]);
    push(@sample, $list[2]);
}

#collect data from gatk log files for each sample for each chunk 
foreach my $sample (@sample){
    foreach my $chunk (1..$total_chunks){
	my $formatted_chunk = `printf "%.3d\n" $chunk`;
	chomp($formatted_chunk);
	my $realignment_log = "$project_dir_analysis/$date/$sample/run/RR_$sample"."_$formatted_chunk.log";
	if (-s $realignment_log){
	    my $realignment_grep = "";
	    $realignment_grep = `grep -P 'ERROR|WARNING' $realignment_log`;

	    if ($realignment_grep eq ""){
                $sum{$sample}{$chunk}{'realigned_bam'} = "PASS";
	    }else{
                $sum{$sample}{$chunk}{'realigned_bam'} = "FAIL";
	    }

	    my $recalibration_report_chunk = "$project_dir_analysis/$date/$sample/recalibration/reports/pre/$sample.chunk_$chunk.realigned.recal_data.grp";
	    if (-s $recalibration_report_chunk && $realignment_grep eq ""){
                $sum{$sample}{$chunk}{'recalibration_report_chunk'} = "PASS";
	    }else{
                $sum{$sample}{$chunk}{'recalibration_report_chunk'} = "FAIL";
	    }
	}

        my $recalibration_log = "$project_dir_analysis/$date/$sample/run/PR_$sample"."_$formatted_chunk.log";
        if (-s $recalibration_log){
	    my $recalibrated_grep = "";
	    $recalibrated_grep = `grep -P 'ERROR|WARNING' $recalibration_log`;

	    if ($recalibrated_grep eq ""){
                $sum{$sample}{$chunk}{'recalibrated_bam'} = "PASS";
	    }else{
                $sum{$sample}{$chunk}{'recalibrated_bam'} = "FAIL";
	    }
        }
    }

    my $recalibration_merge_log = "$project_dir_analysis/$date/$sample/run/CS_$sample"."_000.log";
    if (-s $recalibration_merge_log){
	my $recalibration_merge_grep = "";
	$recalibration_merge_grep = `grep ERROR $recalibration_merge_log`;
        my $recalibration_report = "$project_dir_results/$date/$sample/recalibration/reports/pre/$sample.realigned.recal_data.grp";
	if (-s $recalibration_report && $recalibration_merge_grep eq ""){
            $sum{$sample}{'recalibration_report'} = "PASS";
	}else{
            $sum{$sample}{'recalibration_report'} = "FAIL";
	}
    }

    my $merge_bam_log = "$project_dir_analysis/$date/$sample/run/MB_$sample"."_000.log";
    if (-s $merge_bam_log){
	my $merge_bam_grep = "";
	$merge_bam_grep = `grep -vP '^gatk' $merge_bam_log | grep ERROR`;
	my $merge_bam_file = "$project_dir_results/$date/$sample/recalibration/$sample.bam";
	if (-s $merge_bam_file && $merge_bam_grep eq ""){
            $sum{$sample}{'merge_bam'} = "PASS";
	}else{
            $sum{$sample}{'merge_bam'} = "FAIL";
	}
    }
}

#collect data from mutect log files for somatic variant calling
my $i = 0;
foreach my $sample_norm (@sample){
    $i ++;
    next if $i % 2 != 1;
    my $sample_tumor = "$sample[$i]";
    my $sample_pair = "$sample_norm.vs.$sample_tumor";
    foreach my $chunk (1..$total_chunks){
	my $chunk_formatted=`printf "%.3d\n" $chunk`;
	chomp($chunk_formatted);
	my $mutect_log = "$project_dir_analysis/$date/$sample_pair/run/MU_$sample_pair"."_$chunk_formatted.log";
	$chunk = "chunk_$chunk";
	if (-s $mutect_log){
	    my $mutect_grep = "";
	    $mutect_grep = `grep -vP '^WARN' $mutect_log | grep ERROR`;
	    my $mutect_report = "$project_dir_analysis/$date/$sample_pair/mutect/$sample_pair.$chunk.stats";
	    if (-s $mutect_report && $mutect_grep eq ""){
		$sum{$sample_pair}{$chunk}{'mutect_stats'} = "PASS";
	    }else{
		$sum{$sample_pair}{$chunk}{'mutect_stats'} = "FAIL";
	    }
	}
    }
    my $merge_log = "$project_dir_analysis/$date/$sample_pair/run/MS_$sample_pair"."_000.log";
    if (-s $merge_log){
	my $merge_grep = "";
	$merge_grep = `grep ERROR $merge_log`;
	my $merge_stats = "$project_dir_results/$date/$sample_pair/$sample_pair.stats";
	if (-s $merge_stats && $merge_grep eq ""){
	    $sum{$sample_pair}{'merge_stats'} = "PASS";
	}else{
	    $sum{$sample_pair}{'merge_stats'} = "FAIL";
	}
    }
}

#print extracted data into summary file
#print data from gatk log files for each sample for each chunk 
##############################################################
open (OUT, ">$summary_results/index.html");
print OUT "<HTML>";
print OUT "<HEAD><META HTTP-EQUIV='refresh' CONTENT='60'></HEAD>";
print OUT "<TABLE CELLPADDING=5><TR>";
print OUT "<TH><CENTER>Sample";
print OUT "<TH><CENTER>Chunk";
print OUT "<TH><CENTER>Realigned bam";
print OUT "<TH><CENTER>Recalibration report<BR>chunk";
print OUT "<TH><CENTER>Recalibration report<BR>full";
print OUT "<TH><CENTER>Recalibrated bam";
print OUT "<TH><CENTER>Merged bam";

my @tags = qw(realigned_bam recalibration_report_chunk recalibration_report recalibrated_bam merge_bam);
foreach my $sample (@sample){
	
	foreach my $chunk (1..$total_chunks){
	
		if ($chunk == 1){
			print OUT "<TR><TD>$sample<TD><CENTER>$chunk\n";
		}else{
			print OUT "<TR><TD> <TD><CENTER>$chunk\n";
		}
        
		foreach my $tag (@tags){
		
			if ($tag eq "recalibration_report"){

				if ($chunk == 1){
					
					my $recalibration_merge_plot_log = "$project_dir_analysis/$date/$sample/run/CS_$sample"."_000.log";
					
					if (-s $recalibration_merge_plot_log){
						
						my $pdf = "$project_dir_results/$date/$sample/recalibration/plots/post/$sample.realigned.recalibrated.recal_plot-1.jpeg";
						my $pdf_deployment = "$summary_deployment/$sample.jpeg";
						
						if (-s $pdf){

							system("scp -r $pdf $deployment_server:$pdf_deployment > /dev/null 2>&1");
							print OUT "<TD><CENTER><A HREF = http://$deployment_server/$summary_link/$sample.jpeg><IMG BORDER='5' SRC=tick.png ALT=PASS></A>\n" if $sum{$sample}{$tag} eq "PASS";
							print OUT "<TD><CENTER><A HREF = http://$deployment_server/$summary_link/$sample.jpeg><IMG BORDER='5' SRC=error.png ALT=FAIL></A>\n" if $sum{$sample}{$tag} eq "FAIL";
						
						}else{
							
							print OUT "<TD><CENTER><IMG SRC=error.png ALT=FAIL></A>\n";
		        
						}
						
						
					}else{
		                
						print OUT "<TD> \n";
			    
					} 
						
				}else{
		            
					print OUT "<TD> \n";
		        
				} 

			}elsif ($tag eq "merge_bam"){

				if ($chunk == 1){

				        my $merge_bam_log = "$project_dir_analysis/$date/$sample/run/MB_$sample"."_000.log";

					if (-s $merge_bam_log){

					        print OUT "<TD><CENTER><IMG SRC=tick.png ALT=PASS>\n" if $sum{$sample}{$tag} eq "PASS";
						print OUT "<TD><CENTER><IMG SRC=error.png ALT=FAIL>\n" if $sum{$sample}{$tag} eq "FAIL";

					}

				}else{

				        print OUT "<TD> \n";

				}
				
			}else{

			        if (defined($sum{$sample}{$chunk}{$tag})){

				    print OUT "<TD><CENTER><IMG SRC=tick.png ALT=PASS>\n" if $sum{$sample}{$chunk}{$tag} eq "PASS";
				    print OUT "<TD><CENTER><IMG SRC=error.png ALT=FAIL>\n" if $sum{$sample}{$chunk}{$tag} eq "FAIL";

				}else{

				    print OUT "<TD> \n";

				}
		    
			} #end of if ($tag eq ...)

		} #end of foreach $tag (@tags)

	} #end of foreach $chunk (1..$total_chunks)

} #end of foreach $sample (@sample)

print OUT "</TABLE>\n";


#horizontal line
print OUT "<HR>";


#print mutect status and deploy variant stats
##########################################################
print OUT "<TABLE CELLPADDING=5><TR>";
print OUT "<TH><CENTER>Sample";
print OUT "<TH><CENTER>Chunk";
print OUT "<TH><CENTER>MuTect run";

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
	my $formatted_chunk =`printf "%.3d\n" $chunk`;
	chomp($formatted_chunk);
	my $mutect_log = "$project_dir_analysis/$date/$sample_pair/run/MU_$sample_pair"."_$formatted_chunk.log";

        $chunk = "chunk_$chunk";
	if (-s $mutect_log){
	    print OUT "<TD><CENTER><IMG SRC=tick.png ALT=PASS>\n" if $sum{$sample_pair}{$chunk}{'mutect_stats'} eq "PASS";
	    print OUT "<TD><CENTER><IMG SRC=error.png ALT=FAIL>\n" if $sum{$sample_pair}{$chunk}{'mutect_stats'} eq "FAIL";
	}
    }
    print OUT "<TR><TD> <TD><CENTER>merged";
    my $merge_log = "$project_dir_analysis/$date/$sample_pair/run/MS_$sample_pair"."_000.log";
    if (-s $merge_log){
	if ($sum{$sample_pair}{'merge_stats'} eq "PASS"){
            my $merge_stats = "$project_dir_results/$date/$sample_pair/$sample_pair.stats.keep";
            my $stats_deployment = "$summary_deployment/stats";
            system("ssh $deployment_server mkdir -p -m 775 $stats_deployment > /dev/null 2>&1");
            system("scp -r $merge_stats $deployment_server:$stats_deployment/$sample_pair.stats.keep > /dev/null 2>&1");
            system("scp -r $merge_stats.php $deployment_server:$stats_deployment/$sample_pair.stats.keep.php > /dev/null 2>&1");
            print OUT "<TD><CENTER><A HREF = http://$deployment_server/$summary_link/stats/$sample_pair.stats.keep.php><IMG BORDER='5' SRC=tick.png ALT=PASS></A>\n";
	}elsif ($sum{$sample_pair}{'merge_stats'} eq "FAIL"){
	    print OUT "<TD><CENTER><IMG SRC=error.png ALT=FAIL>\n"
	}
    }
} 

print OUT "</TABLE>\n";


#deploy summary metrics and print links on summary page
############################################################

#horizontal line
print OUT "<HR>";

#deploy files
my $sample_summary = "$project_dir_results/$date/multisample/metrics/$project.$date.sample_summary";
my $sample_interval_summary = "$project_dir_results/$date/multisample/metrics/$project.$date.sample_interval_summary";
my $sample_cumulative_coverage_proportions = "$project_dir_results/$date/multisample/metrics/$project.$date.sample_cumulative_coverage_proportions";

my $sample_summary_php = "$project_dir_analysis/$date/multisample/run/sample_summary.php";
my $sample_interval_summary_php = "$project_dir_analysis/$date/multisample/run/sample_interval_summary.php";
my $sample_cumulative_coverage_proportions_php = "$project_dir_analysis/$date/multisample/run/sample_cumulative_coverage_proportions.php";

#create target directory on server
system("ssh $deployment_server mkdir -p -m 775 $summary_deployment/metrics > /dev/null 2>&1");

#copy files to server
system("scp -r $sample_summary $deployment_server:$summary_deployment/metrics > /dev/null 2>&1");
system("scp -r $sample_summary_php $deployment_server:$summary_deployment/metrics > /dev/null 2>&1");
system("scp -r $sample_summary.xlsx $deployment_server:$summary_deployment/metrics > /dev/null 2>&1");

system("scp -r $sample_interval_summary $deployment_server:$summary_deployment/metrics > /dev/null 2>&1");
system("scp -r $sample_interval_summary_php $deployment_server:$summary_deployment/metrics > /dev/null 2>&1");
system("scp -r $sample_interval_summary.png $deployment_server:$summary_deployment/metrics > /dev/null 2>&1");
system("scp -r $sample_interval_summary.xlsx $deployment_server:$summary_deployment/metrics > /dev/null 2>&1");

system("scp -r $sample_cumulative_coverage_proportions $deployment_server:$summary_deployment/metrics > /dev/null 2>&1");
system("scp -r $sample_cumulative_coverage_proportions_php $deployment_server:$summary_deployment/metrics > /dev/null 2>&1");
system("scp -r $sample_cumulative_coverage_proportions.png $deployment_server:$summary_deployment/metrics > /dev/null 2>&1");
system("scp -r $sample_cumulative_coverage_proportions.xlsx $deployment_server:$summary_deployment/metrics > /dev/null 2>&1");

#print links
my $sample_summary_url = "http://$deployment_server/$summary_link/metrics/$project.$date.sample_summary";
my $sample_summary_php_url = "http://$deployment_server/$summary_link/metrics/sample_summary.php";

my $sample_interval_summary_url = "http://$deployment_server/$summary_link/metrics/$project.$date.sample_interval_summary";
my $sample_interval_summary_php_url = "http://$deployment_server/$summary_link/metrics/sample_interval_summary.php";

my $sample_cumulative_coverage_proportions_url = "http://$deployment_server/$summary_link/metrics/$project.$date.sample_cumulative_coverage_proportions";
my $sample_cumulative_coverage_proportions_php_url = "http://$deployment_server/$summary_link/metrics/sample_cumulative_coverage_proportions.php";

if (-s $sample_summary){
print OUT "<P><FONT SIZE = '+1'><A HREF = '$sample_summary_php_url'>Coverage Summary</A> ".
							   "<A HREF = '$sample_summary_url'>[TSV]</A> ".
							   "<A HREF = '$sample_summary_url.xlsx'>[XLS]</A>".
			  "</FONT><BR>";
}

my $plot = "$sample_interval_summary_url".".png";
if (-s $sample_interval_summary){
print OUT "<P><FONT SIZE = '+1'><A HREF = '$sample_interval_summary_php_url'>Interval Coverage Summary</A> ".
							   "<A HREF = '$sample_interval_summary_url'>[TSV]</A> ".
							   "<A HREF = '$sample_interval_summary_url.xlsx'>[XLS]</A> ";
print OUT						   "<A HREF = '$plot'>[Plot]</A>" if (-s $plot);
print OUT "</FONT><BR>";
}

$plot = "$sample_cumulative_coverage_proportions_url".".png";
if (-s $sample_cumulative_coverage_proportions){			
print OUT "<P><FONT SIZE = '+1'><A HREF = '$sample_cumulative_coverage_proportions_php_url'>Cumulative Coverage Proportions</A> ".
							   "<A HREF = '$sample_cumulative_coverage_proportions_url'>[TSV]</A> ".
							   "<A HREF = '$sample_cumulative_coverage_proportions_url.xlsx'>[XLS]</A> ";
print OUT						   "<A HREF = '$plot'>[Plot]</A>" if (-s $plot);
print OUT "</FONT><BR>";
}

#usage statistics
#################

my $usage_file = "$project_dir_results/$date/multisample/usage.$date.txt";
my $usage_php = "$project_dir_analysis/$date/multisample/run/usage.php";
my $usage_url = "http://$deployment_server/$summary_link/usage.php";
if (-s $usage_file){
    system("scp -r $usage_file $deployment_server:$summary_deployment/usage.txt > /dev/null 2>&1");
    system("scp -r $usage_php $deployment_server:$summary_deployment/usage.php > /dev/null 2>&1");
    print OUT "<HR><P><FONT SIZE = '+1'>Usage of computational resources can be monitored <A HREF = '$usage_url'>here</A></FONT><BR>";
}

system("scp -r $summary_results/index.html $deployment_server:$summary_deployment/index.html > /dev/null 2>&1");
system("ssh $deployment_server chmod -R 775 $summary_deployment/* > /dev/null 2>&1");

