#!/usr/bin/perl -w

#This script generates on-line summary that reflects jobs progress. 
#Once job is finished and non-empty log file is present in the run directory, 
#its content is tested for the presence of error messages from gatk. 
#If (i) there are no error messages in the log file and 
#(ii) expected output file is not empty then green tick appears in the summary table.
#If at least one of these conditions is not satisfied then red cross appears in the summary table.

$project_dir_analysis = "projectDirAnalysis";
$project_dir_results = "projectDirResults";
$project = "#project";
$date = "#today";
$summary_results = "summaryResults";
$deployment_server = "deploymentServer";
$summary_deployment = "summaryDeployment";
$data_deployment = "dataDeployment";
$sample_list = "sampleList";
$ref_chunks = "refChunks";
$ref_intervals = "refIntervals";
$type = "#type";
$summary_link = $1 if  $summary_deployment =~ /\/www\/html\/(.*)/;

%sum = ();

open (CHUNKS, "$ref_chunks");
$total_chunks = 0;
while (<CHUNKS>){
    chomp();
    $chunk = $1 if /^\S+\t\S+\t\S+\t\S+\t(\S+)/;
    $total_chunks = $chunk if $chunk > $total_chunks;
}

if ($type eq "TARGETED" && $ref_intervals ne "null"){
    open (INTERVALS, "$ref_intervals");
    $total_intervals = 0;
    while (<INTERVALS>){
        chomp();
        $interval = $1 if /^\S+\t\S+\t\S+\t\S+\t(\S+)/;
        $total_intervals = $interval if $interval > $total_intervals;
    }
}

@sample = (); @list = ();
open (LIST, "$sample_list");
while (<LIST>){    
    @list = split(/\t/, $_);
    push(@sample, $list[0]);
    push(@sample, $list[2]);
}

#collect data from log files
foreach $sample (@sample){
    foreach $chunk (1..$total_chunks){
	$realignment_log = "$project_dir_analysis/$date/$sample/run/gatk_realign_recal.$sample".".chunk_"."$chunk".".log";
	if (-s $realignment_log){
	    $realignment_grep = "";
	    $realignment_grep = `grep ERROR $realignment_log`;
   	    $realignment_target = "$project_dir_analysis/$date/$sample/realignment/$sample".".chunk_"."$chunk".".RTC.intervals";
	    if (-s $realignment_target && $realignment_grep eq ""){
                $sum{$sample}{$chunk}{'realignment_target'} = "PASS";
	    }else{
                $sum{$sample}{$chunk}{'realignment_target'} = "FAIL";
	    }

   	    $realigned_bam = "$project_dir_analysis/$date/$sample/realignment/$sample".".chunk_"."$chunk".".realigned.bam";
	    if (-s $realigned_bam && $realignment_grep eq ""){
                $sum{$sample}{$chunk}{'realigned_bam'} = "PASS";
	    }else{
                $sum{$sample}{$chunk}{'realigned_bam'} = "FAIL";
	    }

	    $recalibration_report_chunk = "$project_dir_analysis/$date/$sample/recalibration/reports/pre/$sample".".chunk_"."$chunk".".realigned.recal_data.grp";
	    if (-s $recalibration_report_chunk && $realignment_grep eq ""){
                $sum{$sample}{$chunk}{'recalibration_report_chunk'} = "PASS";
	    }else{
                $sum{$sample}{$chunk}{'recalibration_report_chunk'} = "FAIL";
	    }
	}

        $recalibration_log = "$project_dir_analysis/$date/$sample/run/gatk_print_reads.$sample".".chunk_"."$chunk".".log";
        if (-s $recalibration_log){
	    $recalibrated_grep = "";
	    $recalibrated_grep = `grep ERROR $recalibration_log`;
	    $recalibrated_bam = "$project_dir_analysis/$date/$sample/recalibration/$sample".".chunk_"."$chunk.realigned.recalibrated.bam";
	    if (-s $recalibrated_bam && $recalibrated_grep eq ""){
                $sum{$sample}{$chunk}{'recalibrated_bam'} = "PASS";
	    }else{
                $sum{$sample}{$chunk}{'recalibrated_bam'} = "FAIL";
	    }
	    $reduced_bam = "$project_dir_analysis/$date/$sample/recalibration/$sample".".chunk_"."$chunk.realigned.recalibrated.reduced.bam";
	    if (-s $reduced_bam && $recalibrated_grep eq ""){
                $sum{$sample}{$chunk}{'reduced_bam'} = "PASS";
	    }else{
                $sum{$sample}{$chunk}{'reduced_bam'} = "FAIL";
	    }
        }
    }

    $recalibration_merge_log = "$project_dir_analysis/$date/$sample/run/gatk_merge_recalibration_reports.$sample".".log";
    if (-s $recalibration_merge_log){
	$recalibration_merge_grep = "";
	$recalibration_merge_grep = `grep ERROR $recalibration_merge_log`;
        $recalibration_report = "$project_dir_results/$date/$sample/recalibration/reports/pre/$sample".".realigned.recal_data.grp";
	if (-s $recalibration_report && $recalibration_merge_grep eq ""){
            $sum{$sample}{'recalibration_report'} = "PASS";
	}else{
            $sum{$sample}{'recalibration_report'} = "FAIL";
	}
    }
}

$i = 0;
foreach $sample_norm (@sample){
    next if $i % 2 == 1;
    $sample_tumor = "$sample[$i+1]";
    $sample_pair = "$sample_norm.vs.$sample_tumor";
    foreach $chunk (1..$total_chunks){
    $chunk = "chunk_$chunk";
	$mutect_log = "$project_dir_analysis/$date/$sample_pair/run/mutect.$sample_pair.$chunk.log";
	if (-s $mutect_log){
	    $mutect_grep = "";
	    $mutect_grep = `grep ERROR $mutect_log`;
	    $mutect_stats = "$project_dir_analysis/$date/$sample_pair/mutect/$sample_pair.$chunk.stats";
	    if (-s $mutect_stats && $mutect_grep eq ""){
		$sum{$sample_pair}{$chunk}{'mutect_stats'} = "PASS";
	    }else{
		$sum{$sample_pair}{$chunk}{'mutect_stats'} = "FAIL";
	    }
	}
    }
    $merge_log = "$project_dir_analysis/$date/$sample_pair/run/merge_stats.$sample_pair.log";
    if (-s $merge_log){
	$merge_grep = "";
	$merge_grep = `grep ERROR $merge_log`;
	$merge_stats = "$project_dir_results/$date/$sample_pair/$sample_pair.stats";
	if (-s $merge_stats && $merge_grep eq ""){
	    $sum{$sample_pair}{'merge_stats'} = "PASS";
	}else{
	    $sum{$sample_pair}{'merge_stats'} = "FAIL";
	}
    }
    $i ++;
}

#print extracted data into summary file
open (OUT, ">$summary_results/index.html");
print OUT "<HTML>";
print OUT "<HEAD><META HTTP-EQUIV='refresh' CONTENT='60'></HEAD>";
print OUT "<TABLE CELLPADDING=5><TR>";
print OUT "<TH><CENTER>Sample";
print OUT "<TH><CENTER>Chunk";
print OUT "<TH><CENTER>Realignment<BR>intervals";
print OUT "<TH><CENTER>Realigned bam";
print OUT "<TH><CENTER>Recalibration report<BR>chunk";
print OUT "<TH><CENTER>Recalibration report<BR>full";
print OUT "<TH><CENTER>Recalibrated bam";
print OUT "<TH><CENTER>Reduced bam";

@tags = qw(realignment_target realigned_bam recalibration_report_chunk recalibration_report recalibrated_bam reduced_bam);
foreach $sample (@sample){
	
	foreach $chunk (1..$total_chunks){
	
		if ($chunk == 1){
			print OUT "<TR><TD>$sample<TD><CENTER>$chunk\n";
		}else{
			print OUT "<TR><TD> <TD><CENTER>$chunk\n";
		}
        
		foreach $tag (@tags){
		
			if ($tag eq "recalibration_report"){

				if ($chunk == 1){
					
					$recalibration_merge_plot_log = "$project_dir_analysis/$date/$sample/run/gatk_collect_summary_metrics.$sample".".log";
					
					if (-s $recalibration_merge_plot_log){
						
						$pdf = "$project_dir_results/$date/$sample/recalibration/plots/post/$sample.realigned.recalibrated.recal_plot-1.jpeg";
						$pdf_deployment = "$summary_deployment/$sample.jpeg";
						
						if (-s $pdf){

							system("scp -r $pdf $deployment_server:$pdf_deployment > /dev/null 2>&1");
							print OUT "<TD><CENTER><A HREF = http://$deployment_server/$summary_link/$sample.jpeg><IMG BORDER='5' SRC=tick.png ALT=PASS></A>\n" if $sum{$sample}{$tag} eq "PASS";
							print OUT "<TD><CENTER><A HREF = http://$deployment_server/$summary_link/$sample.jpeg><IMG BORDER='5' SRC=error.png ALT=FAIL></A>\n" if $sum{$sample}{$tag} eq "FAIL";
							print OUT "<TD> \n" if $sum{$sample}{$tag} ne "FAIL" && $sum{$sample}{$tag} ne "PASS";
						
						}else{
							
							print OUT "<TD><CENTER><IMG SRC=tick.png ALT=PASS></A>\n" if $sum{$sample}{$tag} eq "PASS";
							print OUT "<TD><CENTER><IMG SRC=error.png ALT=FAIL></A>\n" if $sum{$sample}{$tag} eq "FAIL";
							print OUT "<TD> \n" if $sum{$sample}{$tag} ne "FAIL" && $sum{$sample}{$tag} ne "PASS";
		        
						}
						
						
					}else{
		                
						print OUT "<TD> \n";
			    
					}
					#end of if (-s $recalibration_merge_plot_log){
						
				}else{
		            
					print OUT "<TD> \n";
		        
				}
				#end of if ($chunk == 1){
				
			}else{

			        if (defined($sum{$sample}{$chunk}{$tag})){
				    print OUT "<TD><CENTER><IMG SRC=tick.png ALT=PASS>\n" if $sum{$sample}{$chunk}{$tag} eq "PASS";
				    print OUT "<TD><CENTER><IMG SRC=error.png ALT=FAIL>\n" if $sum{$sample}{$chunk}{$tag} eq "FAIL";
				}else{
				    print OUT "<TD> \n";
				}
		    
			}
			#end of if ($tag eq "recalibration_report"){

		}
		# end of foreach $tag (@tags){

	}
	#end of foreach $chunk (1..$total_chunks){

}
#end of foreach $sample (@sample){

print OUT "</TABLE>\n";


#horizontal line
print OUT "<HR>";


#print mutect status and variant stats
#############
print OUT "<TABLE CELLPADDING=5><TR>";
print OUT "<TH><CENTER>Sample";
print OUT "<TH><CENTER>Chunk";
print OUT "<TH><CENTER>MuTect run";

$i = 0;
foreach $sample_norm (@sample){
    next if $i % 2 == 1;
    $sample_tumor = "$sample[$i+1]";
    $sample_pair = "$sample_norm.vs.$sample_tumor";
    foreach $chunk (1..$total_chunks){
	if ($chunk == 1){
	    print OUT "<TR><TD>$sample_pair<TD><CENTER>$chunk\n";
	}else{
	    print OUT "<TR><TD> <TD><CENTER>$chunk\n";
	}
        $chunk = "chunk_$chunk";
	$mutect_log = "$project_dir_analysis/$date/$sample_pair/run/mutect.$sample_pair.$chunk.log";
	if (-s $mutect_log){
	    print OUT "<TD><CENTER><IMG SRC=tick.png ALT=PASS>\n" if $sum{$sample_pair}{$chunk}{'mutect_stats'} eq "PASS";
	    print OUT "<TD><CENTER><IMG SRC=error.png ALT=FAIL>\n" if $sum{$sample_pair}{$chunk}{'mutect_stats'} eq "FAIL";
	}
    }
    print OUT "<TR><TD> <TD><CENTER>merged";
    $merge_log = "$project_dir_analysis/$date/$sample_pair/run/merge_stats.$sample_pair.log";
    if (-s $merge_log){
	if ($sum{$sample_pair}{'merge_stats'} eq "PASS"){
            $merge_stats = "$project_dir_results/$date/$sample_pair/$sample_pair.stats.keep";
            $stats_deployment = "$summary_deployment/stats";
            system("ssh $deployment_server mkdir -p -m 775 $stats_deployment > /dev/null 2>&1");
            system("scp -r $merge_stats $deployment_server:$stats_deployment/$sample_pair.stats.keep > /dev/null 2>&1");
            system("scp -r $merge_stats.php $deployment_server:$stats_deployment/$sample_pair.stats.keep.php > /dev/null 2>&1");
            print OUT "<TD><CENTER><A HREF = http://$deployment_server/$summary_link/stats/$sample_pair.stats.keep.php><IMG BORDER='5' SRC=tick.png ALT=PASS></A>\n";
	}elsif ($sum{$sample_pair}{'merge_stats'} eq "FAIL"){
	    print OUT "<TD><CENTER><IMG SRC=error.png ALT=FAIL>\n"
	}
    }
    $i++;
} 

print OUT "</TABLE>\n";

#summary metrics
################

#horizontal line
print OUT "<HR>";

#deploy files
$sample_summary = "$project_dir_results/$date/multisample/metrics/$project.$date.sample_summary";
$sample_interval_summary = "$project_dir_results/$date/multisample/metrics/$project.$date.sample_interval_summary";
$sample_cumulative_coverage_proportions = "$project_dir_results/$date/multisample/metrics/$project.$date.sample_cumulative_coverage_proportions";

$sample_summary_php = "$project_dir_analysis/$date/multisample/run/sample_summary.php";
$sample_interval_summary_php = "$project_dir_analysis/$date/multisample/run/sample_interval_summary.php";
$sample_cumulative_coverage_proportions_php = "$project_dir_analysis/$date/multisample/run/sample_cumulative_coverage_proportions.php";

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

$sample_summary_url = "http://$deployment_server/$summary_link/metrics/$project.$date.sample_summary";
$sample_summary_php_url = "http://$deployment_server/$summary_link/metrics/sample_summary.php";

$sample_interval_summary_url = "http://$deployment_server/$summary_link/metrics/$project.$date.sample_interval_summary";
$sample_interval_summary_php_url = "http://$deployment_server/$summary_link/metrics/sample_interval_summary.php";

$sample_cumulative_coverage_proportions_url = "http://$deployment_server/$summary_link/metrics/$project.$date.sample_cumulative_coverage_proportions";
$sample_cumulative_coverage_proportions_php_url = "http://$deployment_server/$summary_link/metrics/sample_cumulative_coverage_proportions.php";

if (-s $sample_summary){
print OUT "<P><FONT SIZE = '+1'><A HREF = '$sample_summary_php_url'>Coverage Summary</A> ".
							   "<A HREF = '$sample_summary_url'>[TSV]</A> ".
							   "<A HREF = '$sample_summary_url.xlsx'>[XLS]</A>".
			  "</FONT><BR>";
}

if (-s $sample_interval_summary){
print OUT "<P><FONT SIZE = '+1'><A HREF = '$sample_interval_summary_php_url'>Interval Coverage Summary</A> ".
							   "<A HREF = '$sample_interval_summary_url'>[TSV]</A> ".
							   "<A HREF = '$sample_interval_summary_url.xlsx'>[XLS]</A> ".
							   "<A HREF = '$sample_interval_summary_url.png'>[Plot]</A>".
			"</FONT><BR>";
}

if (-s $sample_cumulative_coverage_proportions){			
print OUT "<P><FONT SIZE = '+1'><A HREF = '$sample_cumulative_coverage_proportions_php_url'>Cumulative Coverage Proportions</A> ".
							   "<A HREF = '$sample_cumulative_coverage_proportions_url'>[TSV]</A> ".
							   "<A HREF = '$sample_cumulative_coverage_proportions_url.xlsx'>[XLS]</A> ".
							   "<A HREF = '$sample_cumulative_coverage_proportions_url.png'>[Plot]</A> ".
			"</FONT><BR>";
}

#usage statistics
#################

#horizontal line
print OUT "<HR>";

$usage_file = "$project_dir_results/$date/multisample/usage.$date.txt";
$usage_url = "http://$deployment_server/$summary_link/usage.txt";
if (-s $usage_file){
    system("scp -r $usage_file $deployment_server:$summary_deployment/usage.txt > /dev/null 2>&1");
    print OUT "<HR><P><FONT SIZE = '+1'>Usage of computational resources can be monitored <A HREF = '$usage_url'>here</A></FONT><BR>";
}

system("scp -r $summary_results/index.html $deployment_server:$summary_deployment/index.html > /dev/null 2>&1");
system("ssh $deployment_server chmod -R 775 $summary_deployment/* > /dev/null 2>&1");

