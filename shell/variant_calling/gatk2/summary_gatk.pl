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
$variant_caller = "variantCaller";
$enc_dir = "encryptedDir";
$summary_link = $1 if  $summary_deployment =~ /\/www\/html\/(.*)/;
$data_link = $1 if $data_deployment =~ /\/www\/html\/(.*)/;

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

@sample = ();
open (LIST, "$sample_list");
while (<LIST>){
    push(@sample, $1) if /^(\S+)\s/;
}

#collect data from bwa log files
foreach $sample (@sample){
    foreach $chunk (1..$total_chunks){
	$realignment_log = "$project_dir_analysis/$date/$sample/run/gatk_realign_recal.$sample".".chunk_"."$chunk".".log";
	if (-s $realignment_log){
   	    $realignment_target = "$project_dir_analysis/$date/$sample/realignment/$sample".".chunk_"."$chunk".".RTC.intervals";
	    if (-s $realignment_target){
                $sum{$sample}{$chunk}{'realignment_target'} = "PASS";
	    }else{
                $sum{$sample}{$chunk}{'realignment_target'} = "FAIL";
	    }

   	    $realigned_bam = "$project_dir_analysis/$date/$sample/realignment/$sample".".chunk_"."$chunk".".realigned.bam";
	    if (-s $realigned_bam){
                $sum{$sample}{$chunk}{'realigned_bam'} = "PASS";
	    }else{
                $sum{$sample}{$chunk}{'realigned_bam'} = "FAIL";
	    }

	    $recalibration_report_chunk = "$project_dir_analysis/$date/$sample/recalibration/reports/pre/$sample".".chunk_"."$chunk".".realigned.recal_data.grp";
	    if (-s $recalibration_report_chunk){
                $sum{$sample}{$chunk}{'recalibration_report_chunk'} = "PASS";
	    }else{
                $sum{$sample}{$chunk}{'recalibration_report_chunk'} = "FAIL";
	    }
	}

        $recalibration_log = "$project_dir_analysis/$date/$sample/run/gatk_print_reads.$sample".".chunk_"."$chunk".".log";
        if (-s $recalibration_log){
	    $recalibrated_bam = "$project_dir_analysis/$date/$sample/recalibration/$sample".".chunk_"."$chunk.realigned.recalibrated.bam";
	    if (-s $recalibrated_bam){
                $sum{$sample}{$chunk}{'recalibrated_bam'} = "PASS";
	    }else{
                $sum{$sample}{$chunk}{'recalibrated_bam'} = "FAIL";
	    }
	    $reduced_bam = "$project_dir_analysis/$date/$sample/recalibration/$sample".".chunk_"."$chunk.realigned.recalibrated.reduced.bam";
	    if (-s $reduced_bam){
                $sum{$sample}{$chunk}{'reduced_bam'} = "PASS";
	    }else{
                $sum{$sample}{$chunk}{'reduced_bam'} = "FAIL";
	    }
        }
    }

    $recalibration_merge_log = "$project_dir_analysis/$date/$sample/run/gatk_merge_recalibration_reports.$sample".".log";
    if (-s $recalibration_merge_log){
        $recalibration_report = "$project_dir_results/$date/$sample/recalibration/reports/pre/$sample".".realigned.recal_data.grp";
	if (-s $recalibration_report){
            $sum{$sample}{'recalibration_report'} = "PASS";
	}else{
            $sum{$sample}{'recalibration_report'} = "FAIL";
	}
    }
}

foreach $chunk (1..$total_chunks){
    $unified_genotyper_log = "$project_dir_analysis/$date/multisample/run/gatk_unified_genotyper.chunk_"."$chunk".".log";
    if (-s $unified_genotyper_log){
        $unified_genotyper_grep = "";
        $unified_genotyper_grep = `grep ERROR $unified_genotyper_log`;
	$raw_ug_vcf_chunk = "$project_dir_analysis/$date/multisample/unifiedgenotyper/$project"."_multisample.chunk_"."$chunk".".raw.vcf";
	if (-s $raw_ug_vcf_chunk && $unified_genotyper_grep eq ""){
            $sum{'multisample'}{$chunk}{'raw_ug_vcf_chunk'} = "PASS";
	}else{
            $sum{'multisample'}{$chunk}{'raw_ug_vcf_chunk'} = "FAIL";
	}
    }

	if ($type eq "TARGETED" && $ref_intervals ne "null"){
    
		foreach $interval (1..$total_intervals){

			$fragment = "$chunk.interval_$interval";
            $haplotype_caller_log = "$project_dir_analysis/$date/multisample/run/gatk_haplotype_caller.chunk_"."$fragment".".log";

            if (-s $haplotype_caller_log){

                $haplotype_caller_grep = "";
                $haplotype_caller_grep = `grep ERROR $haplotype_caller_log`;
		        $raw_hc_vcf_chunk = "$project_dir_analysis/$date/multisample/haplotypecaller/$project"."_multisample.chunk_"."$fragment".".raw.vcf";

		        if (-s $raw_hc_vcf_chunk && $haplotype_caller_grep eq ""){
	            $sum{'multisample'}{$fragment}{'raw_hc_vcf_chunk'} = "PASS";
	        }else{
                    $sum{'multisample'}{$fragment}{'raw_hc_vcf_chunk'} = "FAIL";
	        }
            }            
	}

    }else{
        $haplotype_caller_log = "$project_dir_analysis/$date/multisample/run/gatk_haplotype_caller.chunk_"."$chunk".".log";
        if (-s $haplotype_caller_log){
            $haplotype_caller_grep = "";
            $haplotype_caller_grep = `grep ERROR $haplotype_caller_log`;
	    $raw_hc_vcf_chunk = "$project_dir_analysis/$date/multisample/haplotypecaller/$project"."_multisample.chunk_"."$chunk".".raw.vcf";
	    if (-s $raw_hc_vcf_chunk && $haplotype_caller_grep eq ""){
                $sum{'multisample'}{$chunk}{'raw_hc_vcf_chunk'} = "PASS";
	    }else{
                $sum{'multisample'}{$chunk}{'raw_hc_vcf_chunk'} = "FAIL";
            }
	}
    }
}

$recalibrate_ug_vcf_log =  "$project_dir_analysis/$date/multisample/run/gatk_recal_ug_vcf.$project"."_multisample.log";
if (-s $recalibrate_ug_vcf_log){
    $raw_ug_vcf = "$project_dir_results/$date/multisample/unifiedgenotyper/$project"."_multisample.raw.vcf.gz";
    if (-s $raw_ug_vcf){
        $sum{'multisample'}{'raw_ug_vcf'} = "PASS"
    }else{
        $sum{'multisample'}{'raw_ug_vcf'} = "FAIL"
    }
    $calibrated_ug_vcf = "$project_dir_results/$date/multisample/unifiedgenotyper/$project"."_multisample.recalibrated.PASS.vcf.gz";
    if (-s $calibrated_ug_vcf){
        $sum{'multisample'}{'calibrated_ug_vcf'} = "PASS"
    }else{
        $sum{'multisample'}{'calibrated_ug_vcf'} = "FAIL"
    }
}

$recalibrate_hc_vcf_log =  "$project_dir_analysis/$date/multisample/run/gatk_recal_hc_vcf.$project"."_multisample.log";
if (-s $recalibrate_hc_vcf_log){
    $raw_hc_vcf = "$project_dir_results/$date/multisample/haplotypecaller/$project"."_multisample.raw.vcf.gz";
    if (-s $raw_hc_vcf){
        $sum{'multisample'}{'raw_hc_vcf'} = "PASS"
    }else{
        $sum{'multisample'}{'raw_hc_vcf'} = "FAIL"
    }
    $calibrated_hc_vcf = "$project_dir_results/$date/multisample/haplotypecaller/$project"."_multisample.recalibrated.PASS.vcf.gz";
    if (-s $calibrated_hc_vcf){
        $sum{'multisample'}{'calibrated_hc_vcf'} = "PASS"
    }else{
        $sum{'multisample'}{'calibrated_hc_vcf'} = "FAIL"
    }
}

%report = ();
foreach $vc (qw(unifiedgenotyper haplotypecaller)){
    $eval_report = "$project_dir_results/$date/multisample/$vc/recalibration/$project"."_multisample.varianteval.report";
    if (-s $eval_report){
        open (REPORT, "$eval_report");
        while (<REPORT>){

            if (/^CompOverlap\s+dbsnp\s+eval\s+raw\s+none\s+all\s+all\s+(\S+)\s+\S+\s+\S+\s+(\S+)\s/){ 
	        $report{$vc}{'raw'}{'total'} = $1;
	        $report{$vc}{'raw'}{'snp'} = $2;
	    }
            if (/^CompOverlap\s+dbsnp\s+eval\s+called\s+none\s+all\s+all\s+(\S+)\s+\S+\s+\S+\s+(\S+)\s/){ 
	        $report{$vc}{'called'}{'total'} = $1;
	        $report{$vc}{'called'}{'snp'} = $2;
	    }
            if (/^CompOverlap\s+dbsnp\s+eval\s+filtered\s+none\s+all\s+all\s+(\S+)\s+\S+\s+\S+\s+(\S+)\s/){ 
	        $report{$vc}{'filtered'}{'total'} = $1;
	        $report{$vc}{'filtered'}{'snp'} = $2;
	    }
	    if (/^CountVariants\s+dbsnp\s+eval\s+called\s+none\s+all\s+all\s+\S+\s+\S+\s+\S+\s+\S+\s+\S+\s+\S+\s+(\S+)\s+\S+\s+(\S+)\s+(\S+)\s/){
	        $report{$vc}{'called'}{'SNP'} = $1;
 	        $report{$vc}{'called'}{'insertion'} = $2;
	        $report{$vc}{'called'}{'deletion'} = $3;
	    }
	    if (/^TiTvVariantEvaluator\s+dbsnp\s+eval\s+raw\s+none\s+all\s+all\s+\S+\s+\S+\s+(\S+)\s/){
	        $report{$vc}{'raw'}{'TiTvall'} = $1;
	    }
	    if (/^TiTvVariantEvaluator\s+dbsnp\s+eval\s+called\s+none\s+(\w+)\s+all\s+\S+\s+\S+\s+(\S+)\s/){
	        $report{$vc}{'called'}{"TiTv$1"} = $2;
	    }
	    if (/^TiTvVariantEvaluator\s+dbsnp\s+eval\s+filtered\s+none\s+all\s+all\s+\S+\s+\S+\s+(\S+)\s/){
	        $report{$vc}{'filtered'}{'TiTvall'} = $1;
	    }
        }
    }
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
		                
						print OUT "<TD><CENTER><IMG SRC=tick.png ALT=PASS></A>\n" if $sum{$sample}{$tag} eq "PASS";
						print OUT "<TD><CENTER><IMG SRC=error.png ALT=FAIL></A>\n" if $sum{$sample}{$tag} eq "FAIL";
						print OUT "<TD> \n" if $sum{$sample}{$tag} ne "FAIL" && $sum{$sample}{$tag} ne "PASS";
			    
					}
					#end of if (-s $recalibration_merge_plot_log){
						
				}else{
		            
					print OUT "<TD> \n";
		        
				}
				#end of if ($chunk == 1){
				
			}else{
	                
				print OUT "<TD><CENTER><IMG SRC=tick.png ALT=PASS>\n" if $sum{$sample}{$chunk}{$tag} eq "PASS";
				print OUT "<TD><CENTER><IMG SRC=error.png ALT=FAIL>\n" if $sum{$sample}{$chunk}{$tag} eq "FAIL";
				print OUT "<TD> \n" if $sum{$sample}{$chunk}{$tag} ne "FAIL" && $sum{$sample}{$chunk}{$tag} ne "PASS";
		    
			}
			#end of if ($tag eq "recalibration_report"){

		}
		# end of foreach $tag (@tags){

	}
	#end of foreach $chunk (1..$total_chunks){

}
#end of foreach $sample (@sample){

print OUT "</TABLE>\n";

#print UnifiedGenotyper statistics
##################################

if ($variant_caller =~ /U/){

#horizontal line
print OUT "<HR>";


#table header
#############
print OUT "<TABLE CELLPADDING=5>";
print OUT "<TR><TH ROWSPAN='2'><CENTER>MULTISAMPLE<BR>UnifiedGenotyper";
print OUT "<TH ROWSPAN='2'><CENTER>Chunk";
print OUT "<TH ROWSPAN='2'><CENTER>Raw vcf<BR>chunk";
print OUT "<TH ROWSPAN='2'><CENTER>Raw vcf<BR>full";
print OUT "<TH ROWSPAN='2'><CENTER>Calibrated vcf" unless $type eq "TARGETED";

print OUT "<TH COLSPAN='3'><CENTER>Raw Variants";
print OUT "<TH COLSPAN='8'><CENTER>Called Variants";
print OUT "<TH COLSPAN='3'><CENTER>Filtered Variants";

print OUT "<TR><TH><CENTER>Total";
print OUT "<TH><CENTER>Known,%";
print OUT "<TH><CENTER>TiTv<BR>total";

print OUT "<TH><CENTER>Total";
print OUT "<TH><CENTER>Known,%";
print OUT "<TH><CENTER>TiTv<BR>total";
print OUT "<TH><CENTER>TiTv<BR>known";
print OUT "<TH><CENTER>TiTv<BR>novel";
print OUT "<TH><CENTER>SNPs<BR>total";
print OUT "<TH><CENTER>Insertions<BR>total";
print OUT "<TH><CENTER>Deletions<BR>total";

print OUT "<TH><CENTER>Total";
print OUT "<TH><CENTER>Known,%";
print OUT "<TH><CENTER>TiTv<BR>total";

#table rows
#variant stats by chunk
#######################

foreach $chunk (1..$total_chunks){
    print OUT "<TR><TD><TD><CENTER>$chunk\n";
    print OUT "<TD><CENTER><IMG SRC=tick.png ALT=PASS>\n" if $sum{'multisample'}{$chunk}{'raw_ug_vcf_chunk'} eq "PASS";
    print OUT "<TD><CENTER><IMG SRC=error.png ALT=FAIL>\n" if $sum{'multisample'}{$chunk}{'raw_ug_vcf_chunk'} eq "FAIL";
    
    if ($chunk == 1){
		print OUT "<TD><CENTER><IMG SRC=tick.png ALT=PASS>\n" if $sum{'multisample'}{'raw_ug_vcf'} eq "PASS";
    	print OUT "<TD><CENTER><IMG SRC=error.png ALT=FAIL>\n" if $sum{'multisample'}{'raw_ug_vcf'} eq "FAIL";

		unless ($type eq "TARGETED"){
	
			#deploy variant recalibration plots
        	$pdf_snp = "$project_dir_results/$date/multisample/unifiedgenotyper/recalibration/$project"."_multisample.SNP.plots.R.pdf";
        	$pdf_ind = "$project_dir_results/$date/multisample/unifiedgenotyper/recalibration/$project"."_multisample.INDEL.plots.R.pdf";
	    	$pdf_deployment_ug = "$summary_deployment/ug";
	    	system("ssh $deployment_server mkdir -p -m 775 $pdf_deployment_ug > /dev/null 2>&1");
        	system("scp -r $pdf_snp $deployment_server:$pdf_deployment_ug/UG.SNP.plots.R.pdf > /dev/null 2>&1");
        	system("scp -r $pdf_ind $deployment_server:$pdf_deployment_ug/UG.INDEL.plots.R.pdf > /dev/null 2>&1");

			print OUT "<TD><CENTER><A HREF = http://$deployment_server/$summary_link/ug/><IMG BORDER='5' SRC=tick.png 	ALT=PASS></A>\n" if $sum{'multisample'}{'calibrated_ug_vcf'} eq "PASS";
			print OUT "<TD><CENTER><A HREF = http://$deployment_server/$summary_link/ug/><IMG BORDER='5' SRC=error.png ALT=FAIL></A>\n" if $sum{'multisample'}{'calibrated_ug_vcf'} eq "FAIL";
	
		}

	    $eval_report = "$project_dir_results/$date/multisample/unifiedgenotyper/recalibration/$project"."_multisample.varianteval.report";

		if (-s $eval_report){

		    print OUT "<TD><CENTER>$report{'unifiedgenotyper'}{'raw'}{'total'}";
		    print OUT "<TD><CENTER>$report{'unifiedgenotyper'}{'raw'}{'snp'}";
	        print OUT "<TD><CENTER>$report{'unifiedgenotyper'}{'raw'}{'TiTvall'}";

		    print OUT "<TD><CENTER>$report{'unifiedgenotyper'}{'called'}{'total'}";
		    print OUT "<TD><CENTER>$report{'unifiedgenotyper'}{'called'}{'snp'}";
	        print OUT "<TD><CENTER>$report{'unifiedgenotyper'}{'called'}{'TiTvall'}";
		    print OUT "<TD><CENTER>$report{'unifiedgenotyper'}{'called'}{'TiTvnovel'}";
		    print OUT "<TD><CENTER>$report{'unifiedgenotyper'}{'called'}{'TiTvknown'}";
	        print OUT "<TD><CENTER>$report{'unifiedgenotyper'}{'called'}{'SNP'}";
		    print OUT "<TD><CENTER>$report{'unifiedgenotyper'}{'called'}{'insertion'}";
		    print OUT "<TD><CENTER>$report{'unifiedgenotyper'}{'called'}{'deletion'}";

	        print OUT "<TD><CENTER>$report{'unifiedgenotyper'}{'filtered'}{'total'}";
		    print OUT "<TD><CENTER>$report{'unifiedgenotyper'}{'filtered'}{'snp'}";
		    print OUT "<TD><CENTER>$report{'unifiedgenotyper'}{'filtered'}{'TiTvall'}";

		} #end of if (-s $eval_report)

    } #end of if ($chunk == 1)

} #end of foreach $chunk (1..$total_chunks)

print OUT "</TABLE>\n";

#deploy VCF files

if ($type eq "TARGETED"){
    $vcf = "$project_dir_results/$date/multisample/unifiedgenotyper/$project"."_multisample.raw.vcf.gz";
    $tbi = "$project_dir_results/$date/multisample/unifiedgenotyper/$project"."_multisample.raw.vcf.gz.tbi";
}else{
    $vcf = "$project_dir_results/$date/multisample/unifiedgenotyper/$project"."_multisample.recalibrated.PASS.vcf.gz";
    $tbi = "$project_dir_results/$date/multisample/unifiedgenotyper/$project"."_multisample.recalibrated.PASS.vcf.gz.tbi";
}

system("scp -r $vcf $deployment_server:$data_deployment/$enc_dir/ug.$project.$date.vcf.gz > /dev/null 2>&1");
system("scp -r $tbi $deployment_server:$data_deployment/$enc_dir/ug.$project.$date.vcf.gz.tbi > /dev/null 2>&1");

#build and append UCSC Genome Browser link to VCF file 
$customTrackURL = "http://$deployment_server/$data_link/$enc_dir/ug.$project.$date.vcf.gz";
$directoryURL = "http://$deployment_server/$data_link/$enc_dir/";
$customTrack = "track type=vcfTabix name=$project visibility=full bigDataUrl=$customTrackURL";
$url = "http://$deployment_server/ucsc/cgi-bin/hgTracks?org=human&db=hg19&position=chr21:33,031,597-33,041,570&hgct_customText=$customTrack";

if ($sum{'multisample'}{'raw_ug_vcf'} eq "PASS"){
    print OUT "<P><FONT SIZE = '+1'>Results can be viewed as a <A HREF = '$url'>custom track</A> in the UCSC Genome Browser</FONT><BR>";
    print OUT "<P><FONT SIZE = '+1'>Download vcf file <A HREF = '$directoryURL'>here</A></FONT><BR>";
}
}


#print HaplotypeCaller statistics
#################################

if ($variant_caller =~ /H/){

	#horizontal line
	print OUT "<HR>";

	#table header
	#############

	print OUT "<TABLE CELLPADDING=5><TR>";
	print OUT "<TR><TH ROWSPAN='2'><CENTER>MULTISAMPLE<BR>HaplotypeCaller";
	print OUT "<TH ROWSPAN='2'><CENTER>Chunk";
	print OUT "<TH ROWSPAN='2'><CENTER>Raw vcf<BR>chunk";
	print OUT "<TH ROWSPAN='2'><CENTER>Raw vcf<BR>full";
	print OUT "<TH ROWSPAN='2'><CENTER>Calibrated vcf" unless $type eq "TARGETED";

	print OUT "<TH COLSPAN='3'><CENTER>Raw Variants";
	print OUT "<TH COLSPAN='8'><CENTER>Called Variants";
	print OUT "<TH COLSPAN='3'><CENTER>Filtered Variants";

	print OUT "<TR><TH><CENTER>Total";
	print OUT "<TH><CENTER>Known,%";
	print OUT "<TH><CENTER>TiTv<BR>total";

	print OUT "<TH><CENTER>Total";
	print OUT "<TH><CENTER>Known,%";
	print OUT "<TH><CENTER>TiTv<BR>total";
	print OUT "<TH><CENTER>TiTv<BR>known";
	print OUT "<TH><CENTER>TiTv<BR>novel";
	print OUT "<TH><CENTER>SNPs<BR>total";
	print OUT "<TH><CENTER>Insertions<BR>total";
	print OUT "<TH><CENTER>Deletions<BR>total";

	print OUT "<TH><CENTER>Total";
	print OUT "<TH><CENTER>Known,%";
	print OUT "<TH><CENTER>TiTv<BR>total";

	#table rows
	#variant stats by chunk
	#######################

	foreach $chunk (1..$total_chunks){

    	if ($type eq "TARGETED"){

			foreach $interval (1..$total_intervals){

				$fragment = "$chunk.interval_$interval";
				print OUT "<TR><TD><TD><CENTER>$interval\n";
				print OUT "<TD><CENTER><IMG SRC=tick.png ALT=PASS>\n" if $sum{'multisample'}{$fragment}{'raw_hc_vcf_chunk'} eq "PASS";
				print OUT "<TD><CENTER><IMG SRC=error.png ALT=FAIL>\n" if $sum{'multisample'}{$fragment}{'raw_hc_vcf_chunk'} eq "FAIL";

				if ($interval == 1){

					print OUT "<TD><CENTER><IMG SRC=tick.png ALT=PASS>\n" if $sum{'multisample'}{'raw_hc_vcf'} eq "PASS";
					print OUT "<TD><CENTER><IMG SRC=error.png ALT=FAIL>\n" if $sum{'multisample'}{'raw_hc_vcf'} eq "FAIL";

					$eval_report = "$project_dir_results/$date/multisample/unifiedgenotyper/recalibration/$project"."_multisample.varianteval.report";

			        if (-s $eval_report){
	       
			            print OUT "<TD><CENTER>$report{'haplotypecaller'}{'raw'}{'total'}";
			            print OUT "<TD><CENTER>$report{'haplotypecaller'}{'raw'}{'snp'}";
			            print OUT "<TD><CENTER>$report{'haplotypecaller'}{'raw'}{'TiTvall'}";
		
			            print OUT "<TD><CENTER>$report{'haplotypecaller'}{'called'}{'total'}";
			            print OUT "<TD><CENTER>$report{'haplotypecaller'}{'called'}{'snp'}";
			            print OUT "<TD><CENTER>$report{'haplotypecaller'}{'called'}{'TiTvall'}";
			            print OUT "<TD><CENTER>$report{'haplotypecaller'}{'called'}{'TiTvnovel'}";
			            print OUT "<TD><CENTER>$report{'haplotypecaller'}{'called'}{'TiTvknown'}";
			            print OUT "<TD><CENTER>$report{'haplotypecaller'}{'called'}{'SNP'}";
			            print OUT "<TD><CENTER>$report{'haplotypecaller'}{'called'}{'insertion'}";
			            print OUT "<TD><CENTER>$report{'haplotypecaller'}{'called'}{'deletion'}";

			            print OUT "<TD><CENTER>$report{'haplotypecaller'}{'filtered'}{'total'}";
			            print OUT "<TD><CENTER>$report{'haplotypecaller'}{'filtered'}{'snp'}";
			            print OUT "<TD><CENTER>$report{'haplotypecaller'}{'filtered'}{'TiTvall'}";

		            } # end of if (-s $eval_report){

		        } #end of if ($interval == 1){

    		} #end of foreach $interval (1..$total_intervals)

		}else{

	        print OUT "<TR><TD><TD><CENTER>$chunk\n";
	        print OUT "<TD><CENTER><IMG SRC=tick.png ALT=PASS>\n" if $sum{'multisample'}{$chunk}{'raw_hc_vcf_chunk'} eq "PASS";
	        print OUT "<TD><CENTER><IMG SRC=error.png ALT=FAIL>\n" if $sum{'multisample'}{$chunk}{'raw_hc_vcf_chunk'} eq "FAIL";

	        if ($chunk == 1){

		   	    print OUT "<TD><CENTER><IMG SRC=tick.png ALT=PASS>\n" if $sum{'multisample'}{'raw_hc_vcf'} eq "PASS";
			    print OUT "<TD><CENTER><IMG SRC=error.png ALT=FAIL>\n" if $sum{'multisample'}{'raw_hc_vcf'} eq "FAIL";

	            $pdf_snp = "$project_dir_results/$date/multisample/haplotypecaller/recalibration/$project"."_multisample.SNP.plots.R.pdf";
	            $pdf_ind = "$project_dir_results/$date/multisample/haplotypecaller/recalibration/$project"."_multisample.INDEL.plots.R.pdf";
			    $pdf_deployment_hc = "$summary_deployment/hc";

			    system("ssh $deployment_server mkdir -p -m 775 $pdf_deployment_hc > /dev/null 2>&1");
	            system("scp -r $pdf_snp $deployment_server:$pdf_deployment_hc/HC.SNP.plots.R.pdf > /dev/null 2>&1");
	            system("scp -r $pdf_ind $deployment_server:$pdf_deployment_hc/HC.INDEL.plots.R.pdf > /dev/null 2>&1");

			    print OUT "<TD><CENTER><A HREF = http://$deployment_server/$summary_link/hc/><IMG BORDER='5' SRC=tick.png ALT=PASS></A>\n" if $sum{'multisample'}{'calibrated_hc_vcf'} eq "PASS";
			    print OUT "<TD><CENTER><A HREF = http://$deployment_server/$summary_link/hc/><IMG BORDER='5' SRC=error.png ALT=FAIL></A>\n" if $sum{'multisample'}{'calibrated_hc_vcf'} eq "FAIL";

	            $eval_report = "$project_dir_results/$date/multisample/haplotypecaller/recalibration/$project"."_multisample.varianteval.report";

			    if (-s $eval_report){
     
			        print OUT "<TD><CENTER>$report{'haplotypecaller'}{'raw'}{'total'}";
			        print OUT "<TD><CENTER>$report{'haplotypecaller'}{'raw'}{'snp'}";
			        print OUT "<TD><CENTER>$report{'haplotypecaller'}{'raw'}{'TiTvall'}";

			        print OUT "<TD><CENTER>$report{'haplotypecaller'}{'called'}{'total'}";
			        print OUT "<TD><CENTER>$report{'haplotypecaller'}{'called'}{'snp'}";
			        print OUT "<TD><CENTER>$report{'haplotypecaller'}{'called'}{'TiTvall'}";
			        print OUT "<TD><CENTER>$report{'haplotypecaller'}{'called'}{'TiTvnovel'}";
			        print OUT "<TD><CENTER>$report{'haplotypecaller'}{'called'}{'TiTvknown'}";
			        print OUT "<TD><CENTER>$report{'haplotypecaller'}{'called'}{'SNP'}";
			        print OUT "<TD><CENTER>$report{'haplotypecaller'}{'called'}{'insertion'}";
			        print OUT "<TD><CENTER>$report{'haplotypecaller'}{'called'}{'deletion'}";

			        print OUT "<TD><CENTER>$report{'haplotypecaller'}{'filtered'}{'total'}";
			        print OUT "<TD><CENTER>$report{'haplotypecaller'}{'filtered'}{'snp'}";
			        print OUT "<TD><CENTER>$report{'haplotypecaller'}{'filtered'}{'TiTvall'}";

	            }

	        }

	    }

	}
	print OUT "</TABLE>\n";


	#deploy VCF files

	if ($type eq "TARGETED"){

	    $vcf = "$project_dir_results/$date/multisample/haplotypecaller/$project"."_multisample.raw.vcf.gz";
	    $tbi = "$project_dir_results/$date/multisample/haplotypecaller/$project"."_multisample.raw.vcf.gz.tbi";

	}else{

	    $vcf = "$project_dir_results/$date/multisample/haplotypecaller/$project"."_multisample.recalibrated.PASS.vcf.gz";
	    $tbi = "$project_dir_results/$date/multisample/haplotypecaller/$project"."_multisample.recalibrated.PASS.vcf.gz.tbi";

	}

	system("scp -r $vcf $deployment_server:$data_deployment/$enc_dir/hc.$project.$date.vcf.gz > /dev/null 2>&1");
	system("scp -r $tbi $deployment_server:$data_deployment/$enc_dir/hc.$project.$date.vcf.gz.tbi > /dev/null 2>&1");

	#build and append UCSC Genome Browser link to VCF file

	$customTrackURL = "http://$deployment_server/$data_link/$enc_dir/hc.$project.$date.vcf.gz";
	$directoryURL = "http://$deployment_server/$data_link/$enc_dir/";
	$customTrack = "track type=vcfTabix name=$project visibility=full bigDataUrl=$customTrackURL";
	$url = "http://$deployment_server/ucsc/cgi-bin/hgTracks?org=human&db=hg19&position=chr21:33,031,597-33,041,570&hgct_customText=$customTrack";

	if ($sum{'multisample'}{'raw_hc_vcf'} eq "PASS"){

	    print OUT "<P><FONT SIZE = '+1'>Results can be viewed as a <A HREF = '$url'>custom track</A> in the UCSC Genome Browser</FONT><BR>";
	    print OUT "<P><FONT SIZE = '+1'>Download vcf file <A HREF = '$directoryURL'>here</A></FONT><BR>";

	}

}

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

print OUT "<P><FONT SIZE = '+1'><A HREF = '$sample_summary_php_url'>Coverage Summary</A> ".
							   "<A HREF = '$sample_summary_url'>[TSV]</A> ".
							   "<A HREF = '$sample_summary_url.xlsx'>[XLS]</A>".
			  "</FONT><BR>";
			  
print OUT "<P><FONT SIZE = '+1'><A HREF = '$sample_interval_summary_php_url'>Interval Coverage Summary</A> ".
							   "<A HREF = '$sample_interval_summary_url'>[TSV]</A> ".
							   "<A HREF = '$sample_interval_summary_url.xlsx'>[XLS]</A> ".
							   "<A HREF = '$sample_interval_summary_url.png'>[Plot]</A>".
			"</FONT><BR>";
			
print OUT "<P><FONT SIZE = '+1'><A HREF = '$sample_cumulative_coverage_proportions_php_url'>Cumulative Coverage Proportions</A> ".
							   "<A HREF = '$sample_cumulative_coverage_proportions_url'>[TSV]</A> ".
							   "<A HREF = '$sample_cumulative_coverage_proportions_url.xlsx'>[XLS]</A> ".
							   "<A HREF = '$sample_cumulative_coverage_proportions_url.png'>[Plot]</A> ".
			"</FONT><BR>";

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
system("ssh $deployment_server chmod 775 $data_deployment/$enc_dir/* > /dev/null 2>&1");




