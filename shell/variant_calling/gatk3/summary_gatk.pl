#!/usr/bin/perl -w
use strict;

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
my $data_deployment = "dataDeployment";
my $sample_list = "sampleList";
my $ref_chunks = "refChunks";
my $ref_intervals = "refIntervals";
my $type = "#type";
my $enc_dir = "encryptedDir";
my $summary_link = $1 if  $summary_deployment =~ /\/www\/html\/(.*)/;
my $data_link = $1 if $data_deployment =~ /\/www\/html\/(.*)/;

my %sum = ();

open (CHUNKS, "$ref_chunks");
my $total_chunks = 0;
while (<CHUNKS>){
    chomp();
    next if /^\s*$/;   # skip blank lines
    my $chunk = $1 if /^\S+\t\S+\t\S+\t\S+\t(\S+)/;
    $total_chunks = $chunk if $chunk > $total_chunks;
}

my $total_intervals = 0;
if ($type eq "TARGETED" && $ref_intervals ne "null"){
    open (INTERVALS, "$ref_intervals");
    while (<INTERVALS>){
        chomp();
	next if /^\s*$/; # skip blank lines
        my $interval = $1 if /^\S+\t\S+\t\S+\t\S+\t(\S+)/;
        $total_intervals = $interval if $interval > $total_intervals;
    }
}

my @sample_all = ();

## TO DO - ensure that sample IDs are unique
open (LIST, "$sample_list");
while (<LIST>){
    push(@sample_all, $1) if /^(\S+)\s/;
}

my @sample = Uniq(@sample_all);

#collect data from gatk log files
foreach my $sample (@sample){
    foreach my $chunk (1..$total_chunks){
	my $realignment_log = "$project_dir_analysis/$date/$sample/run/gatk3_realign_recal.$sample".".chunk_"."$chunk".".log";
	if (-s $realignment_log){
   	    my $realignment_target = "$project_dir_analysis/$date/$sample/realignment/$sample".".chunk_"."$chunk".".RTC.intervals";
	    if (-s $realignment_target){
                $sum{$sample}{$chunk}{'realignment_target'} = "PASS";
	    }else{
                $sum{$sample}{$chunk}{'realignment_target'} = "FAIL";
	    }

   	    my $realigned_bam = "$project_dir_analysis/$date/$sample/realignment/$sample".".chunk_"."$chunk".".realigned.bam";
	    if (-s $realigned_bam){
                $sum{$sample}{$chunk}{'realigned_bam'} = "PASS";
	    }else{
                $sum{$sample}{$chunk}{'realigned_bam'} = "FAIL";
	    }

	    my $recalibration_report_chunk = "$project_dir_analysis/$date/$sample/recalibration/reports/pre/$sample".".chunk_"."$chunk".".realigned.recal_data.grp";
	    if (-s $recalibration_report_chunk){
                $sum{$sample}{$chunk}{'recalibration_report_chunk'} = "PASS";
	    }else{
                $sum{$sample}{$chunk}{'recalibration_report_chunk'} = "FAIL";
	    }
	}

        my $recalibration_hap_caller_log = "$project_dir_analysis/$date/$sample/run/gatk3_print_reads_hap_caller.$sample".".chunk_"."$chunk".".log";
        if (-s $recalibration_hap_caller_log){
	    my $recalibrated_bam = "$project_dir_analysis/$date/$sample/recalibration/$sample".".chunk_"."$chunk.realigned.recalibrated.bam";
	    if (-s $recalibrated_bam){
		$sum{$sample}{$chunk}{'recalibrated_bam'} = "PASS";
	    }else{
		$sum{$sample}{$chunk}{'recalibrated_bam'} = "FAIL";
	    }
	    my $haplotype_caller_grep = "";
	    $haplotype_caller_grep = `grep ERROR $recalibration_hap_caller_log`;
	    my $hc_gvcf_chunk = "$project_dir_analysis/$date/$sample/haplotypecaller/$sample".".chunk_"."$chunk.genomic.vcf";
	    if (-s $hc_gvcf_chunk && $haplotype_caller_grep eq ""){
		$sum{$sample}{$chunk}{'hc_gvcf_chunk'} = "PASS";
	    }else{
		$sum{$sample}{$chunk}{'hc_gvcf_chunk'} = "FAIL";
	    }
	}

	my $recalibration_merge_log = "$project_dir_analysis/$date/$sample/run/gatk3_merge_recalibration_reports.$sample".".log";
	if (-s $recalibration_merge_log){
	    my $recalibration_report = "$project_dir_results/$date/$sample/recalibration/reports/pre/$sample".".realigned.recal_data.grp";
	    if (-s $recalibration_report){
		$sum{$sample}{'recalibration_report'} = "PASS";
	    }else{
		$sum{$sample}{'recalibration_report'} = "FAIL";
	    }
	}

	my $merge_bam_gvcf_log = "$project_dir_analysis/$date/$sample/run/gatk3_merge_recal_chunk_bams_gvcfs.$sample".".log";
	if (-s $merge_bam_gvcf_log){
	    my $merged_gvcf = "$project_dir_results/$date/$sample/haplotypecaller/$sample".".genomic.vcf";
	    if (-s $merged_gvcf){
		$sum{$sample}{'merged_gvcf'} = "PASS";
	    }else{
		$sum{$sample}{'merged_gvcf'} = "FAIL";
	    }
	    my $merged_bam = "$project_dir_results/$date/$sample/recalibration/$sample".".bam";
	    if (-s $merged_bam){
		$sum{$sample}{'merged_bam'} = "PASS";
	    }else{
		$sum{$sample}{'merged_bam'} = "FAIL";
	    }
	}

    }
}

foreach my $chunk (1..$total_chunks){
    my $genotypeGVCFs_log = "$project_dir_analysis/$date/multisample/run/gatk3_genotypeGVCFs.chunk_"."$chunk".".log";
    if (-s $genotypeGVCFs_log){
	my $genotypeGVCFs_grep = "";
	$genotypeGVCFs_grep = `grep ERROR $genotypeGVCFs_log`;
	my $raw_vcf_chunk = "$project_dir_analysis/$date/multisample/genotypeGVCFs/$project"."_multisample.chunk_"."$chunk".".raw.vcf";
	if (-s $raw_vcf_chunk && $genotypeGVCFs_grep eq ""){
	    $sum{'multisample'}{$chunk}{'raw_vcf_chunk'} = "PASS";
	}else{
	    $sum{'multisample'}{$chunk}{'raw_vcf_chunk'} = "FAIL";
	}
    }
}

my $recalibrate_vcf_log =  "$project_dir_analysis/$date/multisample/run/gatk3_recal_vcf.$project"."_multisample.log";
if (-s $recalibrate_vcf_log){
    my $raw_vcf = "$project_dir_results/$date/multisample/genotypeGVCFs/$project"."_multisample.raw.vcf.gz";
    if (-s $raw_vcf){
        $sum{'multisample'}{'raw_vcf'} = "PASS"
    }else{
        $sum{'multisample'}{'raw_vcf'} = "FAIL"
    }
    my $calibrated_vcf = "$project_dir_results/$date/multisample/genotypeGVCFs/$project"."_multisample.recalibrated.PASS.vcf.gz";
    if (-s $calibrated_vcf){
        $sum{'multisample'}{'calibrated_vcf'} = "PASS"
    }else{
        $sum{'multisample'}{'calibrated_vcf'} = "FAIL"
    }
}

my %report = ();
my $eval_report = "$project_dir_results/$date/multisample/genotypeGVCFs/recalibration/$project"."_multisample.varianteval.report";
if (-s $eval_report){
    open (REPORT, "$eval_report");
    while (<REPORT>){
	if (/^CompOverlap\s+dbsnp\s+eval\s+raw\s+none\s+all\s+all\s+(\S+)\s+\S+\s+\S+\s+(\S+)\s/){ 
	    $report{'raw'}{'total'} = $1;
	    $report{'raw'}{'snp'} = $2;
	}
	if (/^CompOverlap\s+dbsnp\s+eval\s+called\s+none\s+all\s+all\s+(\S+)\s+\S+\s+\S+\s+(\S+)\s/){ 
	    $report{'called'}{'total'} = $1;
	    $report{'called'}{'snp'} = $2;
	}
	if (/^CompOverlap\s+dbsnp\s+eval\s+filtered\s+none\s+all\s+all\s+(\S+)\s+\S+\s+\S+\s+(\S+)\s/){ 
	    $report{'filtered'}{'total'} = $1;
	    $report{'filtered'}{'snp'} = $2;
	}
	if (/^CountVariants\s+dbsnp\s+eval\s+called\s+none\s+all\s+all\s+\S+\s+\S+\s+\S+\s+\S+\s+\S+\s+\S+\s+(\S+)\s+\S+\s+(\S+)\s+(\S+)\s/){
	    $report{'called'}{'SNP'} = $1;
	    $report{'called'}{'insertion'} = $2;
	    $report{'called'}{'deletion'} = $3;
	}
	if (/^TiTvVariantEvaluator\s+dbsnp\s+eval\s+raw\s+none\s+all\s+all\s+\S+\s+\S+\s+(\S+)\s/){
	    $report{'raw'}{'TiTvall'} = $1;
	}
	if (/^TiTvVariantEvaluator\s+dbsnp\s+eval\s+called\s+none\s+(\w+)\s+all\s+\S+\s+\S+\s+(\S+)\s/){
	    $report{'called'}{"TiTv$1"} = $2;
	}
	if (/^TiTvVariantEvaluator\s+dbsnp\s+eval\s+filtered\s+none\s+all\s+all\s+\S+\s+\S+\s+(\S+)\s/){
	    $report{'filtered'}{'TiTvall'} = $1;
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
print OUT "<TH><CENTER>Genomic VCF";
print OUT "<TH><CENTER>Merged<BR>genomic VCF";
print OUT "<TH><CENTER>Merged<BR>recalibrated bam";

#print extracted data into run log text file
open (LOG, ">$summary_results/run_log.txt");
print LOG "SAMPLE\tPROCESS\tCHUNK\tRESULT\n";


my @tags = qw(realignment_target realigned_bam recalibration_report_chunk recalibration_report recalibrated_bam hc_gvcf_chunk merged_gvcf merged_bam);
foreach my $sample (sort @sample){
    foreach my $chunk (1..$total_chunks){	
	if ($chunk == 1){
	    print OUT "<TR><TD>$sample<TD><CENTER>$chunk\n";
	}else{
	    print OUT "<TR><TD> <TD><CENTER>$chunk\n";
	}
        
	foreach my $tag (@tags){
	    if ($tag eq "recalibration_report"){
		if ($chunk == 1){
		    if (defined $sum{$sample}{$tag}){
			my $recalibration_merge_plot_log = "$project_dir_analysis/$date/$sample/run/gatk3_collect_summary_metrics.$sample".".log";
			if (-s $recalibration_merge_plot_log){
			    my $pdf = "$project_dir_results/$date/$sample/recalibration/plots/post/$sample.realigned.recalibrated.recal_plot-1.jpeg";
			    my $pdf_deployment = "$summary_deployment/$sample.jpeg";
			    if (-s $pdf){
				system("scp -r $pdf $deployment_server:$pdf_deployment > /dev/null 2>&1");
				print OUT "<TD><CENTER><A HREF = http://$deployment_server/$summary_link/$sample.jpeg><IMG BORDER='5' SRC=tick.png ALT=PASS></A>\n" if $sum{$sample}{$tag} eq "PASS";
				print OUT "<TD><CENTER><A HREF = http://$deployment_server/$summary_link/$sample.jpeg><IMG BORDER='5' SRC=error.png ALT=FAIL></A>\n" if $sum{$sample}{$tag} eq "FAIL";
				print OUT "<TD> \n" if $sum{$sample}{$tag} ne "FAIL" && $sum{$sample}{$tag} ne "PASS";
				print LOG $sample, "\t", $tag, "\tALL\t", $sum{$sample}{$tag}, "\n";
			    }else{
				print OUT "<TD><CENTER><IMG SRC=tick.png ALT=PASS></A>\n" if $sum{$sample}{$tag} eq "PASS";
				print OUT "<TD><CENTER><IMG SRC=error.png ALT=FAIL></A>\n" if $sum{$sample}{$tag} eq "FAIL";
				print OUT "<TD> \n" if $sum{$sample}{$tag} ne "FAIL" && $sum{$sample}{$tag} ne "PASS";
				print LOG $sample, "\t", $tag, "\tALL\t", $sum{$sample}{$tag}, "\n";
			    }
			}else{
			    print OUT "<TD><CENTER><IMG SRC=tick.png ALT=PASS></A>\n" if $sum{$sample}{$tag} eq "PASS";
			    print OUT "<TD><CENTER><IMG SRC=error.png ALT=FAIL></A>\n" if $sum{$sample}{$tag} eq "FAIL";
			    print OUT "<TD> \n" if $sum{$sample}{$tag} ne "FAIL" && $sum{$sample}{$tag} ne "PASS";
			    print LOG $sample, "\t", $tag, "\tALL\t", $sum{$sample}{$tag}, "\n";
			}
                    #end of if (-s $recalibration_merge_plot_log){
		    }else{
			print OUT "<TD> \n";
			print LOG $sample, "\t", $tag, "\tALL\tUNPROCESSED\n";
		    }
		}else{
		    print OUT "<TD> \n";
		}
		#end of if ($chunk == 1){
## add here if tag marged_bam and merged_gvcf
	    }elsif ($tag eq "merged_gvcf" || $tag eq "merged_bam"){
		if ($chunk == 1){
		    if (defined $sum{$sample}{$tag}){
			print OUT "<TD><CENTER><IMG SRC=tick.png ALT=PASS></A>\n" if $sum{$sample}{$tag} eq "PASS";
			print OUT "<TD><CENTER><IMG SRC=error.png ALT=FAIL></A>\n" if $sum{$sample}{$tag} eq "FAIL";
			print OUT "<TD> \n" if $sum{$sample}{$tag} ne "FAIL" && $sum{$sample}{$tag} ne "PASS";
			print LOG $sample, "\t", $tag, "\tALL\t", $sum{$sample}{$tag}, "\n";
		    }else{
			print OUT "<TD> \n";
			print LOG $sample, "\t", $tag, "\tALL\tUNPROCESSED\n";
		    }
		}else{
		    print OUT "<TD> \n";
		}
	    }else{
		if (defined $sum{$sample}{$chunk}{$tag}){
		    print OUT "<TD><CENTER><IMG SRC=tick.png ALT=PASS>\n" if $sum{$sample}{$chunk}{$tag} eq "PASS";
		    print OUT "<TD><CENTER><IMG SRC=error.png ALT=FAIL>\n" if $sum{$sample}{$chunk}{$tag} eq "FAIL";
		    print OUT "<TD> \n" if $sum{$sample}{$chunk}{$tag} ne "FAIL" && $sum{$sample}{$chunk}{$tag} ne "PASS";
		    print LOG $sample, "\t", $tag, "\t", $chunk, "\t", $sum{$sample}{$chunk}{$tag}, "\n";
		}else{
		    print OUT "<TD> \n";
		    print LOG $sample, "\t", $tag, "\t", $chunk, "\tUNPROCESSED\n";		   
		}
	    }
	    #end of if ($tag eq "recalibration_report"){
	}
	# end of foreach my $tag (@tags){
    }
    #end of foreach my $chunk (1..$total_chunks){
}
#end of foreach my $sample (@sample){

print OUT "</TABLE>\n";

#print GenotypeGVCFs statistics
##################################

#horizontal line
print OUT "<HR>";

#table header
#############
print OUT "<TABLE CELLPADDING=5>";
print OUT "<TR><TH ROWSPAN='2'><CENTER>MULTISAMPLE<BR>GenotypeGVCFs";
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

foreach my $chunk (1..$total_chunks){
    print OUT "<TR><TD><TD><CENTER>$chunk\n";
    if (defined $sum{'multisample'}{$chunk}{'raw_vcf_chunk'}){
	print OUT "<TD><CENTER><IMG SRC=tick.png ALT=PASS>\n" if $sum{'multisample'}{$chunk}{'raw_vcf_chunk'} eq "PASS";
	print OUT "<TD><CENTER><IMG SRC=error.png ALT=FAIL>\n" if $sum{'multisample'}{$chunk}{'raw_vcf_chunk'} eq "FAIL";
	print LOG "multisample\traw_vcf_chunk\t", $chunk, "\t", $sum{'multisample'}{$chunk}{'raw_vcf_chunk'}, "\n";
    }else{
	print LOG "multisample\traw_vcf_chunk\t", $chunk, "\tUNPROCESSED\n";
    }
    
    if ($chunk == 1){
	if (defined $sum{'multisample'}{'raw_vcf'}){
	    print OUT "<TD><CENTER><IMG SRC=tick.png ALT=PASS>\n" if $sum{'multisample'}{'raw_vcf'} eq "PASS";
	    print OUT "<TD><CENTER><IMG SRC=error.png ALT=FAIL>\n" if $sum{'multisample'}{'raw_vcf'} eq "FAIL";
	    print LOG "multisample\traw_vcf_full\tALL\t", $sum{'multisample'}{'raw_vcf'}, "\n";
	}else{
	    print LOG "multisample\traw_vcf_full\tALL\tUNPROCESSED\n";
	}

	unless ($type eq "TARGETED"){
	
        #deploy variant recalibration plots
	    my $pdf_snp = "$project_dir_results/$date/multisample/genotypeGVCFs/recalibration/$project"."_multisample.SNP.plots.R.pdf";
	    my $pdf_ind = "$project_dir_results/$date/multisample/genotypeGVCFs/recalibration/$project"."_multisample.INDEL.plots.R.pdf";
	    my $pdf_deployment_vr = "$summary_deployment/variant_recalibration";
	    system("ssh $deployment_server mkdir -p -m 775 $pdf_deployment_vr > /dev/null 2>&1");
	    system("scp -r $pdf_snp $deployment_server:$pdf_deployment_vr/SNP.plots.R.pdf > /dev/null 2>&1");
	    system("scp -r $pdf_ind $deployment_server:$pdf_deployment_vr/INDEL.plots.R.pdf > /dev/null 2>&1");

	    if (defined $sum{'multisample'}{'calibrated_vcf'}){
		print OUT "<TD><CENTER><A HREF = http://$deployment_server/$summary_link/variant_recalibration/><IMG BORDER='5' SRC=tick.png 	ALT=PASS></A>\n" if $sum{'multisample'}{'calibrated_vcf'} eq "PASS";
		print OUT "<TD><CENTER><A HREF = http://$deployment_server/$summary_link/variant_recalibration/><IMG BORDER='5' SRC=error.png ALT=FAIL></A>\n" if $sum{'multisample'}{'calibrated_vcf'} eq "FAIL";
		print LOG "multisample\tcalibrated_vcf\tALL\t", $sum{'multisample'}{'calibrated_vcf'}, "\n";
	    }else{
		print LOG "multisample\tcalibrated_vcf_full\tALL\tUNPROCESSED\n";
	    }
	
	}

#	my $eval_report = "$project_dir_results/$date/multisample/genotypeGVCFs/recalibration/$project"."_multisample.varianteval.report"; #was defined earlier

	if (-s $eval_report){

	    print OUT "<TD><CENTER>$report{'raw'}{'total'}";
	    print OUT "<TD><CENTER>$report{'raw'}{'snp'}";
	    print OUT "<TD><CENTER>$report{'raw'}{'TiTvall'}";

	    print OUT "<TD><CENTER>$report{'called'}{'total'}";
	    print OUT "<TD><CENTER>$report{'called'}{'snp'}";
	    print OUT "<TD><CENTER>$report{'called'}{'TiTvall'}";
	    print OUT "<TD><CENTER>$report{'called'}{'TiTvnovel'}";
	    print OUT "<TD><CENTER>$report{'called'}{'TiTvknown'}";
	    print OUT "<TD><CENTER>$report{'called'}{'SNP'}";
	    print OUT "<TD><CENTER>$report{'called'}{'insertion'}";
	    print OUT "<TD><CENTER>$report{'called'}{'deletion'}";

	    print OUT "<TD><CENTER>$report{'filtered'}{'total'}";
	    print OUT "<TD><CENTER>$report{'filtered'}{'snp'}";
	    print OUT "<TD><CENTER>$report{'filtered'}{'TiTvall'}";

	} #end of if (-s $eval_report)

    } #end of if ($chunk == 1)

} #end of foreach $chunk (1..$total_chunks)

print OUT "</TABLE>\n";

#deploy VCF files

my $vcf = "";
my $tbi = "";

if ($type eq "TARGETED"){
    $vcf = "$project_dir_results/$date/multisample/genotypeGVCFs/$project"."_multisample.raw.vcf.gz";
    $tbi = "$project_dir_results/$date/multisample/genotypeGVCFs/$project"."_multisample.raw.vcf.gz.tbi";
}else{
    $vcf = "$project_dir_results/$date/multisample/genotypeGVCFs/$project"."_multisample.recalibrated.PASS.vcf.gz";
    $tbi = "$project_dir_results/$date/multisample/genotypeGVCFs/$project"."_multisample.recalibrated.PASS.vcf.gz.tbi";
}

system("scp -r $vcf $deployment_server:$data_deployment/$enc_dir/$project.$date.vcf.gz > /dev/null 2>&1");
system("scp -r $tbi $deployment_server:$data_deployment/$enc_dir/$project.$date.vcf.gz.tbi > /dev/null 2>&1");

#build and append UCSC Genome Browser link to VCF file 
my $customTrackURL = "http://$deployment_server/$data_link/$enc_dir/$project.$date.vcf.gz";
my $directoryURL = "http://$deployment_server/$data_link/$enc_dir/";
my $customTrack = "track type=vcfTabix name=$project visibility=full bigDataUrl=$customTrackURL";
my $url = "http://$deployment_server/ucsc/cgi-bin/hgTracks?org=human&db=hg19&position=chr21:33,031,597-33,041,570&hgct_customText=$customTrack";

if ($sum{'multisample'}{'raw_vcf'} eq "PASS"){
    print OUT "<P><FONT SIZE = '+1'>Results can be viewed as a <A HREF = '$url'>custom track</A> in the UCSC Genome Browser</FONT><BR>";
    print OUT "<P><FONT SIZE = '+1'>Download vcf file <A HREF = '$directoryURL'>here</A></FONT><BR>";
}


#summary metrics
################

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

my $sample_summary_url = "http://$deployment_server/$summary_link/metrics/$project.$date.sample_summary";
my $sample_summary_php_url = "http://$deployment_server/$summary_link/metrics/sample_summary.php";

my $sample_interval_summary_url = "http://$deployment_server/$summary_link/metrics/$project.$date.sample_interval_summary";
my $sample_interval_summary_php_url = "http://$deployment_server/$summary_link/metrics/sample_interval_summary.php";

my $sample_cumulative_coverage_proportions_url = "http://$deployment_server/$summary_link/metrics/$project.$date.sample_cumulative_coverage_proportions";
my $sample_cumulative_coverage_proportions_php_url = "http://$deployment_server/$summary_link/metrics/sample_cumulative_coverage_proportions.php";

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

my $usage_file = "$project_dir_results/$date/multisample/usage.$date.txt";
my $usage_url = "http://$deployment_server/$summary_link/usage.txt";
if (-s $usage_file){
    system("scp -r $usage_file $deployment_server:$summary_deployment/usage.txt > /dev/null 2>&1");
    print OUT "<HR><P><FONT SIZE = '+1'>Usage of computational resources can be monitored <A HREF = '$usage_url'>here</A></FONT><BR>";
}

system("scp -r $summary_results/index.html $deployment_server:$summary_deployment/index.html > /dev/null 2>&1");
system("ssh $deployment_server chmod -R 775 $summary_deployment/* > /dev/null 2>&1");
system("ssh $deployment_server chmod 775 $data_deployment/$enc_dir/* > /dev/null 2>&1");

close (OUT);
close (LOG);

## Subroutine to get unique sample names

sub Uniq{
	my %temp_hash = map { $_, 0 } @_;
	return keys %temp_hash;
}



