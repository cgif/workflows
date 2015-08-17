#!/usr/bin/perl -w

use strict;
use warnings;

my $project_dir_analysis = "projectDirAnalysis";
my $project_dir_results = "projectDirResults";
my $date = "Today";
my $project = "Project";
my $deployment_server = "deploymentServer";
my $summary_results = "summaryResults";
my $summary_deployment = "summaryDeployment";
my $sample_list = "sampleList";
my $mark_duplicates = "markDuplicates";
my $metric_level="metricLevel";
my $collect_metric="collectMetric";
my $url = "http://$deployment_server/$1" if $summary_deployment =~ /html\/(.*)/;
my $multisample = "multisample_no";


#create deployment directory for pdf files
system("ssh $deployment_server mkdir $summary_deployment/pdf > /dev/null 2>&1");

#collect list of all samples, libraries and read groups
my %data = ();
my $head = -1;
open (LIST, "$sample_list");
while (<LIST>){
    $head++;
    next unless $head;
    /^(\S+)\t(\S+)\t(\S+)\t/;
    $data{$2}{$3}{$1}++;
}

#collect data
my ($log, $dupmark, $merge_sample_bam, $flagstat);
my %sum = ();

#for each sample
foreach my $sample (keys %data){
    #check if reads were clipped and the number of reads before and after clipping is the same at read group level
    foreach my $lib (keys %{$data{$sample}}){
	foreach my $rg (keys %{$data{$sample}{$lib}}){
	    $log = "$project_dir_analysis/$date/$sample/ClipReads.$sample.$rg.log";
	    if (-s $log){
		$sum{$sample}{$lib}{$rg}{'clip_reads'} = "PASS";
		$sum{$sample}{$lib}{$rg}{'clip_reads'} = "FAIL" if `grep 'Number of reads before and after clipping is not the same' $log`;
	    }
	}
    }

    #check if duplicates are marked at library level
    $log = "$project_dir_analysis/$date/$sample/samtoolsMergeAndTag.$sample.log";
    if (-s $log){
        foreach my $library (keys %{$data{$sample}}){
            if ("$mark_duplicates" eq "TRUE" ){
	        $dupmark = "$project_dir_results/$date/$sample/$sample"."_$library.dupmark.stats";
	        if (-s $dupmark){
                    $sum{$sample}{$library}{'0'}{'remove_dupl'} = "PASS";
	        }else{
                    $sum{$sample}{$library}{'0'}{'remove_dupl'} = "FAIL";
	        }
            }
        }

	#check if bam files are merged at sample level and the number of reads is the same before and after merging
	$merge_sample_bam = "$project_dir_results/$date/$sample/$sample.bam";
	if (-s $merge_sample_bam){
            $sum{$sample}{'0'}{'0'}{'merge_sample'} = "PASS";
	}else{
            $sum{$sample}{'0'}{'0'}{'merge_sample'} = "FAIL";
	}  

        $sum{$sample}{'0'}{'0'}{'merge_sample'} = "FAIL" if `grep 'Number of reads before and after merging is not the same' $log`;

	#collect alignment statistics from the flagstat file
	$flagstat="$project_dir_results/$date/$sample/$sample".".bam.flagstat";
	if (-s $flagstat){

	    open (FLAGSTAT, "$flagstat");
	    while(<FLAGSTAT>){
	        $sum{$sample}{'0'}{'0'}{'total'} = $1 if /^(\d+) \+ \d+ in total \(/;
	        $sum{$sample}{'0'}{'0'}{'duplicates'} = $1 if /^(\d+) \+ \d+ duplicates$/;
	        $sum{$sample}{'0'}{'0'}{'mapped'} = $1 if /^(\d+) \+ \d+ mapped \(/;
	        $sum{$sample}{'0'}{'0'}{'paired'} = $1 if /^(\d+) \+ \d+ properly paired \(/;
            }

	    $sum{$sample}{'0'}{'0'}{'duplicates_pct'} = ($sum{$sample}{'0'}{'0'}{'duplicates'}/$sum{$sample}{'0'}{'0'}{'total'})*100;
	    $sum{$sample}{'0'}{'0'}{'mapped_pct'} = ($sum{$sample}{'0'}{'0'}{'mapped'}/$sum{$sample}{'0'}{'0'}{'total'})*100;
	    $sum{$sample}{'0'}{'0'}{'paired_pct'} = ($sum{$sample}{'0'}{'0'}{'paired'}/$sum{$sample}{'0'}{'0'}{'total'})*100;

	}else{

	    $sum{$sample}{'0'}{'0'}{'total'} = "NA";
	    $sum{$sample}{'0'}{'0'}{'duplicates'} = "NA";
	    $sum{$sample}{'0'}{'0'}{'mapped'} = "NA";
	    $sum{$sample}{'0'}{'0'}{'paired'} = "NA";

	    $sum{$sample}{'0'}{'0'}{'duplicates_pct'} = "NA";
	    $sum{$sample}{'0'}{'0'}{'mapped_pct'} = "NA";
	    $sum{$sample}{'0'}{'0'}{'paired_pct'} = "NA";

	}

    }

}

#print summary
open (OUT, ">$summary_results/index.html");
print OUT "<HTML>";
print OUT "<HEAD><META HTTP-EQUIV='refresh' CONTENT='60'></HEAD>";
print OUT "<BODY><TABLE CELLPADDING=5><TR>";
print OUT "<TH><CENTER>Sample<TH><CENTER>Library<TH>Read group";
print OUT "<TH><CENTER>Clip reads";
print OUT "<TH><CENTER>Mark Duplicates" if ("$mark_duplicates" eq "TRUE");
print OUT "<TH><CENTER>Merge Bam";
print OUT "<TH>Total<BR>reads<TH>Duplicated<BR>reads<TH>Mapped<BR>reads<TH>Reads in<BR>concordant pairs";
print OUT "<TH>Sample<BR>metrics" if $metric_level =~ /S/;
print OUT "<TH>Library<BR>metrics" if $metric_level =~ /L/;
print OUT "<TH>Read group<BR>metrics" if $metric_level =~ /RG/;

my $f1 = 0;
my $f2 = 0;
my ($source, $destination);

#for each sample, library, read group
foreach my $sample (sort {$a cmp $b} keys %data){
    print OUT "<TR><TD>$sample";
    $log = "$project_dir_analysis/$date/$sample/samtoolsMergeAndTag.$sample.log";
    next unless (-s $log);
    $f1 = 1; #f1 is flagging for the first read group for each sample

    foreach my $library (sort {$a cmp $b} keys %{$data{$sample}}){
	print OUT "<TR><TD><TD>" unless $f1;
	print OUT "<TD>$library";
	$f2 = 1; #f2 is flagging for the first read group for each library

	foreach my $read (sort {$a cmp $b} keys %{$data{$sample}{$library}}){
	    print OUT "<TR><TD><TD>" unless $f2;
	    print OUT "<TD>$read";

	    #print results of reads clipping
	    if ("$sum{$sample}{$library}{$read}{'clip_reads'}" eq "PASS"){
		print OUT "<TD><CENTER><IMG SRC=tick.png ALT=PASS>";
	    }elsif ("$sum{$sample}{$library}{$read}{'clip_reads'}" eq "FAIL"){
		print OUT "<TD><CENTER><IMG SRC=error.png ALT=FAIL>"
	    }else{
		print OUT "<TD>";
	    }

	    #print results of duplicates marking
	    if ("$mark_duplicates" eq "TRUE"){
		if ($f2){
		    print OUT "<TD><CENTER><IMG SRC=tick.png ALT=PASS>" if $sum{$sample}{$library}{'0'}{'remove_dupl'} eq "PASS";
		    print OUT "<TD><CENTER><IMG SRC=error.png ALT=FAIL>" if $sum{$sample}{$library}{'0'}{'remove_dupl'} eq "FAIL";
		}else{
		    print OUT "<TD>";
		}
	    }

	    #print results of merging bam files and alignment statistics
	    if ($f1){
		print OUT "<TD><CENTER><IMG SRC=tick.png ALT=PASS>" if $sum{$sample}{'0'}{'0'}{'merge_sample'} eq "PASS";
		print OUT "<TD><CENTER><IMG SRC=error.png ALT=FAIL>" if $sum{$sample}{'0'}{'0'}{'merge_sample'} eq "FAIL";
		print OUT "<TD><CENTER>$sum{$sample}{'0'}{'0'}{'total'}";
		printf OUT ("<TD><CENTER>$sum{$sample}{'0'}{'0'}{'duplicates'} (%.1f%%)", $sum{$sample}{'0'}{'0'}{'duplicates_pct'});
		printf OUT ("<TD><CENTER>$sum{$sample}{'0'}{'0'}{'mapped'} (%.1f%%)", $sum{$sample}{'0'}{'0'}{'mapped_pct'});
		printf OUT ("<TD><CENTER>$sum{$sample}{'0'}{'0'}{'paired'} (%.1f%%)", $sum{$sample}{'0'}{'0'}{'paired_pct'});
	    }else{
		print OUT "<TD><TD><TD><TD><TD>";
	    }

	    #add links to plots
	    if ($metric_level =~ /S/){
		print OUT "<TD>";
		if ($f1){
		    foreach my $ext (qw (gcBias insert_size_histogram quality_by_cycle quality_distribution)){
			$source = "$project_dir_results/$date/$sample/$sample".".$ext".".pdf";
			$destination = "$sample".".$ext".".pdf";
			system("scp -r $source $deployment_server:$summary_deployment/pdf/$destination > /dev/null 2>&1");
		        print OUT "<A HREF = '$url/pdf/$destination'>$ext</A><BR>";
		    }
		}
	    }

	    if ($metric_level =~ /L/){
		print OUT "<TD>";
		if ($f2){
		    foreach my $ext (qw (gcBias insert_size_histogram)){
			$source = "$project_dir_results/$date/$sample/$sample"."_$library".".$ext".".pdf";
			$destination = "$sample"."_$library".".$ext".".pdf";
			system("scp -r $source $deployment_server:$summary_deployment/pdf/$destination > /dev/null 2>&1");
			print OUT "<A HREF = '$url/pdf/$destination'>$ext</A><BR>";
		    }
		}
	    }

            if ($metric_level =~ /RG/){
		print OUT "<TD>";
		foreach my $ext (qw (quality_by_cycle quality_distribution)){
		    $source = "$project_dir_results/$date/$sample/$sample"."_$read".".$ext".".pdf";
		    $destination = "$sample"."_$read".".$ext".".pdf";
		    system("scp -r $source $deployment_server:$summary_deployment/pdf/$destination > /dev/null 2>&1");
		    print OUT "<A HREF = '$url/pdf/$destination'>$ext</A><BR>";
		}
	    }
	    $f1 = 0;
	    $f2 = 0;
	    print OUT "\n";
        }
    }
}

print OUT "</TABLE>";

#deploying metrics once mergeandtag run is completed
if ("$multisample" eq "multisample_yes"){

system("ssh $deployment_server mkdir $summary_deployment/metrics > /dev/null 2>&1");
my $metrics_path = "$project_dir_results/$date/multisample";
my $html_path = "$project_dir_analysis/$date/multisample";

#create merged pdf with metrics at different levels, deploy them to the Web server and print links on the summary page
my ($metrics_name, $metrics_file, $counts_name, $counts_file, $read_groups, $cur_file, $libraries, $samples, $html_file, $lines, $png_file, $chart_name, $chart_file);
if ($metric_level =~ /RG/){
    print OUT "<HR><TABLE><TR><TD><FONT SIZE = '+1'>Read Group metrics";
    foreach my $category (qw(quality_by_cycle quality_distribution)){
	$metrics_name = "$project.$date.$category.read_group.pdf";
	$metrics_file = "$metrics_path/$category.read_group.pdf";
	$read_groups = "";
	foreach my $sample (sort {$a cmp $b} keys %data){
	    foreach my $library (sort {$a cmp $b} keys %{$data{$sample}}){
		foreach my $read (sort {$a cmp $b} keys %{$data{$sample}{$library}}){
		    $cur_file = "$project_dir_results/$date/$sample/$sample"."_"."$read.$category.pdf";
		    $read_groups .= "$cur_file " if (-s $cur_file);
		}
	    }
	}
	chop($read_groups);
	system("gs -dBATCH -dNOPAUSE -q -sDEVICE=pdfwrite -sOutputFile=$metrics_file $read_groups");
	system("scp -r $metrics_file $deployment_server:$summary_deployment/metrics/$metrics_name");
	print OUT "<TR><TD><TD><A HREF = '$url/metrics/$metrics_name'>$category</A>";
    }
    print OUT "</FONT></TABLE><BR>";
}

if ($metric_level =~ /L/){
    print OUT "<HR><TABLE><TR><TD><FONT SIZE = '+1'>Library metrics";
    foreach my $category (qw(gcBias insert_size_histogram)){
	$metrics_name = "$project.$date.$category.library.pdf";
	$metrics_file = "$metrics_path/$category.library.pdf";
	$libraries = "";
	foreach my $sample (sort {$a cmp $b} keys %data){
	    foreach my $library (sort {$a cmp $b} keys %{$data{$sample}}){
		$cur_file = "$project_dir_results/$date/$sample/$sample"."_"."$library.$category.pdf";
		$libraries .= "$cur_file " if (-s $cur_file);
	    }
	}
	chop($libraries);
	system("gs -dBATCH -dNOPAUSE -q -sDEVICE=pdfwrite -sOutputFile=$metrics_file $libraries");
	system("scp -r $metrics_file $deployment_server:$summary_deployment/metrics/$metrics_name");
	print OUT "<TR><TD><TD><A HREF = '$url/metrics/$metrics_name'>$category</A>";
    }
    print OUT "</FONT></TABLE><BR>";
}

if ($metric_level =~ /S/){
    print OUT "<HR><TABLE><TR><TD><FONT SIZE = '+1'>Sample metrics";
    foreach my $category (qw(quality_by_cycle quality_distribution gcBias insert_size_histogram)){
	$metrics_name = "$project.$date.$category.sample.pdf";
	$metrics_file = "$metrics_path/$category.sample.pdf";
	$samples = "";
	foreach my $sample (sort {$a cmp $b} keys %data){
	    $cur_file = "$project_dir_results/$date/$sample/$sample.$category.pdf";
	    $samples .= "$cur_file " if (-s $cur_file);
	}
	chop($samples);
	system("gs -dBATCH -dNOPAUSE -q -sDEVICE=pdfwrite -sOutputFile=$metrics_file $samples");
	system("scp -r $metrics_file $deployment_server:$summary_deployment/metrics/$metrics_name");
	print OUT "<TR><TD><TD><A HREF = '$url/metrics/$metrics_name'>$category</A>";
    }
    print OUT "</FONT></TABLE><BR>";
}

if ($metric_level =~ /S/){
    print OUT "<HR><TABLE><TR><TD><FONT SIZE = '+1'>Alignment summary metrics";
    foreach my $category (qw(FIRST_OF_PAIR SECOND_OF_PAIR PAIR UNPAIRED)){
	$metrics_name = "$project.$date.alignment_summary_metrics.$category";
	$metrics_file = "$metrics_path/$metrics_name";
	$html_file = "$html_path/$metrics_name.php";
        $lines = `wc -l $metrics_file|cut -f 1 -d ' '`;
	chomp($lines);
	if ($lines > 1){
	    system("scp -r $metrics_file $deployment_server:$summary_deployment/metrics/$metrics_name");
	    system("scp -r $html_file $deployment_server:$summary_deployment/metrics/$metrics_name.php");
	    print OUT "<TR><TD><TD><A HREF = '$url/metrics/$metrics_name.php'>$category</A>";
	}
    }
    print OUT "</FONT></TABLE><BR>";
}

#deploy files with summary metrics for different types of analysis to the Web server and print links on the summary page
if ($collect_metric =~ /TP/){
    $metrics_name = "$project.$date.targetedPcrMetrics";
    $metrics_file = "$metrics_path/$metrics_name";
    $html_file = "$html_path/$metrics_name.php";
    if (-s $metrics_file){
        system("scp -r $metrics_file* $deployment_server:$summary_deployment/metrics");
	system("scp -r $html_file $deployment_server:$summary_deployment/metrics/$metrics_name.php");
        print OUT "<P><FONT SIZE = '+1'><A HREF = '$url/metrics/$metrics_name.php'>Targeted PCR metrics</A>".
                                        "<A HREF = '$url/metrics/$metrics_name'> [TSV]</A>".
                                        "<A HREF = '$url/metrics/$metrics_name.readsAlligned.png'> [Plot aligned reads]</A>".
                                        "<A HREF = '$url/metrics/$metrics_name.basesOnBait.png'> [Plot reads ON and OFF amplicon]</A>".
                                        "<A HREF = '$url/metrics/$metrics_name.meanTargetCoverage.png'> [Plot mean target coverage]</A>".
                                        "<A HREF = '$url/metrics/$metrics_name.foldEnrichment.png'> [Plot fold enrichment]</A>".
                                        "<A HREF = '$url/metrics/$metrics_name.cumulativeCoverage.png'> [Plot cumulative coverage]</A>".
					"</FONT><BR>";
    }
    $metrics_name = "$project.$date.perTargetCoverage";
    $metrics_file = "$metrics_path/$metrics_name";
    $html_file = "$html_path/$metrics_name.php";
    if (-s $metrics_file){
        system("scp -r $metrics_file $deployment_server:$summary_deployment/metrics/$metrics_name");
	system("scp -r $html_file $deployment_server:$summary_deployment/metrics/$metrics_name.php");
        print OUT "<P><FONT SIZE = '+1'><A HREF = '$url/metrics/$metrics_name.php'>Target coverage (Picard)</A>".
                                        "<A HREF = '$url/metrics/$metrics_name'> [TSV]</A>".
					"</FONT><BR>";
    }
    $metrics_name = "$project.$date.non_overlapping.perTargetCoverage";
    $metrics_file = "$metrics_path/$metrics_name";
    $html_file = "$html_path/$metrics_name.php";
    if (-s $metrics_file){
        system("scp -r $metrics_file $deployment_server:$summary_deployment/metrics/$metrics_name");
	system("scp -r $html_file $deployment_server:$summary_deployment/metrics/$metrics_name.php");
        print OUT "<P><FONT SIZE = '+1'><A HREF = '$url/metrics/$metrics_name.php'>Target coverage non-overlapping amplicon regions (Picard)</A>".
                                        "<A HREF = '$url/metrics/$metrics_name'> [TSV]</A>".
					"</FONT><BR>";
    }
    $metrics_name = "$project.$date.sample_interval_summary";
    $metrics_file = "$metrics_path/$metrics_name";
    $html_file = "$html_path/$metrics_name.php";
    $png_file = "$metrics_path/$metrics_name.png";
    if (-s $metrics_file){
        system("scp -r $metrics_file $deployment_server:$summary_deployment/metrics/$metrics_name");
	system("scp -r $html_file $deployment_server:$summary_deployment/metrics/$metrics_name.php");
	system("scp -r $png_file $deployment_server:$summary_deployment/metrics/$metrics_name.png");
        print OUT "<P><FONT SIZE = '+1'><A HREF = '$url/metrics/$metrics_name.php'>Target coverage (GATK)</A>".
                                        "<A HREF = '$url/metrics/$metrics_name'> [TSV]</A>".
                                        "<A HREF = '$url/metrics/$metrics_name.png'> [Plot]</A>".
					"</FONT><BR>";
    }
    $metrics_name = "$project.$date.non_overlapping.sample_interval_summary";
    $metrics_file = "$metrics_path/$metrics_name";
    $html_file = "$html_path/$metrics_name.php";
    $png_file = "$metrics_path/$metrics_name.png";
    if (-s $metrics_file){
        system("scp -r $metrics_file $deployment_server:$summary_deployment/metrics/$metrics_name");
	system("scp -r $html_file $deployment_server:$summary_deployment/metrics/$metrics_name.php");
	system("scp -r $png_file $deployment_server:$summary_deployment/metrics/$metrics_name.png");
        print OUT "<P><FONT SIZE = '+1'><A HREF = '$url/metrics/$metrics_name.php'>Target coverage non-overlapping amplicon regions (GATK)</A>".
                                        "<A HREF = '$url/metrics/$metrics_name'> [TSV]</A>".
                                        "<A HREF = '$url/metrics/$metrics_name.png'> [Plot]</A>".
					"</FONT><BR>";
    }
}

if ($collect_metric =~ /HS/){
    $metrics_name = "$project.$date.hybridMetrics";
    $metrics_file = "$metrics_path/$metrics_name";
    $html_file = "$html_path/$metrics_name.php";
    if (-s $metrics_file){
        system("scp -r $metrics_file* $deployment_server:$summary_deployment/metrics");
	system("scp -r $html_file $deployment_server:$summary_deployment/metrics/$metrics_name.php");
        print OUT "<P><FONT SIZE = '+1'><A HREF = '$url/metrics/$metrics_name.php'>Hybrid metrics</A>".
                                        "<A HREF = '$url/metrics/$metrics_name.readsAlligned.png'> [Plot aligned reads]</A>".
                                        "<A HREF = '$url/metrics/$metrics_name.basesOnBait.png'> [Plot reads ON and OFF bait]</A>".
                                        "<A HREF = '$url/metrics/$metrics_name.meanTargetCoverage.png'> [Plot mean target coverage]</A>".
                                        "<A HREF = '$url/metrics/$metrics_name.foldEnrichment.png'> [Plot fold enrichment]</A>".
                                        "<A HREF = '$url/metrics/$metrics_name.cumulativeCoverage.png'> [Plot cumulative coverage]</A>".
					"</FONT><BR>";
    }
    $metrics_name = "$project.$date.perTargetCoverage";
    $metrics_file = "$metrics_path/$metrics_name";
    $html_file = "$html_path/$metrics_name.php";
    if (-s $metrics_file){
        system("scp -r $metrics_file $deployment_server:$summary_deployment/metrics/$metrics_name");
	system("scp -r $html_file $deployment_server:$summary_deployment/metrics/$metrics_name.php");
        print OUT "<P><FONT SIZE = '+1'><A HREF = '$url/metrics/$metrics_name.php'>Target coverage</A></FONT><BR>";
    }
}

if ($collect_metric =~ /RS/){
    $metrics_name = "$project.$date.RnaSeqMetrics";
    $metrics_file = "$metrics_path/$metrics_name";
    $counts_name = "$project.$date.readCounts";
    $counts_file = "$metrics_path/$counts_name";
    $html_file = "$html_path/$metrics_name.php";
    if (-s $metrics_file){
        system("scp -r $metrics_file* $deployment_server:$summary_deployment/metrics");
        system("scp -r $counts_file* $deployment_server:$summary_deployment/metrics");
	system("scp -r $html_file $deployment_server:$summary_deployment/metrics/$metrics_name.php");
        print OUT "<P><FONT SIZE = '+1'><A HREF = '$url/metrics/$metrics_name.php'>RNA-Seq metrics</A>".
					"<A HREF = '$url/metrics/$metrics_name.basesAlignment.png'> [Plot genomic distribution of aligned bases]</A>".
					"<A HREF = '$url/metrics/$metrics_name.median5primeTo3primeBias.png'> [Plot 5' to 3' bias]</A><BR>".
					"<A HREF = '$url/metrics/$counts_name.png'> [Plot chromosomal distribution]</A><BR>".
					"<A HREF = '$url/metrics/$counts_name.normalized.png'> [Plot chromosomal distribution normalized by the total length of transcripts]</A><BR>".
					"<A HREF = '$url/metrics/$counts_name.normalizedNoMTGL.png'> [Plot chromosomal distribution normalized by the total length of transcripts, MT and GL contigs removed]</A><BR>".
					"</FONT><BR>";
    }
    $chart_name = "$project.$date.chartOutput.pdf";
    $chart_file = "$metrics_path/$chart_name";
    if (-s $chart_file){
        system("scp -r $chart_file $deployment_server:$summary_deployment/metrics/$chart_name");
        print OUT "<P><FONT SIZE = '+1'><A HREF = '$url/metrics/$chart_name'>RNA integrity chart</A></FONT><BR>";
    }
}
}

print OUT "</BODY></HTML>";

#deploy metrics and plots
system("scp -r $summary_results/index.html $deployment_server:$summary_deployment/index.html > /dev/null 2>&1");
system("ssh $deployment_server chmod 0664 $summary_deployment/* > /dev/null 2>&1");
system("ssh $deployment_server chmod 0664 $summary_deployment/pdf/* > /dev/null 2>&1");
system("ssh $deployment_server chmod -R 0775 $summary_deployment/pdf > /dev/null 2>&1");
system("ssh $deployment_server chmod -R 0775 $summary_deployment/metrics > /dev/null 2>&1");
