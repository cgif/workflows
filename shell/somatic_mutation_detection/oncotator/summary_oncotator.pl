#!/usr/bin/perl -w

use strict;
use warnings;

#my $analysis_dir = "#analysisDir";
#my $results_dir = "#resultsDir";
#my $project = "#project";
#my $deployment_server = "#deploymentServer";
#my $summary_deployment = "#summaryDeployment";

my $analysis_dir = "/groupvol/cgi/analysis/johnson_glioma/oncotator/2014-11-27";
my $results_dir = "/groupvol/cgi/results/johnson_glioma/oncotator/2014-11-27";
my $project = "johnson_glioma";
my $deployment_server = "eliot.med.ic.ac.uk";
my $summary_deployment = "/www/html/report/project/johnson_glioma/oncotator/2014-11-27";

#collect data from vcf files
my %sum = (); 
my @samples = ();

#open directory with vcf files for individual samples
opendir (SAMPLES, "$analysis_dir") || print "Can't open $analysis_dir\n";
 
while (defined(my $sample_vcf = readdir(SAMPLES))){

	next unless $sample_vcf =~ /annotated/;

	my $sample = $1 if $sample_vcf =~ /(.*?)\..*/;
	push(@samples, $sample);

	$sample_vcf = "$analysis_dir/$sample_vcf";

	%sum = stats($sample, $sample_vcf, "sample", %sum);

}

#for combined vcf file...
my $project_vcf = "$results_dir/$project.annotated.vcf";
%sum = stats($project, $project_vcf, "project", %sum);

sub stats {

	my ($sample, $sample_vcf, $size, %sum) = @_;

	#calculate total number of variants
	$sum{$sample}{'Total'} = `cat $sample_vcf | grep -v '^#' | wc -l`;
	chomp($sum{$sample}{'Total'});

	open (VCF, "$sample_vcf") || print "Can't open $sample_vcf\n";
	while(my $var = <VCF>){

		chomp($var);
		next if $var =~ /^#/;

		my @data = (); my @info = ();
		@data = split(/\t/, $var);
		@info = split(/;/, $data[7]);

		#flag variants that come up in multiple samples
		my $f_mult = 0; 
		if ($size =~ /project/){
			foreach (@info){
				$f_mult++ if /set=.*-.*/;
			}
		}

		$sum{$sample}{'Total'}++;
		$sum{$sample}{'Total_mult'}++ if $f_mult > 2;

		#filter common variants
		my $f_common = 0;
		foreach (@info){

			if (/1000Genome_AF=(.*)/){
			    	$f_common++ if $1 > 0.01;
			}

			if (/ESP_MAF=(.*)/){
			    	my @AF = split(/\|/, $1);
			    	foreach my $AF (@AF){
					my @AF_all = split(/,/, $AF);
			    		$f_common++ if $AF_all[2] > 1;	
				}
			}

		}
		$sum{$sample}{'Common'}++ if $f_common;
		$sum{$sample}{'Common_mult'}++ if $f_common && $f_mult > 2;

		next if $f_common;

		#calculate number of SNPs and indels
		if ($data[7] =~ /VT=SNP/){
			$sum{$sample}{'SNVs'}++;
			$sum{$sample}{'SNVs_mult'}++ if $f_mult > 2;
		} else {
			$sum{$sample}{'Indels'}++;
			$sum{$sample}{'Indels_mult'}++ if $f_mult > 2;
		}

		my $f_lof_nonsyn = 0;
		my $f_cancer_mut = 0;
		my $f_repair_gene = 0;
		foreach (@info){

			#classify variance as LOF, nonsynonymous or non-coding
			if (/^variant_classification=/) {
				if (/Frame_Shift_Del|Frame_Shift_Ins|De_novo_Start_OutOfFrame|Start_Codon_SNP|Splice_Site|Nonsense_Mutation|Nonstop_Mutation/){
    					$sum{$sample}{'LOF'}++ ;
					$f_lof_nonsyn++;
    					$sum{$sample}{'LOF_mult'}++ if $f_mult > 2;
				} elsif (/In_Frame_Del|In_Frame_Ins|Missense_Mutation/) {
					$sum{$sample}{'Non-synonymous'}++;
					$f_lof_nonsyn++;
					$sum{$sample}{'Non-synonymous_mult'}++ if $f_mult > 2;
				} elsif (/Silent/) {

				} else {
					$sum{$sample}{'Non-coding'}++;
					$sum{$sample}{'Non-coding_mult'}++ if $f_mult > 2;
				}
			}

			#flag mutations seen in Cosmic database of somatic mutations in cancer samples
			#COSMIC can also return tissue type, primary site and total alterations in gene
    			if (/^COSMIC_n_overlapping_mutations=(\d+)/){
				$sum{$sample}{'COSMIC_mutations'}++ if $1 > 0;
				$f_cancer_mut++ if $1 > 0;
				$sum{$sample}{'COSMIC_mutations_mult'}++ if $1 > 0 && $f_mult > 2;
			}

			#flag mutations seen in ClinVar
			if (/ClinVar_SYM/) {
   				$sum{$sample}{'ClinVar_mutations'}++;
 				$f_cancer_mut++;
   				$sum{$sample}{'ClinVar_mutations_mult'}++ if $f_mult > 2; 
			}

			#flag mutations seen in CCLE
			if (/CCLE_By_GP_overlapping_mutations/) {
    				$sum{$sample}{'CCLE_mutations'}++;
  				$f_cancer_mut++;
    				$sum{$sample}{'CCLE_mutations_mult'}++ if $f_mult > 2; 
			}

			#flag mutations seen in the list of DNA repair genes
			if (/HumanDNARepairGenes_Role/) {
    				$sum{$sample}{'DNA_repair_genes'}++;
  				$f_repair_gene++;
    				$sum{$sample}{'DNA_repair_genes_mult'}++ if $f_mult > 2;
			}

		}

		if ($size =~ /project/){
			#look for overlaps between different categories
			$sum{'lof_nonsyn_cancer_mut'}++ if $f_lof_nonsyn && $f_cancer_mut;
			$sum{'lof_nonsyn'}++ if $f_lof_nonsyn;
			$sum{'cancer_mut'}++ if $f_cancer_mut;
			$sum{'lof_nonsyn_repair_genes'}++ if $f_lof_nonsyn && $f_repair_gene;
			$sum{'repair_genes'}++ if $f_repair_gene;
		}

	}

	return %sum;

}

#print results
open (XLS, ">$results_dir/$project.stats.xls");
print XLS "Sample\t";
foreach my $sample (sort @samples){
	print XLS "$sample\t";
}
print XLS "ALL SAMPLES(multiples)\n";

foreach my $key (qw(Total Common SNV Indel LOF Non-synonymous Non-coding COSMIC_mutations ClinVar_mutations CCLE_mutations DNA_repair_genes)) {
	print XLS "$key\t";
	foreach my $sample (sort @samples){
		$sum{$sample}{$key} = 0 unless defined($sum{$sample}{$key});
		print XLS "$sum{$sample}{$key}\t";
	}
	$sum{$project}{$key} = 0 unless defined($sum{$project}{$key});
	my $key_mult = "$key"."_mult";
	$sum{$project}{$key_mult} = 0 unless defined($sum{$project}{$key_mult});
	print XLS "$sum{$project}{$key}($sum{$project}{$key_mult})\n";
}

open (OUT, ">$results_dir/index.html");
print OUT "<HTML><BODY><FONT SIZE = +1>\n";
print OUT "<TABLE CELLPADDING = 5>\n";
print OUT "<TR><TH><CENTER><FONT SIZE = +1>Sample";
foreach my $sample (sort @samples){
	print OUT "<TH><CENTER><FONT SIZE = +1>$sample";
}
print OUT "<TH><CENTER><FONT SIZE = +1>ALL SAMPLES<BR>multiples</TR>\n";

foreach my $key (qw(Total Common SNVs Indels LOF Non-synonymous Non-coding COSMIC_mutations ClinVar_mutations CCLE_mutations DNA_repair_genes)) {
	print OUT "<TR><TH><CENTER><FONT SIZE = +1>$key";
	foreach my $sample (sort @samples){
		$sum{$sample}{$key} = 0 unless defined($sum{$sample}{$key});
		print OUT "<TD><CENTER><FONT SIZE = +1>$sum{$sample}{$key}";
	}
	$sum{$project}{$key} = 0 unless defined($sum{$project}{$key});
	my $key_mult = "$key"."_mult";
	$sum{$project}{$key_mult} = 0 unless defined($sum{$project}{$key_mult});
	print OUT "<TD><CENTER><FONT SIZE = +1>$sum{$project}{$key}($sum{$project}{$key_mult})</TR>\n";
}
print OUT "</TABLE>\n";

#print tables for calculating p-value for overlaps between different categories
my $cell2 = $sum{'cancer_mut'} - $sum{'lof_nonsyn_cancer_mut'};
my $cell3 = $sum{'lof_nonsyn'} - $sum{'lof_nonsyn_cancer_mut'};
my $cell4 = $sum{$project}{'Total'} - $sum{'lof_nonsyn_cancer_mut'} - $cell2 - $cell3;
print OUT "<P><HR>\n";
print OUT "<TABLE CELLPADDING = 7>\n";
print OUT "<TR><TH><CENTER><TH><CENTER><FONT SIZE = +1>LOF & nonsynonymous<TH><CENTER><FONT SIZE = +1>silent & non-coding</TR>\n";
print OUT "<TR><TH><CENTER><FONT SIZE = +1>Seen in cancer<TD><CENTER><FONT SIZE = +1>$sum{'lof_nonsyn_cancer_mut'}<TD><CENTER><FONT SIZE = +1>$cell2</TR>\n";
print OUT "<TR><TH><CENTER><FONT SIZE = +1>Not seen in cancer<TD><CENTER><FONT SIZE = +1>$cell3<TD><CENTER><FONT SIZE = +1>$cell4</TR>\n";
print OUT "</TABLE>\n";

$cell2 = $sum{'repair_genes'} - $sum{'lof_nonsyn_repair_genes'};
$cell3 = $sum{'lof_nonsyn'} - $sum{'lof_nonsyn_repair_genes'};
$cell4 = $sum{$project}{'Total'} - $sum{'lof_nonsyn_repair_genes'} - $cell2 - $cell3;
print OUT "<P><HR>\n";
print OUT "<TABLE CELLPADDING = 10>\n";
print OUT "<TR><TH><CENTER><TH><CENTER><FONT SIZE = +1>LOF & nonsynonymous<TH><CENTER><FONT SIZE = +1>silent & non-coding</TR>\n";
print OUT "<TR><TH><CENTER><FONT SIZE = +1>In repair genes<TD><CENTER><FONT SIZE = +1>$sum{'lof_nonsyn_repair_genes'}<TD><CENTER><FONT SIZE = +1>$cell2</TR>\n";
print OUT "<TR><TH><CENTER><FONT SIZE = +1>Not in repair genes<TD><CENTER><FONT SIZE = +1>$cell3<TD><CENTER><FONT SIZE = +1>$cell4</TR>\n";
print OUT "</TABLE>\n";

print OUT "<P><HR>\n";

print OUT "<P>Variants statistics in <A HREF = $project.stats.xls>XLS</A> format<BR>\n";
print OUT "<P>List of <A HREF = $project.multiples.xls>genes</A> with multipe LOF & non-synonymous mutations across all patients<BR>\n";
print OUT "<P>Annotated <A HREF = $project.annotated.maf>MAF</A> file<BR>\n";
print OUT "<P>Annotated <A HREF = $project.annotated.vcf>VCF</A> file<BR>\n";
print OUT "<P>Non annotated <A HREF = $project.maf>MAF</A> file<BR>\n";
print OUT "<P>Non annotated <A HREF = $project.vcf>VCF</A> file<BR>\n";

#deply results
system("chmod 660 $analysis_dir/index.html");
system("scp -r $results_dir/index.html $deployment_server:$summary_deployment/index.html");
system("scp -r $results_dir/$project.stats.xls $deployment_server:$summary_deployment/$project.stats.xls");
system("scp -r $results_dir/$project.annotated.maf $deployment_server:$summary_deployment/$project.annotated.maf");
system("scp -r $results_dir/$project.annotated.vcf $deployment_server:$summary_deployment/$project.annotated.vcf");
system("scp -r $results_dir/$project.maf $deployment_server:$summary_deployment/$project.maf");
system("scp -r $results_dir/$project.vcf $deployment_server:$summary_deployment/$project.vcf");
system("scp -r $results_dir/$project.multiples.xls $deployment_server:$summary_deployment/$project.multiples.xls");
system("ssh $deployment_server chmod -R 664 $summary_deployment/*");





