#!/usr/bin/perl -w

use strict;
use warnings;

my $bedtools_version = "#bedtoolsVersion";
my $dir = "#bdDir";
my $genes_bed = "#genesBed";
my $project = "#project";
my $deployment_server = "#deploymentServer";
my $summary_deployment = "#summaryDeployment/";

#my $bedtools_version = "2.13.3";
#my $dir = "/groupvol/cgi/results/johnson_glioma/breakDancer/2015-02-11";
#my $genes_bed = "/groupvol/cgi/resources/annotations/eukaryote/human/GRCh37/GRCh37.74.gene.bed";
#my $project = "johnson_glioma";
#my $deployment_server = "eliot.med.ic.ac.uk";
#my $summary_deployment = "/www/html/report/project/johnson_glioma/breakDancer/2015-02-11";

my %blood = (); 
my %tumor = ();

#bed file with all translocations
my $all_bed = "$dir/multisample/$project".".all.bed";
open(ALL_BED, ">$all_bed");
print ALL_BED "browser position chr1\n";
print ALL_BED "browser hide all\n";

#collect data about breakpoints for all girmline samples
my($blood_sample, $tumor_sample, $sample_pair);
my %blood_samples = ();
my %tumor_samples = ();

opendir (BD, "$dir") || print "Can't open $dir\n";
while (defined($sample_pair = readdir(BD))){

        next unless $sample_pair =~ /(.*)\.vs\.(.*)/;
	$blood_sample = $1;
	$tumor_sample = $2;
	if (-e "$dir/$sample_pair/$blood_sample".".bd" && -e "$dir/$sample_pair/$tumor_sample".".bd"){
	        $blood_samples{$sample_pair} = $blood_sample;
	        $tumor_samples{$sample_pair} = $tumor_sample;
	}
}

my($bd_blood, @blood_data, $blood_bp1, $blood_ext1, $blood_bp2, $blood_ext2);

foreach  $sample_pair (keys %blood_samples){

        $blood_sample = $blood_samples{$sample_pair};
	print ALL_BED "track name=$blood_sample type=bed visibility=3 db=hg19 itemRgb=On\n";

	$bd_blood = "$dir/$sample_pair/$blood_sample".".bd";
	open(BLOOD, "$bd_blood");
	while (<BLOOD>){

		next if /^#|MT|GL|NC_007605|hs37d5/;
		@blood_data = split(/\t/,$_);

		#collect all germlime breakpoints into %blood
		$blood{"$blood_data[0]:$blood_data[3]"}{"$blood_data[1]:$blood_data[4]"}++;
		$blood_bp1 = $blood_data[1]; $blood_ext1 = $blood_bp1 + 100;
		$blood_bp2 = $blood_data[4]; $blood_ext2 = $blood_bp2 + 100;

                #print all girmline translocations into bed file
		print ALL_BED "chr$blood_data[0]\t$blood_bp1\t$blood_ext1\tCTX:chr$blood_data[3]:$blood_bp2\t0\t+\t$blood_bp1\t$blood_ext1\t0,0,255\n";
		print ALL_BED "chr$blood_data[3]\t$blood_bp2\t$blood_ext2\tCTX:chr$blood_data[0]:$blood_bp1\t0\t+\t$blood_bp2\t$blood_ext2\t0,0,255\n";

	}
}


#collect data about SOMATIC breakpoints for all samples
my($bd_tumor, @tumor_data, $chrom_pair, $tumor_bp1, $tumor_ext1, $tumor_bp2, $tumor_ext2, $f1, $blood_coord, @blood_coord, $f2, $tumor_coord, @tumor_coord);

foreach  $sample_pair (keys %tumor_samples){

        $tumor_sample = $tumor_samples{$sample_pair};
	print ALL_BED "track name=$tumor_sample type=bed visibility=3 db=hg19 itemRgb=On\n";

	$bd_tumor = "$dir/$sample_pair/$tumor_sample".".bd";
	open(TUMOR, "$bd_tumor");
	while (<TUMOR>){

		next if /^#|MT|GL|NC_007605|hs37d5/;
		@tumor_data = split(/\t/,$_);

		#find out whether this breakpoind is within 1kb of any girmline event from any sample 
		$chrom_pair = "$tumor_data[0]:$tumor_data[3]";
		$tumor_bp1 = $tumor_data[1]; $tumor_ext1 = $tumor_bp1 + 100;
		$tumor_bp2 = $tumor_data[4]; $tumor_ext2 = $tumor_bp2 + 100;

		$f1 = 0;
		foreach $blood_coord (keys %{$blood{$chrom_pair}}){
			@blood_coord = split(/:/, $blood_coord);
			$blood_bp1 = $blood_coord[0]; $blood_bp2 = $blood_coord[1];
			$f1++ if ((($tumor_bp1 > $blood_bp1 - 1000) && ($tumor_bp1 < $blood_bp1 + 1000)) && (($tumor_bp2 > $blood_bp2 - 1000) && ($tumor_bp2 < $blood_bp2 + 1000)));
			last if $f1;
		}

		#if event is somatic
		unless ($f1){

		        #print breakpoint info into bed file
			print ALL_BED "chr$tumor_data[0]\t$tumor_bp1\t$tumor_ext1\tCTX:chr$tumor_data[3]:$tumor_bp2\t0\t+\t$tumor_bp1\t$tumor_ext1\t255,0,0\n";
			print ALL_BED "chr$tumor_data[3]\t$tumor_bp2\t$tumor_ext2\tCTX:chr$tumor_data[0]:$tumor_bp1\t0\t+\t$tumor_bp2\t$tumor_ext2\t255,0,0\n";

			#collect all somatic breakpoints into %tumor
			#if there is a somatic breakpoint within 1 kb in another sample, save it under the same coordinate key (will be usefull while looking for breakpoints that come up in at least 2 samples)
			$f2 = 0;
			foreach $tumor_coord (keys %{$tumor{$chrom_pair}}){
				@tumor_coord = split(/:/, $tumor_coord);
				if ((($tumor_bp1 > $tumor_coord[0] - 1000) && ($tumor_bp1 < $tumor_coord[0] + 1000)) && (($tumor_bp2 > $tumor_coord[1] - 1000) && ($tumor_bp2 < $tumor_coord[1] + 1000))){
					$tumor{$chrom_pair}{$tumor_coord}{$tumor_sample} = "chr$tumor_data[0]\t$tumor_bp1\t$tumor_ext1\tCTX:chr$tumor_data[3]:$tumor_bp2\t0\t+\t$tumor_bp1\t$tumor_ext1\t255,0,0\nchr$tumor_data[3]\t$tumor_bp2\t$tumor_ext2\tCTX:chr$tumor_data[0]:$tumor_bp1\t0\t+\t$tumor_bp2\t$tumor_ext2\t255,0,0\n";
					$f2++;
					last;
				}
			}

			#create new entry only if there is no already another somatic breakpoint within 1 kb
			$tumor{$chrom_pair}{"$tumor_bp1:$tumor_bp2"}{$tumor_sample} = "chr$tumor_data[0]\t$tumor_bp1\t$tumor_ext1\tCTX:chr$tumor_data[3]:$tumor_bp2\t0\t+\t$tumor_bp1\t$tumor_ext1\t255,0,0\nchr$tumor_data[3]\t$tumor_bp2\t$tumor_ext2\tCTX:chr$tumor_data[0]:$tumor_bp1\t0\t+\t$tumor_bp2\t$tumor_ext2\t255,0,0\n" unless $f2;
		}
	}
}

#collect somatic breakpoints that come up in > 2 samples
my %count = ();
my ($coord, $count, $sample);
foreach $chrom_pair (keys %tumor){
	foreach $coord (keys %{$tumor{$chrom_pair}}){
		$count = 0;
		foreach $sample_pair (keys %tumor_samples){
		    $tumor_sample = $tumor_samples{$sample_pair};
			$count++ if defined($tumor{$chrom_pair}{$coord}{$tumor_sample});
		}
		$count{$chrom_pair}{$coord} = $count if $count > 1;
	}
}

#bed file with somatic breakpoints that come up in > 2 samples
my $bed_mult = "$dir/multisample/$project".".mult.bed";
open(MULT_BED, ">$bed_mult");
print MULT_BED "browser position chr1\n";
print MULT_BED "browser hide all\n";

#bed file with somatic breakpoints that come up in > 2 samples for annotation with bedtools
my $bed_ann = "$dir/multisample/$project".".mult_for_ann.bed";
open(ANN_BED, ">$bed_ann");

foreach $sample_pair (keys %tumor_samples){

        $tumor_sample = $tumor_samples{$sample_pair};
	print MULT_BED "track name=$tumor_sample"."_mult type=bed visibility=3 db=hg19 itemRgb=On\n";
	print ANN_BED "track name=$tumor_sample"."_mult type=bed visibility=3 db=hg19 itemRgb=On\n";
	foreach $chrom_pair (keys %count){
		foreach $coord (keys %{$count{$chrom_pair}}){	

		        next unless defined($tumor{$chrom_pair}{$coord}{$tumor_sample});
			print MULT_BED "$tumor{$chrom_pair}{$coord}{$tumor_sample}";

			$tumor{$chrom_pair}{$coord}{$tumor_sample} =~ s/chr//g;
			$tumor{$chrom_pair}{$coord}{$tumor_sample} =~ s/\t0\t/\t$count{$chrom_pair}{$coord}\t/g;
			print ANN_BED "$tumor{$chrom_pair}{$coord}{$tumor_sample}";

		}
	}
}

#find genes that overlap with breakpoints
my $bed_ann_genes = "$dir/multisample/$project".".mult_annotated.bed";
`/apps/bedtools/$bedtools_version/bin/intersectBed -a $bed_ann -b $genes_bed -wa -wb > $bed_ann_genes`;

#make HTML and xls table bp1: chrom, coord, gene name, gene ID, bp2: chrom, coord, gene name, gene ID
my $html = "$dir/multisample/$project".".html";
open(OUT, ">$html");
print OUT "<HTML><BODY>\n";
print OUT "<H2><CENTER>Genes that overlap with the same somatic inter-chromosomal translocation breakpoints in multiple samples</H2>\n";
print OUT "<TABLE CELLPADDING = 5>\n";
print OUT "<TR><TH><CENTER><FONT SIZE = +1>Sample<TH><CENTER><FONT SIZE = +1>Breakpoint 1<TH COLSPAN=3><CENTER><FONT SIZE = +1>Overlapping Genes<TH><CENTER><FONT SIZE = +1>Breakpoint 2<TH COLSPAN=3><CENTER><FONT SIZE = +1>Overlapping Genes\n";
print OUT "<TR><TH><TH><TH><FONT SIZE = +1>Gene ID<TH><FONT SIZE = +1>Gene Name<TH><FONT SIZE = +1>Biotype<TH><TH><FONT SIZE = +1>Gene ID<TH><FONT SIZE = +1>Gene Name<TH><FONT SIZE = +1>Biotype\n";
print OUT "<TR><TD COLSPAN=9><HR>\n";

my $xls = "$dir/multisample/$project".".xls";
open(XLS, ">$xls");
print XLS "Sample\tBreakpoint 1\tOverlapping genes\tBreakpoint 2\tOverlapping genes\n";

my(@bp2, @chrom_pair, @coord_pair, $string1, $genes1, $string2, $genes2);
foreach $chrom_pair (keys %count){
	foreach $coord (keys %{$count{$chrom_pair}}){	
	        foreach $sample_pair (keys %tumor_samples){

		        $tumor_sample = $tumor_samples{$sample_pair};
		        next unless defined($tumor{$chrom_pair}{$coord}{$tumor_sample});
			@tumor_data = (); @bp2 = (); @chrom_pair = (); @coord_pair = (); $string1 = ""; $genes1 = ""; $string2 = ""; $genes2 = ""; 
			@tumor_data = split(/\t/, $tumor{$chrom_pair}{$coord}{$tumor_sample});
			@bp2 = split(/:/,$tumor_data[3]);

			$string1 = "$tumor_data[0]\t$tumor_data[1]";
			$genes1 = `cat $bed_ann_genes|grep '$string1'|cut -f 13|sort|uniq`;
			chop($genes1);
			$genes1 =~ s/\n/<BR>/g;
			$genes1 =~ s/\;/ /g;

			$string2 = "$bp2[1]\t$bp2[2]";
			$genes2 = `cat $bed_ann_genes|grep '$string2'|cut -f 13|sort|uniq`;
			chop($genes2);
			$genes2 =~ s/\n/<BR>/g;
			$genes2 =~ s/\;/ /g;

			print OUT "<TR><TD>$tumor_sample<TD>chr$tumor_data[0]:$tumor_data[1]<TD COLSPAN = 3>$genes1<TD>chr$bp2[1]:$bp2[2]<TD COLSPAN = 3>$genes2\n";


			$genes1 =~ s/<BR>/;/g;
			$genes2 =~ s/<BR>/;/g;
			print XLS "$tumor_sample\tchr$tumor_data[0]:$tumor_data[1]\t$genes1\tchr$bp2[1]:$bp2[2]\t$genes2\n";

		}
		print OUT "<TR><TD COLSPAN=9><HR>\n";
	}
}
print OUT "</TABLE>\n";
print OUT "<P><FONT SIZE = +1>Download table in <A HREF = $project.mult.xls>XLS</A> format<BR>\n";
print OUT "<P><FONT SIZE = +1>Germline and somatic inter-chromosomal translocation events in <A HREF = $project.all.bed>BED</A> format<BR>\n";
print OUT "<P><FONT SIZE = +1>Somatic inter-chromosomal translocation events discovered in multiple samples in <A HREF = $project.mult.bed>BED</A> format<BR>\n";
print OUT "</BODY></HTML>\n";

system("chmod 660 $dir/multisample/*");
system("scp -r $html $deployment_server:$summary_deployment/index.html");
system("scp -r $xls $deployment_server:$summary_deployment/$project.mult.xls");
system("scp -r $all_bed $deployment_server:$summary_deployment/$project.all.bed");
system("scp -r $bed_mult $deployment_server:$summary_deployment/$project.mult.bed");
system("ssh $deployment_server chmod -R 664 $summary_deployment/*");

