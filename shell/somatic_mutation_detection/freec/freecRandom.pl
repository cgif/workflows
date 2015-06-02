#!/usr/bin/perl -w

use strict;
use warnings;

my $freec_dir = "/groupvol/cgi/results/johnson_glioma/FREEC/2015-04-22";
my $chrom_length = "/home/asoskins/workflows/shell/somatic_mutation_detection/freec/GRCh37.len";
my $exon_bed = "/groupvol/cgi/resources/annotations/eukaryote/human/GRCh37/GRCh37.74.exon.bed";
my $deployment_server = "eliot.med.ic.ac.uk";
my $summary_deployment = "/www/html/report/project/johnson_glioma/FREEC/2015-02-22";

#my $freec_dir = "#resultsDir";
#my $chrom_length = "#chrLenFile";
#my $exon_bed = "#exonBed";
#my $deployment_server = "#deploymentServer";
#my $summary_deployment = "#summaryDeployment";

my $index_file = "$freec_dir/multisample/index.html"
#my $randomizations = 1000;
my $randomizations = 5;
my %rand = ();
my %real = ();

open (COUNTS,"$index_file");
while(<COUNTS>){ 

	chmod();
	$real{$1} = $2 if /<TR><TD>(\d+)<TD><A HREF = .*>(\d+)<\/A>/;

}

foreach $key (keys %real){
	print "$key - $real{$key}\n";
}

my $chrom_start = 1;
my $chrom_end = 0;
open (CHROM, "$chrom_length");
open (CHROM_RAND, ">$TMPDIR/chrom_coord.txt");
while(<CHROM>){

	chmod();
	@data = split(/\t/, $_);
	next if $data[1] =~ /X|Y/;
	$chrom_end = $chrom_start + $data[2]
	print CHROM_RAND "$data[1]\t$chrom_start\t$chrom_end\n";
	$chrom_start = $chrom_end + 1;

}
my $genome_size = $chrom_end;
print "genome size $genome_size\n";

open (CHROM_RAND, "$TMPDIR/chrom_coord.txt");
while(<CHROM_RAND>){
	print;
}

foreach my $i (1..$randomizations) {

	my $cnv_bed = "$TMPDIR/randomCNVs.bed";
	my $cnv_ann = "$TMPDIR/randomCNVs.ann";
	my $cnv_ann_int1 = "$TMPDIR/randomCNVs.ann1";
	my $cnv_ann_int2 = "$TMPDIR/randomCNVs.ann2";
	my $cnv_ann_int3 = "$TMPDIR/randomCNVs.ann3";
	my $cnv_gene = "$TMPDIR/randomCNVs.gen";

	#pick up random locations for CNVs and print them into bed file
	open(CNV_BED, ">$cnv_bed") || print "Can't open $cnv_bed\n";
	opendir (FREEC_DIR, "$freec_dir") || print "Can't open $freec_dir\n";
	while (defined(my $sample = readdir(FREEC_DIR))){	

		my $FREEC_file = "$freec_dir/$sample/tumor.bam_CNVs.p.value.txt";
		#my $FREEC_file = "$freec_dir/$sample/$sample.CNVs.txt";
		next unless (-e $FREEC_file);

		my $sample_name = $1 if $sample =~ /(.*?)\..*/;
		my %int = ();

		open (FREEC_CNV, "$FREEC_file") || print "Can't open $FREEC_file\n";
		while(<FREEC_CNV>){
			chomp(); 
			next if /^(chr|X|Y)/;

			my @data = split(/\t/,$_); 
			my $start = 0;
			my $end = 0;
			my $length = $data[2] - $data[1]; 
			my $chr = "";
			my $start_bed = 0;
			my $end_bed = 0;
			my $copy = $data[4];

			RAND:while (1){
				OVER:while (1){
					$start = int(rand($genome_size)); 
					$end = $start + $length;
					my $f_over = 0;
					foreach my $int_start (keys %int){
						foreach my $int_end (keys %{$int{$int_start}}){
							$f_over++ if (($start >= $int_start && $start <= $int_end) || ($end >= $int_start && $end <= $int_end));
						}
					}
					last OVER unless $f_over;
				}
				open (DICT, "$TMPDIR/chrom_coord.txt") || print "Can't open $TMPDIR/chrom_coord.txt\n";
				while (my $coord = <DICT>){
					my @coord = split(/\t/, $coord);
					if ($start >= $coord[1] && $start <= $coord[2]){
						$chr = "$coord[0]";
						$start_bed = $start - $coord[1];
						$end_bed = $start_bed + $length;
						if ($end < $coord[2]){
							$int{$start}{$end}++;
							last RAND;
						}
					}
				}
			}
			print CNV_BED "$chr\t$start_bed\t$end_bed\t$sample_name:$copy\n";
		}
	}

	#run bedtools to select genes overlapping CNVs
	`/apps/bedtools/2.13.3/bin/intersectBed -a $cnv_bed -b $exon_bed -wa -wb > $cnv_ann`;

	#print list of genes that were affected in > 2 patients
	open(CNV_ANN, "$cnv_ann") || print "Can't open $cnv_ann\n";
	open(CNV_ANN_INT1, ">$cnv_ann_int1") || print "Can't open $cnv_ann_int1\n";
	while(<CNV_ANN>){
		chomp();
		my @data = split(/\t/,$_);
		$data[7] =~ s/.*gene_name \"(.*?)\".*/$1/; 
		print CNV_ANN_INT1 "$data[3]\t$data[7]\n";
	}
	
	`cat $cnv_ann_int1| sort | uniq > $cnv_ann_int2`;

	open(CNV_ANN_INT2, "$cnv_ann_int2") || print "Can't open $cnv_ann_int2\n";
	open(CNV_ANN_INT3, ">$cnv_ann_int3") || print "Can't open $cnv_ann_int3\n";
	while(<CNV_ANN_INT2>){
		chomp();
		my @data = split(/\t/,$_);
		$data[0] =~ s/(.*):(.*)/$2/; 
		print CNV_ANN_INT3 "$data[0]\t$data[1]\n";
	} 

	`cat $cnv_ann_int3| sort | uniq -c  > $cnv_gene`;

	#collect randomization data into table
	open(CNV_GEN, "$cnv_gene") || print "Can't open $cnv_gene\n";
	while(<CNV_GEN>){
		chomp();
		s/^\s+(.*)/$1/;
		my @data = split(/\s/,$_);
		$rand{$data[0]}{$i} ++;
	}
}

`rm $TMPDIR/randomCNVs*`;

my %pval = ();
foreach my $count (sort {$a <=> $b} keys %real) {
	print "$count:";
	$pval{$count} = 0;
	foreach my $i (1..$randomizations) {
		$rand{$count}{$i} = 0 unless exists($rand{$count}{$i});
		$pval{$count}++ if $rand{$count}{$i} >= $real{$count};
		print "$real{$count}|$rand{$count}{$i}|$pval{$count}\t";
	}
	print "\n";
}

print "\n\n";

my $pval;
my $f1 = 0;
open (INDEX, ">$index_file");
open (INDEX_PVAL, "$TMPDIR/index.html");
while (<INDEX>){

	$f1++ if /#Samples/; 
	$f1 = 0 if /TABLE/;
	if ($f1){
		if (/#Samples/){
			print "<TR><TH><FONT SIZE = +1>#Samples<TH><FONT SIZE = +1>#Genes<TH><FONT SIZE = +1>P-value\n";
		}elsif (/(<TR><TD>(\d+)<TD><A HREF = .*>\d+<\/A>)/){
			$count = $2;
			$pval = $pval{$count}/$randomizations;
			print "$1<TD>$pval\n";
		}
	}else{
		print INDEX_PVAL;
	}

}

