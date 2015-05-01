#!/usr/bin/perl -w

use strict;
use warnings;

my $dir = "#bdDir";
my $results_dir = "#resultsDir";

my %blood = (); 
my %tumor = ();

#bed file with all translocations
my $all_bed = "$results_dir/all".".bed";
open(ALL_BED, ">$all_bed");
print ALL_BED "browser position chr1\n";
print ALL_BED "browser hide all\n";

#collect data about breakpoints for all girmline samples
my($blood_sample, $bd_blood, @blood_data, $blood_bp1, $blood_ext1, $blood_bp2, $blood_ext2);

opendir (BD, "$dir") || print "Can't open $dir\n";
while (defined(my $sample_pair = readdir(BD))){

        next if $sample_pair =~ /multisample/;
	$blood_sample = $1 if $sample_pair =~ /(.*)\.vs\..*/;
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
my($tumor_sample, $bd_tumor, @tumor_data, $chrom_pair, $tumor_bp1, $tumor_ext1, $tumor_bp2, $tumor_ext2, $f1, $blood_coord, @blood_coord, $blood_bp1, $blood_bp2, $f2, $tumor_coord, @tumor_coord);
my @tumor_samples = ();

opendir (BD, "$dir") || print "Can't open $dir\n";
while (defined(my $sample_pair = readdir(BD))){

        next if $sample_pair =~ /multisample/;
	$tumor_sample = $1 if $sample_pair =~ /(.*)\.vs\..*/;
	push(@tumor_samples,$tumor_sample);
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
			#if there is a somatic breakpoint within 1 kb in another sample, save it under the same coordinate key (will be usefull while looking for breakpoints that come un in > 2 samples)
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
my ($chrom_pair, $coord, $count, $sample);
foreach $chrom_pair (keys %tumor){
	foreach $coord (keys %{$tumor{$chrom_pair}}){
		$count = 0;
		foreach $tumor_sample (@tumor_samples){
			$count++ if defined($tumor{$chrom_pair}{$coord}{$tumor_sample});
		}
		$count{$chrom_pair}{$coord} = $count if $count > 1;
	}
}

#bed file with somatic breakpoints that come up in > 2 samples
my $bed_mult = "$results_dir/mult".".bed";
open(MULT_BED, ">$bed_mult");
print MULT_BED "browser position chr1\n";
print MULT_BED "browser hide all\n";

#bed file with somatic breakpoints that come up in > 2 samples for annotation with bedtools
my $bed_ann = "$results_dir/mult_for_ann.bed";
open(ANN_BED, ">$bed_ann");

foreach $tumor_sample (@tumor_samples){
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


