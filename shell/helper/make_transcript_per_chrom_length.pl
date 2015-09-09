#!/usr/bin/perl -w

$tmp_dir = $ARGV[0];
$rf = "$tmp_dir/tmp.rf";
$lf = "$tmp_dir/tmp.lf";

open(RF, "$rf");
open(LF, ">$lf");
print LF "chrom\ttranscript_length\n";

%data=();

while(<RF>){
	@data = split(/\t/, $_);
	$chrom = $data[2];
	@start = split(/\,/, $data[9]);
	@end = split(/\,/, $data[10]);
	$i = 0;
	foreach $start (@start) {
		next unless $start =~ /\d+/;
		$data{$chrom}{$start}{$end[$i]}++;
		$i++;
	}
}

$gl_length = 0;
foreach $chrom (sort {$a cmp $b} keys %data){
	$prev_start = 0;
	$prev_end = 0;
	$length = 0;
	
	foreach $start (sort {$a <=> $b} keys %{$data{$chrom}}){
		foreach $end (sort {$a <=> $b} keys %{$data{$chrom}{$start}}){
			if ($prev_end < $start) {
				$length += $prev_end - $prev_start if $prev_end > 0;
				$prev_start = $start;
				$prev_end = $end;
			} else {
				$prev_end = $end;
			}
		}
	}
	$length += $prev_end - $prev_start;

	if ($chrom =~ /GL/) {
		$gl_length += $length;
	} else {
		print LF "$chrom\t$length\n";
	}
}
print LF "GL\t$gl_length\n";
