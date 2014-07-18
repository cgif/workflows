#!/usr/bin/perl -w

$tmp_dir = $ARGV[0];
$gff = "$tmp_dir/tmp.gff";
$out = "$tmp_dir/tmp.length";
open(OUT, ">$out");

#collect gene ids from gff file
$ids = `cut -f 9 $gff|cut -f 1 -d ";"|sort|uniq`;
@ids = split(/\n/, $ids);

foreach $id (@ids){ 
#collect exon coordinates for the gene
    $exons = ""; @exons = (); %exons = ();
    $exons = `grep '$id' $gff|awk '\$3 == \"exon\"'|cut -f 4,5`;
    @exons = split(/\n/, $exons);

#remove shorter exons if they heve the same start coordinates
    foreach (@exons){
	/(\d+)\s+(\d+)/;
	$exon_start = $1; $exon_end = $2;
	$exons{$exon_start} = $exon_end unless defined($exons{$exon_start}) && $exons{$exon_start} > $exon_end;
    }

#merge overlapping exons and calculate gene length
    %remove = (); $sum = 0;
    foreach $start (sort {$a <=> $b} keys %exons){
	next if defined($remove{$start});
	$end = $exons{$start};
	foreach $start2 (sort {$a <=> $b} keys %exons){
	    if ($start2 > $start && $start2 <= $end){
		$remove{$start2}++;
		$end = $exons{$start2} if $exons{$start2} > $end;
	    }elsif ($start2 > $end){
		last;
	    }
	}
	$sum = $sum + ($end - $start);
    }

#print results to the file
    $id =~ s/.*\"(\w+)\"/$1/;
    print OUT "$id\t$sum\n";
}

