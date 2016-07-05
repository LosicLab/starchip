#!/usr/bin/perl
use strict;
use warnings;

#check if the positions in file1 are within Z bp of file1.  if so, prints just the first to appear.  

my $distancemax = 10;

my %cRNA_hash;
# cRNA_hash{file1ID}{chrom}[lower/upper==0/1] = position

BIGLOOP: while (my $x = <>) {
	chomp $x;
	my ($chrom, $pos1, $pos2) = split(/[:\-\t]/, $x);
	foreach my $cRNAIDs (keys %cRNA_hash) {
		if (exists $cRNA_hash{$cRNAIDs}{$chrom}[0]) {
			if (abs($cRNA_hash{$cRNAIDs}{$chrom}[0] - $pos1) <= $distancemax && abs($cRNA_hash{$cRNAIDs}{$chrom}[1] - $pos2) <= $distancemax ) {
				next BIGLOOP;
			}
		}
	}
	$cRNA_hash{$x}{$chrom}[0]=$pos1;
	$cRNA_hash{$x}{$chrom}[1]=$pos2;
	print "$x\n"; 
}


