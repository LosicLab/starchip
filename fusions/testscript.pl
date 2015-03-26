#!/usr/bin/perl
use strict;
use warnings;


open ABS, "<data/hg19.abparts" or die $!;
my @AbParts;
my $index=0; 
while (my $x = <ABS>) {
	my @line = split(/\s+/, $x); #AbParts lines are formatted: chrm	pos1 pos2
	$AbParts[$index][0]{$line[0]} = $line[1]; ##data format: $AbParts[index][chrom][0=lower/1=upper]=position
	$AbParts[$index][1]{$line[0]} = $line[2]; ##data format: $AbParts[index][chrom][0=lower/1=upper]=position
	$index++;
}
print "$index total ab lines\n";
open TEST, "<testfusions.txt" or die $!;
EXITHERE: while (my $x = <TEST>) {
	my ($chr1, $pos1, $chr2, $pos2) = split(/\s+/, $x);
	print "$chr1:$pos1 -- $chr2:$pos2\n";
	foreach my $z (0..$index) {
		if ($AbParts[$z][0]{$chr1} <= $pos1 && $pos1 <= $AbParts[$z][1]{$chr1}) {
			print "pos1 match\n";
			foreach my $y (0..$index) {
				if ($AbParts[$y][0]{$chr2} <= $pos2 && $pos2<= $AbParts[$y][1]{$chr2}) {
					print "pos2 match\n";
					next EXITHERE;  
				}
			}
		}
	}
	print "continuing with $chr1:$pos1 -- $chr2:$pos2\n";
}
