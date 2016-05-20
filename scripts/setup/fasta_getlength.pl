#!/usr/bin/perl
use warnings;
use strict;


my $length = 0;
my $ID; 

while (my $x = <>) {
	if ($x =~ m/^>/) {
		if ($length > 0) {
			print "$ID\t$length\n";
		}
		chomp $x; 
		$ID=$x;
		$ID =~ s/>//;
		$ID =~ s/\s.*//; 
		$length = 0;

	}
	else {
		$length += length($x);
	}
}
chomp $ID;
	print "$ID\t$length\n";

