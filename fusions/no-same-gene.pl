#! /usr/bin/perl
use warnings;
use strict;

#usage no-same-gene.pl input >output

#remember perl counts from 0
#my $a1col = 19;
#my $a1col = 20;
#my $a2col = 22;
#my $a2col = 23;

open FILE, "<$ARGV[0]" or die $!;

while (my $x = <FILE>) {
	my @line=split(/\s+/, $x);

	my @geneannot=grep(/gene_id/, @line);  

	my @anno1=split(/;/, $geneannot[0]);
	my @anno2=split(/;/, $geneannot[1]);
	#my @anno1=split(/;/, $line[$a1col]);
	#my @anno2=split(/;/, $line[$a2col]);
	my @gene1name=grep(/gene_name/, @anno1);
	my @gene2name=grep(/gene_name/, @anno2);
	if ($gene1name[0] eq $gene2name[0]) {
		#print "@gene1name\t@gene2name\n";
		next;
	}
	print "$x";
}
