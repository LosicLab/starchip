#!/usr/bin/perl
use warnings;
use strict;

open MSA, $ARGV[0] or die $!;
my $count = 0;
my $seqcount =0;
my $sequence ='';
my @matrix;  

#load fasta into a matrix[rows][sequence]
while (my $msaline = <MSA> ){
	if ($msaline =~ m/^>/ ) { 
		if ( $count == 0 ) {
			next;
		}
		else {
			foreach my $y (0 ..length($sequence)-1) {
				$matrix[$seqcount][$y] = substr $sequence, $y, 1 ;
			}
			$seqcount++;  
			$sequence = ''; 
			next;
		}
	}
	else {
		chomp $msaline ; 
		$sequence .= $msaline; 
	}
	#print $msaline;
	$count++;
}
foreach my $y (0..length($sequence)-1) {
	$matrix[$seqcount][$y] = substr $sequence, $y, 1 ; 
}

#cycle through each position, top to bottom.  count the frequencies.  
foreach my $y ( 0..length($sequence)-1 ) {
	my %counthash = ();
	foreach my $x ( 0..$seqcount) {
		next if $matrix[$x][$y] eq '-';
		if (defined $counthash{$matrix[$x][$y]}) {
			$counthash{$matrix[$x][$y]} = $counthash{$matrix[$x][$y]} +1;
		}
		else {
			$counthash{$matrix[$x][$y]} = 1;
		}
	}
	my $max_key=();
	my $key2=();
	$max_key =(reverse sort {$counthash{$a} <=> $counthash{$b}} keys %counthash)[0];
	$key2 = (reverse sort {$counthash{$a} <=> $counthash{$b}} keys %counthash)[1];
	if (defined $key2) {
		if (($counthash{$max_key}+1) / ($counthash{$key2}+1) > 5 ){
		print "$max_key";
		}
		else { print "*";}
	}
	else {print "$max_key";}
	#print "$y: $max_key : $counthash{$max_key}\n   $key2 : $counthash{$key2}\n";
}
print "\n";
