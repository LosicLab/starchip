#! /usr/bin/perl
use strict;
use warnings; 

## Take in a list of cRNA, and file.  Then search that file for those cRNA, outputing either the value found or a 0 if not found.
#inputs
my $cRNAfile = $ARGV[0];
my $searchFile = $ARGV[1];
my $matrixType = $ARGV[2]; 
my $outfile = $searchFile . ".grepcircles.temp" ; 
my $outfile2 = $searchFile . ".grepcircles.temp2" ;
#initiate vars
my @crna; 
my $columnindex;
my $columnindex2; 
my $id = $searchFile ;
$id =~ s/rawdata\///; 
$id =~ s/backsplices.//;
$id =~ s/\.[0-9]+\.spliced//; 

#differentiate count matrix vs linear splices
if ($matrixType eq "count") {
	$columnindex = 0 ; 
}
else {
	$columnindex = 26;
	$columnindex2 = 27; 
}

#store crna in memory
open cRNA, "<$cRNAfile" or die $!; 
while (my $x = <cRNA> ) {
	my ($chr, $p1, $p2) = split(/\s+/, $x); 
	my $crna_line = $chr . ":" . $p1 . "-" . $p2 ; 
	push @crna, $crna_line;
}
close cRNA; 

#create output file(s)
open OUT, ">$outfile" or die $!;
	print OUT "$id\n";
if ($matrixType ne "count") {
	open OUT2, ">$outfile2" or die $!;
	print OUT2 "$id\n"; 
}

#store file to search in memory
open SEARCH, "<$searchFile" or die $!; 
my @searchlines = <SEARCH> ;

#cycle through crna
foreach my $circ (@crna) {
	my @matches = grep { /$circ/ } @searchlines ; 
	if (scalar @matches == 1) {
		my @matchsplit = split(/\s+/, $matches[0]); 
		print OUT "$matchsplit[$columnindex]\n";
		if ($matrixType ne "count") {
			print OUT2 "$matchsplit[$columnindex2]\n";
		}
	}
	elsif (scalar @matches == 0) {
		print OUT "0\n"; 
		if ($matrixType ne "count") {
                        print OUT2 ".\n";
                }
        }
}
	

