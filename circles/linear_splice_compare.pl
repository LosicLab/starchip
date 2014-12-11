#!/usr/bin/perl
use strict;
use warnings;

#usage : linear_splice_compare.pl cRNA SampleID SJ.out.tab joinstrands > output.file

open (cRNA, "<$ARGV[0]") or die $!;
my $SpJu = "star_" . $ARGV[1] . "/SJ.out.tab";
open (SJ, "<$SpJu") or die $!;
	my @sjout = <SJ>;
close SJ or die $!;
my $joinedstrands = "/sc/orga/projects/losicb01a/common_folder/starnet/100bp/cRNA/joinstrands.star_" . $ARGV[1] ; 
print "CircleID\tLinearCount\tType\tCircularCount\n";
while (my $circ = <cRNA> ) {
#go through each circle candidate, calculate the linear splices at that location, the circular splices, then gene annotation, the circles at that gene, and the total linear splices at that gene. 

###Linear splices at the cRNA sites###
	my @circle = split(/\s+/, $circ); # 0: chrm 1:p1 (smaller) 2:p2 (larger) 3:other
	my $chr = $circle[0]; my $p1 = $circle[1]; my $p2 = $circle[2];
	#search for linear splices that end with p1 or start with p2. 
	#because SJ.out and circs are both  organized low to high, we don't have to worry about strand.
	my @p1hits = grep(/$chr\t.*\t$p1\t/, @sjout);
	my @p2hits = grep(/$chr\t$p2\t/, @sjout);
	#print @p1hits;print "p2hits:\n";print @p2hits;
	my $p1unique = 0; my $p1multi = 0;my $p2unique = 0;my $p2multi = 0;
	if (@p1hits) {
		foreach (@p1hits) {
			my @hitx = split(/\s+/, $_);
			$p1multi += $hitx[7];$p1unique += $hitx[6];
		}
	}
	if (@p2hits) {
               foreach (@p2hits)	{
                        my @hitx = split(/\s+/, $_);
                        $p2multi += $hitx[7];$p2unique += $hitx[6];
                }
        }
	my $p1sum=$p1multi + $p1unique;
	my $p2sum=$p2multi + $p2unique;
	#print "p1u:$p1unique p1m:$p1multi p2u:$p2unique p2m:$p2multi\n";
	my $type;
	my $linearcount; 
	if ($p1sum+1 > (3*$p2sum)+1) { #if only hits that end in p1 are there, probably the last exon
		$type ="E" ;
		$linearcount = $p1sum;
	}
	elsif ($p2sum+1 > (3*$p1sum)+1) { #if only hits that start in p2 are there, probably the first exon
		$type ="S" ;
		$linearcount = $p2sum;
	}
	else {
		$type ="M" ;
		$linearcount = ($p1sum + $p2sum)/2;
	}
	print "$circle[0]-$circle[1]-$circle[2]-$circle[3]\t$linearcount\t$type\t"; 
####Cirucular Splices at the site.  
	my $circlereadscommand = "grep -P \"\\t$chr\\t$p2\\t.\\t$chr\\t$p1\" $joinedstrands |cut -f1";
	#print "$circlereadscommand\n";
	my $circlecount =`$circlereadscommand`;
	chomp($circlecount);
	if ($circlecount eq "") {
		$circlecount =0;
	}			
	print "$circlecount\n";
}



