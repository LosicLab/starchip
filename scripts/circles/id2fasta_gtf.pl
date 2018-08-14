#!/usr/bin/perl
use warnings;
use strict;

my $crnasource = "circRNA";
my $exon_halfsize= 1; 
my $ref = $ARGV[1] ; 
my $chrflag = 0; #1 if remove chr from cRNA ids.
my $offsetval = 10000; 
my $nstring = "n"x$offsetval ;
my $sequence = $nstring;
my $chrmpos = 1 + length($nstring); 
open CIRCLES, "<$ARGV[0]" or die $!;
#my $saffile = $ARGV[0] . ".SAF" ; 
my $gtffile = $ARGV[0] . ".gtf" ; 
open GTF, ">$gtffile" or die $!;
print ">cRNA_custom\n";
while (my $y = <CIRCLES>) {
	next if ( $y eq 'ID');
	chomp $y; 
	my @splitted = split (/[:-]/, $y);
	##print "$splitted[0]\t$splitted[1]\t$splitted[2]\n";
	my $chr = $splitted[0];
	if ($chrflag == 1) {
		$chr =~ s/chr//;
	}
#	my $pos1 = $splitted[1] +1 ;#including this to account for the way junctions are annotated
#	my $pos2 = $splitted[2] -1 ;# ie fusing point A to point B doesn't actually include either of those points.  
	my $pos1 = $splitted[1] ; 
	my $pos2 = $splitted[2] ; 
	my $circlesize = $pos2 - $pos1;
	my $offset = $offsetval; 
	if ($circlesize < (2*$offsetval) ){
		$offset = int($circlesize/2); 
	}  
	my $pos2a = $pos2 - $offset ; 
	my $pos1a = $pos1 + $offset ;
	my $cmdA = "samtools faidx $ref $chr:$pos2a-$pos2";
	my $cmdB = "samtools faidx $ref $chr:$pos1-$pos1a";
	my @fastaA=`$cmdA`;
	my @fastaB=`$cmdB`;
	my $circleseq =""; 
	#contstruct a sequnce that goes: pos2a-pos2pos1-pos1a
	while (@fastaA) {
                my $x = shift(@fastaA);
                chomp $x ;
                if ($x =~ m/^[actgACTG]/) {
                        $sequence .=$x ;
			$circleseq .= $x; 
                }
        }
	while (@fastaB) {
                my $x = shift(@fastaB);
                chomp $x ;
                if ($x =~ m/^[actgACTG]/) {
                        $sequence .=$x ;
			$circleseq .= $x;
                }
        }
	my $exonpos1 = $chrmpos + $offset - $exon_halfsize; 
	my $exonpos2 = $exonpos1 + 2*$exon_halfsize; 
	#print SAF "$y\tcRNA_custom\t$chrmpos\t"; 
	print GTF "cRNA_custom\t$crnasource\tgene\t$chrmpos\t"; 
	$chrmpos += length($circleseq);
	$chrmpos--;  
	#print SAF "$chrmpos\t+\n";
	print GTF "$chrmpos\t.\t+\t.\tgene_id \"$y\";\n"; 
	print GTF "cRNA_custom\t$crnasource\texon\t$exonpos1\t$exonpos2\t.\t+\t.\tgene_id \"$y\"; gene_name \"$y\"; exon_id \"$y.jxn\"; transcript_id \"$y\"; \n";
	$sequence .= $nstring;	
	$chrmpos += length($nstring); 
	$chrmpos++;
}
while (my $chunk = substr($sequence, 0, 80, "")) {
	print "$chunk\n";
}
