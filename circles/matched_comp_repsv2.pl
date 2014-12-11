#! /usr/bin/perl
use warnings;

#matched_comp_reps.pl by Kipp Akers
#usage: matched_comp_repsv2.pl infile
#where infile is the output from a bedtools window command similar to: bedtools window 200 mycircs.bed hg19repeats.bed >infile
# 
$alucomp=0;
$alunc=0;
$starter=1;
@AoA = NULL;
$chrom = NULL;
$start = NULL;
$end = NULL;
$samestrand = 0;
$complimentary= 0;
$aluhitindex=0;

open (IN, $ARGV[0]) or die $!;
print "chrom\tstart\tend\tflankingRepeats\tflankingAlus\tcomplimentaryRepeats\tcomplimentaryAlu\tnoncompRepeats\tnoncompAlus\n";
while (my $x = <IN>) {
	@line1 = split(/\s/, $x);
	#0-chrom 1-start 2-end 3-chrom 4-start 5-end 6-repname 7-repscore 8-repstrand 9-14:0s 15-type 16-family
	if ($starter ==1) {
		#initiate
		push @AoA, [ @line1 ];
		$starter++;
		$chrom = $line1[0];
		$start = $line1[1];
		$end = $line1[2];	
		$hitindex=1;
		next;
	}
	#if we have a new cRNA
	if ($line1[0] ne $chrom || $line1[1] != $start || $line1[2] != $end ) {
		#Look through the AoA to calculate stats for the old cRNA
		for my $i (1..($hitindex)) {
			#count flanking repeats that are Alu
			if ($AoA[$i][16] eq "Alu") {
				$aluhitindex++;
			}
			for my $h (($i+1)..$hitindex) {#loop through all combindations of slices from the AoA
				#count complimentary repeats:
				#if same family on complimentary strands at opposite ends	# (repeat2 is above circ     &	repeat1 is below circ ) or       ( rep1 is above cir       & rep2 is below circ)
				if ($AoA[$h][16] eq $AoA[$i][16] && $AoA[$h][8] ne $AoA[$i][8] && (($AoA[$h][4] >=$AoA[$h][2] && $AoA[$i][5] <= $AoA[$i][1]) || ($AoA[$i][4] >=$AoA[$i][2] && $AoA[$h][5] <= $AoA[$h][1]))) {
					$complimentary++; 
					if ($AoA[$h][16] eq "Alu") {
						$alucomp++;
					}
					#can use the below line to look at what non-alu comp repeats we're getting
					#else {print "$AoA[$h][15] $AoA[$h][16]\n";}
				}
				#count non-complimentary repeats
				#if same family on same strand at opposite ends
				if ($AoA[$h][16] eq $AoA[$i][16] && $AoA[$h][8] eq $AoA[$i][8] && (($AoA[$h][4] >=$AoA[$h][2] && $AoA[$i][5] <= $AoA[$i][1]) || ($AoA[$i][4] >=$AoA[$i][2] && $AoA[$h][5] <= $AoA[$h][1]))) {
					$samestrand++;
					if ($AoA[$h][16] eq "Alu") {
						$alunc++;
					}
				}
			}
		}
		print "$chrom\t$start\t$end\t$hitindex\t$aluhitindex\t$complimentary\t$alucomp\t$samestrand\t$alunc\n";
		#clear everything
		$hitindex=1;
		$aluhitindex=0;
		$complimentary=0;
		$alucomp=0;
		$samestrand=0;
		$alunc=0;
		@AoA = NULL;
		push @AoA, [ @line1 ];
		$chrom = $line1[0];
		$start = $line1[1];
		$end = $line1[2];	
	}


	else { #is same candidate circ
		$hitindex++;
		push @AoA, [@line1];
		#check for comp family members
		#print "$x";
	}
}

