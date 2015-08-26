#! /usr/bin/perl
#use warnings;

#filter_circs.pl by Kipp Akers
# This script takes in by standard input the output from chimera-score.pl 
$count=0;
$pluscount=0;
$minuscount=0;
$starter=0;
$matchcount=0;
$supportreads=0;
##$wiggle=5;

$exacto="joinstrands." . $ARGV[0];
open (EXACTMATCHES, ">$exacto") or die $!;
$removeme="removed." . $ARGV[0];
open (REMOVED, ">$removeme") or die $!; 
#print EXACTMATCHES "totalreads\t+strand\t-strand\tchrm\tpos\tstrand\tchrm\tpos\tstrand\tjxntype\toverlapL\toverlapR\n";  
#$details="matchdetails.txt";
#open (MATCHDETAILS, ">>$details") or die $!; 
#$bedfile=$ARGV[0] . ".bed";
#open (BEDFILE, ">$bedfile") or die $!;

while (defined (my $queryline = <STDIN>)) {
	chomp($queryline);
	if ($starter == 0) {
		$starter++;
		@line1=split(/\s+/, $queryline); 
		#0: chromosome of donor 1: first base of the intron of the donor 2: strand of the donor 3: chromosome of the acceptor
	        #4: first base of the intron of the acceptor5: strand of the acceptor 6: junction type: -1=junction is between the mates, 1=GT/AG, 2=CT/AC 7: repeat length to the left of the junction 
        	#8: repeat length to the right of the junction 9: alignscore 10: cigar (both locations separated by a Z)
		$chrm=$line1[0];
		$pos1=$line1[1];
		$pos2=$line1[4];
		$count=1;
		$cigar=$line1[10];
		# i check if the reads fit into the circle they're mapped to: 0=no 1=yes
		$legitimate = &check_badpair($pos1, $pos2, $cigar, $queryline);
		push @legits, $legitimate; 
		if ($line1[2] eq "+") {
                        $pluscount =1;
			push @scores, $line1[9];
			push @plusscores, $line1[9];
                }
                elsif ($line1[2] eq "-") {
                        $minuscount =1;
			push @scores, $line1[9];
			push @minusscores, $line1[9];
                }
		next;
	}
	@line2=@line1;
	#line1 is the current line, line2 is the previous read
	@line1=split(/\s+/, $queryline); # 0chrm 1pos 2strand 3chrm 4pos 5strand 6jxn type 78overlap 9score 10cigar
	#if ($line1[0] eq $chrm && $line1[1]==$pos1 && $line1[4]==$pos2) { #match on chrm and positions of introns
	#print "$line[1] - $pos1 <= $line[7] || $line1[8] &
	if ($line1[0] eq $chrm && ( ($line1[1] - $pos1 <= $line1[7]) || ($line1[1] - $pos1 <= $line1[8]) ) && ( (abs($line1[4] - $pos2) <= $line1[7]) || (abs($line1[4] - $pos2) <= $line1[8] ))) { 
	#match on chrm and positions of introns
		$legitimate=&check_badpair($line1[1], $line1[4], $line1[10], $queryline);
		push @legits, $legitimate; 
		$count++;
		$matchcount++;
		if ($line1[2] eq "+") {
			$pluscount++;
			push @scores, $line1[9];
			push @plusscores, $line1[9];
		}
		elsif ($line1[2] eq "-") {
			$minuscount++;
			push @scores, $line1[9];
			push @minusscores, $line1[9];
		}
		#print MATCHDETAILS "$ARGV[0] Matchcount:$matchcount\n@line2\n@line1\n";
	}
	else {	#ie we've found a new backsplice
		if ($count > $supportreads) { #supportreads is set to zero, so this doesn't really filter.  
			$plusmedian=&median(@plusscores);
			$minusmedian=&median(@minusscores);
			$median=&median(@scores);
			$avglegit=&avg(@legits);
			if ($avglegit >= 0.95) {
				#at this point I flip the positions, so it goes low to high.  
				print EXACTMATCHES "$count\t$pluscount\t$minuscount\t$line2[3]\t$line2[4]\t$line2[5]\t$line2[0]\t$line2[1]\t$line2[2]\t$line2[6]\t$line2[7]\t$line2[8]\t$median\t$plusmedian\t$minusmedian\tchr$line2[3]:$line2[4]-$line2[1]\n";
				#			0     1           2               3       4               5       6       7               8       9               10      11      12        13            14
				# tchr$line2[3]:$line2[4]-$line2[1]\n";
				#  15		
				#if ($pluscount >0) { 	
				#	print BEDFILE "$chrm\t$pos1\t$pos2\tIDplaceholder\t$pluscount\t+\t$pos1\t$pos2\t255,0,0\t1\t1\t0\n"
				#}
				#if ($minuscount >0 ) {
				#	print BEDFILE "$chrm\t$pos2\t$pos1\tIDplaceholder\t$minuscount\t-\t$pos2\t$pos1\t0,0,255\t1\t1\t0\n"
				#}	
			}
			else {
				print REMOVED "$avglegit\t$count\t$pluscount\t$minuscount\t$line2[3]\t$line2[4]\t$line2[5]\t$line2[0]\t$line2[1]\t$line2[2]\t$line2[6]\t$line2[7]\t$line2[8]\t$median\t$plusmedian\t$minusmedian\tchr$line2[3]:$line2[4]-$line2[1]\n";
			}
			$matchcount=0;
			$chrm=$line1[0];
	                $pos1=$line1[1];
        	        $pos2=$line1[4];
                	$count=1;
			@scores =();
			@plusscores=();
			@minusscores=();
			@legits=();
			$legitimate=&check_badpair($pos1, $pos2, $line1[10], $queryline);
			push @legits, $legitimate;
			if ($line1[2] eq "+") {
				$pluscount = 1;
				$minuscount = 0;
				push @scores, $line1[9];
				push @plusscores, $line1[9];
			}
			elsif ($line1[2] eq "-") {
				$pluscount = 0;
				$minuscount = 1;
				push @scores, $line1[9];
				push @minusscores, $line1[9];
	}	}	}
}
$plusmedian=&median(@plusscores);
$minusmedian=&median(@minusscores);
$median=&median(@scores);
$avglegit=&avg(@legits);
if ($avglegit >= 0.95) {
	print EXACTMATCHES "$count\t$pluscount\t$minuscount\t$line2[3]\t$line2[4]\t$line2[5]\t$line2[0]\t$line2[1]\t$line2[2]\t$line2[6]\t$line2[7]\t$line2[8]\t$median\t$plusmedian\t$minusmedian\tchr$line2[3]:$line2[4]-$line2[1]\n";
	#if ($pluscount >0) {
	#	print BEDFILE "$chrm\t$pos1\t$pos2\tIDplaceholder\t$pluscount\t+\t$pos1\t$pos2\t255,0,0\t1\t1\t0\n";
	#}
	#if ($minuscount >0 ) {
	#	print BEDFILE "$chrm\t$pos2\t$pos1\tIDplaceholder\t$minuscount\t-\t$pos2\t$pos1\t0,0,255\t1\t1\t0\n"
	#}
}
else {
	print REMOVED "$avglegit\t$count\t$pluscount\t$minuscount\t$line2[3]\t$line2[4]\t$line2[5]\t$line2[0]\t$line2[1]\t$line2[2]\t$line2[6]\t$line2[7]\t$line2[8]\t$median\t$plusmedian\t$minusmedian\tchr$line2[3]:$line2[4]-$line2[1]\n";
}


sub median {
    my (@data) = sort { $a <=> $b } @_;
    if ( scalar(@data) % 2 ) {
        return ( $data[ @data / 2 ] );
    } else {
        my ( $upper, $lower );
        $lower = $data[ @data / 2 ];
        $upper = $data[ @data / 2 - 1 ];
        return (( $lower + $upper )/2 );
    }
}
sub avg {
	my $elements = 0 ;
	my $counter = 0;
	foreach my $x (@_) {
		$elements++;
		$counter += $x; 
	}
	my $average=$counter/$elements;
	return ($average)
} 
sub check_badpair {
	#take in positions, cigar string, check to makes sure the read is OK
	my ($p1, $p2, $bigcigar, $fulline )=@_;
	my ($cigarette1, $cigarette2)=split(/Z/, $bigcigar);
	#star puts locations into a cigar string and merges pairs together.  the 'p' flag indicate the distance between pairs.  
	#we want to check that p makes sense given we think these are circles.  
	my $mycigar; 
	if ($cigarette1 =~ m/p/) {
		$mycigar = $cigarette1 ;
	}
	elsif ($cigarette2 =~ m/p/) {
		$mycigar = $cigarette2 ; 
	}
	else { return("1"); #do something for single stranded
	}
	my $aligned_total=0;
	#step through the cigar:
	my @split = split(//, $mycigar);
	my $count="";
	my $pairbuffer="";
	my $spliceout="";
	my $aligned_total="";
	foreach my $x (0..$#split) {
		if ($split[$x] =~ m/[\-0-9]/ ) {
                	$count .= $split[$x]; #rejoin numbers
                }
		else { #when we have a complete cigar term
                        if ($split[$x] eq "S") { } #do nothing for softclipping
                        elsif ($split[$x] eq "M") { #on match, add to aligned total 
                        	$aligned_total += $count ; 
                        }
                        #for padding and skipped ref (usually intron) and deletions move the support index without adding support
                        elsif ($split[$x] eq "p") { $pairbuffer += $count ; }
                        elsif ($split[$x] eq "N") { $spliceout += $count ; }
                        elsif ($split[$x] eq "D") { $spliceout += $count ; }
                        elsif ($split[$x] eq "I") { }#do nothing for insertion.  should count negatively, but no easy way to do this. 
                        $count="";
		}
	}
	my $totaldistance = $p1 - $p2 + 5; #the 5 is a little wiggle room.  
	my $alignmentsize = $pairbuffer + $aligned_total + $spliceout;
	my $checkval = $totaldistance - $alignmentsize;
	#print "$mycigar\t$checkval\t$p1\t$p2\t$totaldistance\t$pairbuffer\t$aligned_total\t$spliceout\n";	
	if ($checkval >= 0) { return("1");}
	else { return("0")};   
}
