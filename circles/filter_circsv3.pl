#! /usr/bin/perl
#use warnings;

#filter_circs.pl by Kipp Akers
# 
$count=0;
$pluscount=0;
$minuscount=0;
$starter=0;
$matchcount=0;
$supportreads=0;
##$wiggle=5;

$exacto="joinstrands." . $ARGV[0];
open (EXACTMATCHES, ">$exacto") or die $!;
#print EXACTMATCHES "totalreads\t+strand\t-strand\tchrm\tpos\tstrand\tchrm\tpos\tstrand\tjxntype\toverlapL\toverlapR\n";  
#$details="matchdetails.txt";
#open (MATCHDETAILS, ">>$details") or die $!; 
#$bedfile=$ARGV[0] . ".bed";
#open (BEDFILE, ">$bedfile") or die $!;

while (defined (my $queryline = <STDIN>)) {
	if ($starter == 0) {
		$starter++;
		@line1=split(/\s+/, $queryline);
		$chrm=$line1[0];
		$pos1=$line1[1];
		$pos2=$line1[4];
		$count=1;
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
	@line1=split(/\s+/, $queryline); # 0chrm 1pos 2strand 3chrm 4pos 5strand 6jxn type 78overlap 9score
	#if ($line1[0] eq $chrm && $line1[1]==$pos1 && $line1[4]==$pos2) { #match on chrm and positions of introns
	#print "$line[1] - $pos1 <= $line[7] || $line1[8] &
	if ($line1[0] eq $chrm && ( ($line1[1] - $pos1 <= $line1[7]) || ($line1[1] - $pos1 <= $line1[8]) ) && ( (abs($line1[4] - $pos2) <= $line1[7]) || (abs($line1[4] - $pos2) <= $line1[8] ))) { #match on chrm and positions of introns
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
	else {
		if ($count > $supportreads) {
			$plusmedian=&median(@plusscores);
			$minusmedian=&median(@minusscores);
			$median=&median(@scores);
			#at this point I flip the positions, so it goes low to high.  
			print EXACTMATCHES "$count\t$pluscount\t$minuscount\t$line2[3]\t$line2[4]\t$line2[5]\t$line2[0]\t$line2[1]\t$line2[2]\t$line2[6]\t$line2[7]\t$line2[8]\t$median\t$plusmedian\t$minusmedian\tchr$line2[3]:$line2[4]-$line2[1]\n";
			#if ($pluscount >0) {
			#	print BEDFILE "$chrm\t$pos1\t$pos2\tIDplaceholder\t$pluscount\t+\t$pos1\t$pos2\t255,0,0\t1\t1\t0\n"
			#}
			#if ($minuscount >0 ) {
			#	print BEDFILE "$chrm\t$pos2\t$pos1\tIDplaceholder\t$minuscount\t-\t$pos2\t$pos1\t0,0,255\t1\t1\t0\n"
			#}	
			$matchcount=0;
			$chrm=$line1[0];
	                $pos1=$line1[1];
        	        $pos2=$line1[4];
                	$count=1;
			@scores =();
			@plusscores=();
			@minusscores=();
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
print EXACTMATCHES "$count\t$pluscount\t$minuscount\t$line2[3]\t$line2[4]\t$line2[5]\t$line2[0]\t$line2[1]\t$line2[2]\t$line2[6]\t$line2[7]\t$line2[8]\t$median\t$plusmedian\t$minusmedian\tchr$line2[3]:$line2[4]-$line2[1]\n";
#if ($pluscount >0) {
#	print BEDFILE "$chrm\t$pos1\t$pos2\tIDplaceholder\t$pluscount\t+\t$pos1\t$pos2\t255,0,0\t1\t1\t0\n";
#}
#if ($minuscount >0 ) {
#	print BEDFILE "$chrm\t$pos2\t$pos1\tIDplaceholder\t$minuscount\t-\t$pos2\t$pos1\t0,0,255\t1\t1\t0\n"
#}


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
