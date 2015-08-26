#! /usr/bin/perl
use warnings;
use List::Util qw( min max );
#usage: candidate_stats.pl circs.candidates splicesuffix
 #circs.candidates should be in the format: chr<space>pos1<space>pos2<space>whatever else.  output as circs"${cutoff}"."${minSubjLimit}".investigate from circle_star.sh
 #splicesuffix would be the end of the files we're pulling stats from.  ie if files are form: star_AOR_XXX.5.spliced splicesuffix=.5.spliced

#known issue: in rare cases where there is equal support for multiple values, the consensus values don't correlate to eache.
#	ie: 3 subjects with 10 exons and 3 subjects with 11 exons.  Consensus might list 10 exon start values and 11 exon length values.  


open CANDIDATES, "<$ARGV[0]" or die $!;
$tail=$ARGV[1];
opendir DIRECTORY, "./" or die $!;
@files=readdir(DIRECTORY);
$consensusfilename=$ARGV[0] . ".consensus";
$fullinfofilename=$ARGV[0] . ".allvariants";
open CONSENSUS, ">$consensusfilename" or die $!;##             
print CONSENSUS "chrm\tpos1\tpos2\tjxntype\toverlapL\toverlapR\tMedianAS\tGenomicSize\tSplicedSize\tLeftBorder\tRightBorder\tLeftEnvelope\tRightEnvelope\tEnvelopingStrand\tEnvelopingJxnType\tExons\tExonStarts\tExonsizes\n";
open ALLVARIANTS, ">$fullinfofilename" or die $!;
print ALLVARIANTS "chrm\tpos1\tpos2\tjxntype\toverlapL\toverlapR\tMedianAS\tGenomicSize\tSplicedSize\tLeftBorder\tRightBorder\tLeftEnvelope\tRightEnvelope\tEnvelopingStrand\tEnvelopingJxnType\tExons\tExonStarts\tExonsizes\n";

while (my $cRNA=<CANDIDATES>) {
	#initialize some values:
	$jxntype = NULL;
	$overlapL=NULL;
	$overlapR=NULL;
	@medianAS=();
        @splicedsizes=();
        @leftIntrons=();
        @rightIntrons=();
        @leftEnvDist=();
        @rightEnvDist=();
        @envStrand=();
        @envJxn=();
        @exons=();
        @exonstarts=();
        @exonsizes=();
	@spliceinfo=();
	
	#$dumbount=0;
	@cRNAarray=split(/\s+/, $cRNA); #0:chrm 1:pos 2:pos 3:#of inds
	#create the perl style grep query 
	$query=$cRNAarray[0] . '\t' . $cRNAarray[1] . '\t.\t' . $cRNAarray[0] . '\t' . $cRNAarray[2] ;	
	#read all files in cwd
	for my $f (0..$#files) {
		#if the file contains the supplied tail
		if ($files[$f] =~ m/$tail/ ) {
			my $hit=`grep -P \"$query\" $files[$f]`;
			if (length ($hit)) { 
				#$dumbcount++;
				#print "$hit\n";
				@hitsplit=split(/\s+/, $hit); #17:splicedsize 18:leftborder intron distance 19:rightborder intron distance 20:distance to left enveloping splice 21: distance to right enveloping splice 23:enveloping splice strand 25:#of exons 26:exon start positions(relative, comma sep) 27: exon sizes (comma sep)
				$jxntype = $hitsplit[10];
				$overlapL = $hitsplit[11];
				$overlapR = $hitsplit[12];
				push (@spliceinfo, "$hitsplit[18]_$hitsplit[26]_$hitsplit[27]_$hitsplit[28]");
				push (@medianAS, $hitsplit[13]);
				push (@splicedsizes, $hitsplit[18]);
				push (@leftIntrons, $hitsplit[19]);
				push (@rightIntrons, $hitsplit[20]);
				push (@leftEnvDist, $hitsplit[21]);
				push (@rightEnvDist, $hitsplit[22]);
				push (@envStrand, $hitsplit[24]);
				push (@envJxn, $hitsplit[25]);
				push (@exons, $hitsplit[26]);
				push (@exonstarts, $hitsplit[27]);
				push (@exonsizes, $hitsplit[28]);
			}
		#print "$files[$f]\n";
		}
	}
	#print out the 'steady' info about each circ candidate
	#chrom start end jxntype overlapL overlapR 
	print CONSENSUS "$cRNAarray[0]\t$cRNAarray[2]\t$cRNAarray[1]\t$jxntype\t$overlapL\t$overlapR\t"; 
	print ALLVARIANTS "$cRNAarray[0]\t$cRNAarray[2]\t$cRNAarray[1]\t$jxntype\t$overlapL\t$overlapR\t"; #$hitsplit[16]\t";
	#Median alignment score and Genomic Size:
	my $median=&median(@medianAS);
	print CONSENSUS "$median\t$hitsplit[16]\t";
	print ALLVARIANTS "$median\t$hitsplit[16]\t";
	#find the index of the most common splice, using the @spliceinfo array. 
	my %count;
	$count{$_}++ for @spliceinfo;
	# presuming at least one item in @spliceinfo:
	my ($winner, $winner_count) = each %count;
	my @winners=$winner;
	while (my ($maybe, $maybe_count) = each %count) {
        	if ($maybe_count == $winner_count) {
        	        push (@winners, $maybe);
        	}
		if ($maybe_count > $winner_count) {
                	$winner = $maybe;
                	$winner_count = $maybe_count;
                	@winners=$maybe;
        	}		
	}
	#now @winners has the most common splice motifs.  if there are many, select one:
	# select the best one (has the smallest splice size)
	my $topsplice = $winners[0];
	my $smallest_size = $topsplice;
	$smallest_size =~ s/_.+//;
	if (scalar(@winners) > 1) {
        	foreach my $x (@winners) {
                	my $size = $x;
                	$size =~ s/_.+//;
                	if ($size < $smallest_size) {
                        	$topsplice = $x;
                        	$smallest_size = $size;
                	}
        	}
	}	
	#get index of that. 
	my( $index )= grep { $spliceinfo[$_] eq $topsplice } 0..$#spliceinfo;
	## below, I no longer use the arraystats max values, i take the one from $index
	#Spliced Size 
	my ($splicemax, $splicelist) = arraystats(@splicedsizes);
	$splicemax=$splicedsizes[$index];
	print CONSENSUS "$splicemax\t";
	print ALLVARIANTS "$splicelist\t";
	#LeftBorder Introns
	my ($leftmax, $leftlist) = arraystats(@leftIntrons);
	print CONSENSUS "$leftmax\t";
        print ALLVARIANTS "$leftlist\t";
	#RightBorder Introns
        my ($rightmax, $rightlist) = arraystats(@rightIntrons);
        print CONSENSUS "$rightmax\t";
        print ALLVARIANTS "$rightlist\t";
	#Distance to leftside exon of an enveloping splice
        my ($leftEmax, $leftElist) = arraystats(@leftEnvDist);
        print CONSENSUS "$leftEmax\t";
        print ALLVARIANTS "$leftElist\t";
        #Distance to rightside exon of an enveloping splice
        my ($rightEmax, $rightElist) = arraystats(@rightEnvDist);
        print CONSENSUS "$rightEmax\t";
        print ALLVARIANTS "$rightElist\t";
        #Strand of enveloping splice
        my ($Estrandmax, $Estrandlist) = arraystats(@envStrand);
	print CONSENSUS "$Estrandmax\t";
        print ALLVARIANTS "$Estrandlist\t";
	#Enveloping Splice Junction type
	my ($jxnmax, $jxnlist) = arraystats(@envJxn);
	print CONSENSUS "$jxnmax\t";
        print ALLVARIANTS "$jxnlist\t";
	# number of Exons
	my ($exonsmax, $exonslist) = arraystats(@exons);
	$exonsmax = $exons[$index];
        print CONSENSUS "$exonsmax\t";
        print ALLVARIANTS "$exonslist\t";
	#start postitions of exons
	my ($exonstartmax, $exonstartlist) = arraystats(@exonstarts);
        $exonstartmax = $exonstarts[$index];
	print CONSENSUS "$exonstartmax\t";
        print ALLVARIANTS "$exonstartlist\t";
	#Exon Sizes
	my ($exonsizemax, $exonsizelist) = arraystats(@exonsizes);
	$exonsizemax = $exonsizes[$index];
        print CONSENSUS "$exonsizemax\n";
        print ALLVARIANTS "$exonsizelist\n";
	#die;
	my $splicesize2;
	my @splices = split(/,/,$exonsizemax); 
	for (@splices) {
		$splicesize2 += $_;
	}
	if ($splicesize2 != $splicemax ){
		print "warning: for cRNA:$cRNAarray[0]\t$cRNAarray[2]\t$cRNAarray[1] there was a splicing error. $splicesize2 != $splicemax. $hitsplit[28] @hitsplit \n";
	}  
	#die;
}
sub arraystats
{
	my %hash = ();
	for my $x (0..$#_) {
		if (defined $_[$x]) {
                        #print "$x $splicedsizes[$x]\n";
                        if (exists $hash{$_[$x]}) {
                                $hash{$_[$x]}++;
                        }
                        else {
                                $hash{$_[$x]} = 1;
                        }
                }
        }
	my $count = 0;
	my $list = ();
	my $max = ();
        foreach my $name (sort { $hash{$b} <=> $hash{$a} } keys %hash) {
                if ($count==0) { $max=$name; }
                $count++;
                $list .= $name . "(" . $hash{$name} . ");";

        }
        return ($max, $list);
}
sub median {
    my @data = ();
	for my $x (0..$#_) {
		if (defined $_[$x]) {
			push (@data, $_[$x]);
		}
	}
    @data = sort { $a <=> $b } @data;
    if ( scalar(@data) % 2 ) {
        return ( $data[ @data / 2 ] );
    } else {
        my ( $upper, $lower );
        $lower = $data[ @data / 2 ];
        $upper = $data[ @data / 2 - 1 ];
        return (( $lower + $upper )/2 );
    }
}
