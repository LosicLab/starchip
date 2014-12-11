#! /usr/bin/perl
use warnings;
use List::Util qw( min max );
#usage: candidate_stats.pl circs.candidates splicesuffix
 #circs.candidates should be in the format: chr<space>pos1<space>pos2<space>whatever else.  output as circs"${cutoff}"."${minSubjLimit}".investigate from circle_star.sh
 #splicesuffix would be the end of the files we're pulling stats from.  ie if files are form: star_AOR_XXX.5.spliced splicesuffix=.5.spliced

open CANDIDATES, "<$ARGV[0]" or die $!;
$tail=$ARGV[1];
opendir DIRECTORY, "./" or die $!;
@files=readdir(DIRECTORY);
$consensusfilename=$ARGV[0] . $ARGV[1] . ".consensus";
$fullinfofilename=$ARGV[0] . $ARGV[1] . ".allvariants";
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
				#print "$dumbcount\t $hit\n";
				@hitsplit=split(/\s+/, $hit); #17:splicedsize 18:leftborder intron distance 19:rightborder intron distance 20:distance to left enveloping splice 21: distance to right enveloping splice 23:enveloping splice strand 25:#of exons 26:exon start positions(relative, comma sep) 27: exon sizes (comma sep)
				$jxntype = $hitsplit[10];
				$overlapL = $hitsplit[11];
				$overlapR = $hitsplit[12];
				push (@medianAS, $hitsplit[13]);
				push (@splicedsizes, $hitsplit[17]);
				push (@leftIntrons, $hitsplit[18]);
				push (@rightIntrons, $hitsplit[19]);
				push (@leftEnvDist, $hitsplit[20]);
				push (@rightEnvDist, $hitsplit[21]);
				push (@envStrand, $hitsplit[23]);
				push (@envJxn, $hitsplit[24]);
				push (@exons, $hitsplit[25]);
				push (@exonstarts, $hitsplit[26]);
				push (@exonsizes, $hitsplit[27]);
			}
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
	#Spliced Size 
	my ($splicemax, $splicelist) = arraystats(@splicedsizes);
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
        print CONSENSUS "$exonsmax\t";
        print ALLVARIANTS "$exonslist\t";
	#start postitions of exons
	my ($exonstartmax, $exonstartlist) = arraystats(@exonstarts);
        print CONSENSUS "$exonstartmax\t";
        print ALLVARIANTS "$exonstartlist\t";
	#Exon Sizes
	my ($exonsizemax, $exonsizelist) = arraystats(@exonsizes);
        print CONSENSUS "$exonsizemax\n";
        print ALLVARIANTS "$exonsizelist\n";
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
