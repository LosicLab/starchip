#! /usr/bin/perl
use warnings;
use strict;
use Getopt::Std;
use Cwd 'abs_path';

## usage: fusions-from-star.pl  outputname Chimeric.out.junction  

##IMPORTANT NOTE!  this calls 'samtools' and 'bedtools' and 'razers3'  please have these installed and in your path under those aliases.  
	#It shouldn't crash without them, but you will get error messages rather than some of the outputs.
	#$ module load samtools bedtools optitype
	# ^ should do the trick. 
# the goal of this program is to take a file of potential fusions from STAR (Chimeric.out.tab), and determine the read distribution around the fusion.
# to do:    
	# blat output.  blat is REALLY slow, and because of mem loading, probably faster if run in batches.  
		# think the best solution may be to access the webserver?  can't be.
		# run all blat at once.
		# $ blat t=dna q=rna /path/to/database /path/to/query 
		# the other issue is that blat/blast don't lend themselves to computational filtering easily.  Kind of need a human eye to look at the result and see the telltale FP flags.   

#other things to do:
	# a warning if no valid input
	# add uniqueness/read distribution bonus to score

#Very adjustable variables
my $pairedend =1; #1 means paired end data.  any other value means single end. $spancutoff should be 0 if data is single end.   
my $consensus ="TRUE";  # anything but TRUE will make this skip the consensus output. 
my $cutoff = 2; # number of minimum read support at jxn.  reccommend keep this above 2 or your run can take a long time.  
my $cutoff2 = 2; # number of unique read support values (higher indicates more likely to be real. lower is more likely amplification artifact)
		 #this cutoff2 value ranges from 1-70.  Our mapping requires 15bp overhang on a fusion so 84-99 and 0-14 will alwasy be 0 and max.  
		 #may make sense to scale cutoff2 by the read support for a given fusion.  ie if only 15 reads support it, 20 is impossible.  
my $spancutoff = 2; #minimum number of non-split reads support.  
my $wiggle = 500 ; #number of base-pairs of 'wiggle-room' when determining the location of a fusion (for spanning read counts)
my $overlapLimit = 5; #wiggle room for joining very closely called fusion sites. 
my $samechrom_wiggle = 20000; #this is the distance that fusions have to be from each other if on the same chromosome.  Set to 0 if you want no filtering of same-chromosome proximal fusions.
my $lopsidedupper = 5; # (topsidereads + 0.1) / (bottomsidereads + 0.1) must be below this value. set very high to disable.  Reccomended setting 5
my $lopsidedlower = 0.2; # (topsidereads + 0.1) / (bottomsidereads + 0.1) must be above this value. set to 0 to disable. Reccomended setting 0.2
my $refbed="/sc/orga/work/akersn01/ref/Homo_sapiens.GRCh37.74.chr.gtf.bed";
my $repeatbed="/sc/orga/work/akersn01/ref/repeats_hg19plus.bed";
my $ref = "/sc/orga/work/akersn01/ref/star/ENSEMBL.homo_sapiens.release-75_overhang100/Homo_sapiens.GRCh37.75.dna.primary_assembly.fa";

print "Using the following variables:\nPaired-End: $pairedend\nSplit Reads Cutoff: $cutoff\nUnique Support Values Min: $cutoff2\nSpanning Reads Cutoff: $spancutoff\nLocation Wiggle Room (spanning reads): $wiggle bp\nLocation Wiggle Room (split reads) : $overlapLimit bp\nMin-distance : $samechrom_wiggle bp\nRead Distribution upper limit: $lopsidedupper X\nRead Distribution lower limit: $lopsidedlower X\n";

#Scoring Parameters
my $splitscoremod = 10; 
my $spanscoremod = 20;
my $skewpenalty = 4;
my $repeatpenalty = 0.5 ; # score = score*(repeatpenalty^repeats)  --> a fusion can have 0,1,or 2 sites fall into repeat regions. 

#Adjustable Variables (shouldn't change unless the format of Chimeric.out.junction from star changes:
my $col_jxntype=6;
my $col_startposA=10;
my $col_startposB=12;
my $col_cigarA=11;
my $col_cigarB=13;
my $col_chrA=0;
my $col_chrB=3;
my $col_FusionposA=1;
my $col_FusionposB=4;
my $col_strandA=2;
my $col_strandB=5;
my $col_overlapL=7;
my $col_overlapR=8;
my $numbcolumns=14; #need this one in case you junciton.out file/input changes.  this should be a non-existant final column (ie Chimeric.junciton.out has 0-13 columns)
			#could fix this to be a lot more robust...
my $chrflag=0; #1: chromosomes are format: chr1,chr2,chrX  0:chromosomes are format: 1,2,X

#set variables
my $linecount=0;
my $readlength =0;
my %fusions =();
my $script_dir=abs_path($0);
$script_dir =~ s/hashtest.pl//;
my $consensusloc= $script_dir . 'consensus.sh';
my $annotateloc= $script_dir . 'coordinates2genes.sh';
my $troublemakers = $script_dir . 'pseudogenes.txt';


#file management
if (length(@ARGV) != 1) { die "Wrong number of arguments!\n";}
my $outbase = $ARGV[0];
my $outsumm = $outbase . ".summary";
print "your final outputs will be in $outsumm and $outsumm.annotated\n";
my $junction = $ARGV[1];
my $sam = $junction;
$sam =~ s/junction$/sam/ ; 
print "opening $troublemakers\n";
open PSEUDOS, "<$troublemakers" or die $!;
my %pseudogenes;
while (<PSEUDOS>) {
	chomp;
	my ($gene1, $gene2) = split /\s+/; 
	$pseudogenes{$gene1} .= exists $pseudogenes{$gene1} ? ",$gene2" : $gene2; 
}



# primary data format; $chr1_pos1_chr2_pos2_strandA_strandB[0/1/2][0-RL]
	#where strand is + or -
	#and [0/1/2] indicates left-side/right-side of fusion/non-split reads
	#0-Read length starts at the fusion site =0 and expands outwards from there, to pos2+RL on right side, pos1-RL on left side.
	# in most cases, we will use + strand notation for positions.

my $outanno = $ARGV[1];
open ANNOTATED, "$outanno" or die $!;
open (SUMM, ">$outsumm") or die;
print SUMM "Partner1\tPartner2\tScore\tDiscordantReads\tSplitReads\tNearGene1\tDistance1\tNearGene2\tDistance2\n";
while (my $x = <ANNOTATED>) {
        my @line=split(/\s+/, $x);
	my $cols = scalar(@line);
        my @geneannot=grep(/gene_id/, @line);
	my @indices=grep{ $line[$_] =~ /gene_id/ } 0..$#line ; 
	if (@geneannot) {
		my @anno1=split(/;/, $geneannot[0]);
       		my @anno2=split(/;/, $geneannot[1]);
	        my @gene1name=grep(/gene_name/, @anno1);
        	my @gene2name=grep(/gene_name/, @anno2);
	        #skip same-gene 'fusions'
		if ($gene1name[0] eq $gene2name[0]) {
            	    next;  }
		#skip commonly confused gene pairs:
		my $dist2 = $line[($cols-1)];
		my $dist1 = $line[($cols-4)];
		$gene1name[0] =~ s/gene_name:// ;
		$gene1name[0] =~ s/"//g ;  
		$gene2name[0] =~ s/gene_name:// ;
		$gene2name[0] =~ s/"//g ;
		#use a hash to eliminate known troublemaker genes. 
		#if (($gene1name[0] eq "HLA-B" && $gene2name[0] eq "HLA-C") || ($gene1name[0] eq "HLA-C" && $gene2name[0] eq "HLA-B")) { next; }
		#if (($gene1name[0] eq "FTL" && $gene2name[0] eq "FTLP3") || ($gene1name[0] eq "FTLP3" && $gene2name[0] eq "FTL")) { next; }
		#if (($gene1name[0] eq "CES1" && $gene2name[0] eq "CES1P1") || ($gene1name[0] eq "CES1P1" && $gene2name[0] eq "CES1")) { next; }
		#if (($gene1name[0] eq "CES1" && $gene2name[0] eq "CES1P1") || ($gene1name[0] eq "CES1P1" && $gene2name[0] eq "CES1")) { next; }
		no warnings 'uninitialized'; 
		if (($pseudogenes{$gene1name[0]} =~ m/$gene2name[0]/ )|| ($pseudogenes{$gene2name[0]} =~ m/$gene1name[0]/)) { print "skip $gene1name[0] $gene2name[0] \n"; next; }
		else { print "not skipping $gene1name[0] $gene2name[0] \n"; }
		my $score = $line[3]; 
		$score = $score*($repeatpenalty**$line[($indices[0]-2)]) ;   
		print SUMM "$line[0]\t$line[1]\t$score\t$line[3]\t$line[4]\t$gene1name[0]\t$dist1\t$gene2name[0]\t$dist2\n";
	}
}



### BEGIN SUBROUTINES ###
##also some reference info on cigar operators:
##cigar operators:
# M match
# I insertion (into the reference)
# D deletion from ref
# N skipped over ref seq
# S soft clipping
# H hard clipping
# P padding (silent deletion from padded ref)

# below subroutine modified from http://bioinfomative.blogspot.com/2012/07/parsing-cigar-and-md-fields-of-sambam.html
# Given a Cigar string return an estimate of the read length
sub splitCigar {
	no warnings 'uninitialized';
	my $string = $_[0];
	my @split = split(//, $string);
	my $count="";
	my %matches=();
	foreach my $x (0..$#split) {
		#print "$x: $split[$x]\n";
		if ($split[$x] =~ m/[0-9]/ ) {
			#print "number detected\t";
			$count .= $split[$x];
			#print "$count\n";
		}
		else {
			#print "non number $split[$x]\t";
			$matches{$split[$x]} += $count;
			$count=0;
			#print "$split[$x] : $matches{$split[$x]}\n";
		}
	}
	$matches{"S"} + $matches{"M"} + $matches{"I"};
}
sub supportCigar { #input: SJ line 
##this goes through a line and fills in two arrays (which are part of a hash of array of arrays).
# The main goal is to get read support for the readlength around a fusion.  This fills in arrays outside those bounds, but that data is meaningless
# the reason is that only array points defined are global.  those created in the subroutine, stay here, so read support for two fusions at 1000nt outside the fusion will be concatenated. 
# I could fix this by making the above definitions arbitrarily large, but the values still wouldn't hold meaningful information.  
	my $fusionname = $_[$numbcolumns];
##deal with non-split reads, exit
        if ($_[$col_jxntype] eq "-1") {
		if ($_[($numbcolumns+1)] ~~"1") {
			#mates are same orientation as the name
	        	$fusions{$fusionname}[2][0]++;
        		return;
		}
		if ($_[($numbcolumns+1)] ~~"2")	{
			#mates are in the inverse orientation wrt the name
			$fusions{$fusionname}[2][5]++;
			return;
		}
	}
##count read support by side of fusion:
	$fusions{$fusionname}[2][$_[($numbcolumns+1)]]++;
##Left side of the fusion
  ## - strand
	if ($_[$col_strandA] eq "-") {
		my $supportindex = 0;  #support index should start at 0 and move up to read length.  this is the index that traverses the reference seq to add support.
		my $cigarA = $_[$col_cigarA];
		my @split = split(//, $cigarA);
		my $count="";
		EXITER: foreach my $x (0..$#split) { ##separate cigar terms into ind. characters
			if ($split[$x] =~ m/[0-9]/ ) {
				$count .= $split[$x]; #rejoin numbers
			}
			else { #when we have a complete cigar term
				if ($split[$x] eq "S") { } #do nothing for softclipping
				elsif ($split[$x] eq "M") { #on match, add read support for length of cigar M
					foreach my $y ($supportindex .. ($count+$supportindex)) {
						#print "adding 1 to $supportindex + $y which was $fusions{$fusionname}[0][$y]\n";
						$fusions{$fusionname}[0][$y]++; 
					}
					$supportindex = $supportindex + $count ; #move the support index
				}
				#for padding and skipped ref (usually intron) and deletions move the support index without adding support
				#elsif ($split[$x] eq "p") { $supportindex =$supportindex + $count ;}
				elsif ($split[$x] eq "p") { $count=""; $fusions{$fusionname}[2][3]++; last EXITER ;}
				elsif ($split[$x] eq "N") { $supportindex =$supportindex + $count ;}
				elsif ($split[$x] eq "D") { $supportindex =$supportindex + $count ;}
				elsif ($split[$x] eq "I") { }#do nothing for insertion.  should count negatively, but no easy way to do this. 
				$count="";
			}
		}
	}
  ## + strand
        if ($_[$col_strandA] eq "+") {
                my $supportindex = ($_[$col_FusionposA] - $_[$col_startposA]);  #support index should start at read length and move down to zero.  this is the index that traverses the reference seq to add support.
                #print "supportindex: $supportindex\t";
		my $cigarA = $_[$col_cigarA];
                my @split = split(//, $cigarA);
                my $count="";
                my $pskip=0; 
		#pskip helps me deal with split/span support.  when the cigar has a P, it means that we're matching R1 + R2 to reference.  But we're only concerned with the partner on the fusion
		#On left side +/- you want the cigar after/before p.
		#On right side +/- you want the cigar before/after p. 
                if ($cigarA =~ m/p/ ) {
                        $pskip = 1;
                }		
		foreach my $x (0..$#split) { ##separate cigar terms
                        if ($split[$x] =~ m/[0-9,-]/ ) {
                                $count .= $split[$x];
                        }
                        else { #when we have a complete cigar term
                                if ($split[$x] eq "S") { } #do nothing for softclipping
                                elsif ($split[$x] eq "M") { #on match, add read support for length of cigar M
                                        if ($pskip ne "1") {
						foreach my $y (reverse (($supportindex-$count) .. $supportindex)) {
	                                        	$fusions{$fusionname}[0][$y]++;
        	                                }
                			}
		                        $supportindex = $supportindex - $count ; #move the support index
                                }
                                #for padding and skipped ref (usually intron) and deletions move the support index without adding support
                                elsif ($split[$x] eq "p") { $supportindex =$supportindex - $count ; $pskip =0; $fusions{$fusionname}[2][3]++;}
                                elsif ($split[$x] eq "N") { $supportindex =$supportindex - $count ;}
                                elsif ($split[$x] eq "D") { $supportindex =$supportindex - $count ;}
                                elsif ($split[$x] eq "I") { }#do nothing for insertion.  should count negatively, but no easy way to do this. 
                        	#print "$supportindex\t$count\t$split[$x]\n";
				$count="";
                        }
                }
        }
	#print "we have fusion name $fusionname\n";
	#foreach my $f (0..($readlength-1)) {
	#	print "$fusions{$fusionname}[0][$f] ";
	#}
	#print "\n";
##Right side of the fusion
  ## + strand
        if ($_[$col_strandB] eq "+") {
                my $supportindex = 0;  #support index should start at 0 and move up to read length.  this is the index that traverses the reference seq to add support.
                my $cigarB = $_[$col_cigarB];
                my @split = split(//, $cigarB);
                my $count="";
                EXITER: foreach my $x (0..$#split) { ##separate cigar terms
                        if ($split[$x] =~ m/[0-9]/ ) {
                                $count .= $split[$x];
                        }
                        else { #when we have a complete cigar term
                                if ($split[$x] eq "S") { } #do nothing for softclipping
                                elsif ($split[$x] eq "M") { #on match, add read support for length of cigar M
                                        foreach my $y ($supportindex .. ($count+$supportindex)) {
                                                $fusions{$fusionname}[1][$y]++;
                                        }
                                        $supportindex = $supportindex + $count ; #move the support index
                                }
                                #for padding and skipped ref (usually intron) and deletions move the support index without adding support
                                #elsif ($split[$x] eq "p") { $supportindex =$supportindex + $count ;}
                               	elsif ($split[$x] eq "p") { $count=""; $fusions{$fusionname}[2][4]++; last EXITER ;}
				elsif ($split[$x] eq "N") { $supportindex =$supportindex + $count ;}
                                elsif ($split[$x] eq "D") { $supportindex =$supportindex + $count ;}
                                elsif ($split[$x] eq "I") { }#do nothing for insertion.  should count negatively, but no easy way to do this. 
                                $count="";
                        }
                }
        }
  ## - strand
	if ($_[$col_strandB] eq "-") {
                my $supportindex = ($_[$col_FusionposB] - $_[$col_startposB]);  #support index should start at read length and move down to zero.  this is the index that traverses the reference seq to add support.
                my $cigarB = $_[$col_cigarB];
                my @split = split(//, $cigarB);
                my $count="";
		my $pskip=0; 
		if ($cigarB =~ m/p/ ) {
			$pskip = 1;
		}
                foreach my $x (0..$#split) { ##separate cigar terms
                        if ($split[$x] =~ m/[0-9,-]/ ) {
                                $count .= $split[$x];
                        }
                        else { #when we have a complete cigar term
                                if ($split[$x] eq "S") { } #do nothing for softclipping
                                elsif ($split[$x] eq "M") { #on match, add read support for length of cigar M
                                        if ($pskip ne "1" ) {
						foreach my $y (reverse (($supportindex-$count) .. $supportindex)) {
                                	                $fusions{$fusionname}[1][$y]++;
                                       		}
					}
                                        $supportindex = $supportindex - $count ; #move the support index
                                }
                                #for padding and skipped ref (usually intron) and deletions move the support index without adding support
                                elsif ($split[$x] eq "p") { $supportindex =$supportindex - $count ; $pskip = 0;$fusions{$fusionname}[2][4]++;}
                                elsif ($split[$x] eq "N") { $supportindex =$supportindex - $count ;}
                                elsif ($split[$x] eq "D") { $supportindex =$supportindex - $count ;}
                                elsif ($split[$x] eq "I") { }#do nothing for insertion.  should count negatively, but no easy way to do this. 
                                #print "$supportindex\t$count\t$split[$x]\n";
                                $count="";
                        }
                }
        } 
}
sub extractSequence {
	#print SUMM ">$_[0]\n";
	my @splitseq = split(/_/, $_[0]);#0:chrm1 1:pos 2:chr2 3:pos2 4:strandA 5:strandB
###obtain the REFERENCE sequence based on positions.
##Left side sequence
	my $seqpos1; my $seqpos2;
	my $chrA = $splitseq[0]; $chrA =~ s/chr//;
  ##if +
	if ($splitseq[4]~~"+") {
		$seqpos2 = ($splitseq[1]-1);
		$seqpos1 = $seqpos2 - $readlength;
	}
  ##if - 
	if ($splitseq[4]~~"-") {
		$seqpos1 = ($splitseq[1]+1);
		$seqpos2 = $seqpos1 +$readlength;
	}
	my $cmdA = "samtools faidx $ref $chrA:$seqpos1-$seqpos2";
##Right side sequence
	my $seqBpos1 ; my $seqBpos2; 
	my $chrB = $splitseq[2]; $chrB =~ s/chr//;
  ##if + 
	if ($splitseq[5] ~~ "+") {
		$seqBpos1 = ($splitseq[3]+1);
		$seqBpos2 = $seqBpos1 + $readlength ; 
	}
  ##if -
	if ($splitseq[5] ~~ "-") {
		$seqBpos2 = ($splitseq[3]-1);
		$seqBpos1 = $seqBpos2 - $readlength;
	}
	my $cmdB = "samtools faidx $ref $chrB:$seqBpos1-$seqBpos2";
##collect REference fasta with samtools
	my @fastaA=`$cmdA`;
	my @fastaB=`$cmdB`;
	#print "results @fasta\n";
	my $sequenceA;
	my $sequenceB;
  ## FASTA-->string of sequence
	while (@fastaA) {
		my $x = shift(@fastaA);
		chomp $x ;
		if ($x =~ m/^[actgACTG]/) {
			$sequenceA .=$x ;
		}
	}
	while (@fastaB) {
		my $x = shift(@fastaB);
		chomp $x ;
		if ($x =~ m/^[actgACTG]/) {
                        $sequenceB .=$x ;
                }
	}
	my $refseq; 
	if ($splitseq[4] ~~ "-") {
		my @reverseA=&revcompl($sequenceA);
		$refseq = join("", @reverseA);
		print SUMM "@reverseA";
	}
	else { print SUMM "$sequenceA"; $refseq = $sequenceA; }
	if ($splitseq[5] ~~ "-") {
                my @reverseB=&revcompl($sequenceB);
		$refseq = join("", $refseq, @reverseB); 
                print SUMM "@reverseB\t";
        }
	else { print SUMM "$sequenceB\t"; $refseq = join("", $refseq, $sequenceB); }
##obtain the consensus sequence from the reads themselves
	if ($consensus eq "TRUE"){
		my $tempID = int(rand(10000000)); 
		#consensus command is : ~/scripts/fusions/consensus.sh chrom1 pos1 chrom2 pos2 junctionfile samfile fusionID reference_sequence
		my $consensuscmd = "$consensusloc $splitseq[0] $splitseq[1] $splitseq[2] $splitseq[3] $junction $sam $tempID $refseq";
		my @consfastq=`$consensuscmd`; 
		my $consensusSeq;
		while (@consfastq) {
			my $x = shift(@consfastq);
			chomp $x;
			next if ($x =~ m/^@/); #skip readID
			last if ($x =~ m/^\+/); #stop after reads
			$consensusSeq .= $x;
 		}
		$consensusSeq =~ s/^n+//;
		if ($consensusSeq ne ""){
			print SUMM "$consensusSeq\n";		
		}
		else { print SUMM".\n"; }
	}
}
sub reversestrand {
	if ($_[0] ~~ "+"){
		return "-";
	}
	elsif ($_[0] ~~ "-"){
		return "+";
	}
}

#this sub taken from user jkahn on perlmonks: http://www.perlmonks.org/?node_id=197793
sub revcompl { # operates on all elements passed in
  my (@dna) = @_;
  my @done;
  foreach my $segment (@dna) {
    my $revcomp = reverse($segment);
    $revcomp =~ tr/ACGTacgt/TGCAtgca/;
    push @done, $revcomp;
  }
  return @done; # or reverse @done;  # what's best semantics?
}
sub adjustposition { #takes in 0pos, 1strand, 2pos 3strand
	my $position1;
	my $position2;
	if ($_[1] ~~ "+") {
		$position1 = $_[0] -1;
	}
	elsif ($_[1] ~~ "-") {
		$position1 = $_[0] +1;
	}
	if ($_[3] ~~ "+") {
		$position2 = $_[2] +1;
	}
	elsif ($_[3] ~~ "-") {
		$position2 = $_[2] -1;
	}
	return ($position1, $position2); 
}

#credit: http://www.bagley.org/~doug/shootout/  slash http://dada.perl.it/shootout/moments.perl.html
sub kurtosis { 
	my @nums=@_;
	#print "@nums\n"; 
	my $sum = 0;
	foreach (@nums) { $sum += $_ }
	my $n = scalar(@nums);
	my $mean = $sum/$n;
	my $average_deviation = 0;
	my $standard_deviation = 0;
	my $variance = 0;
	my $skew = 0;
	my $kurtosis = 0;
	foreach (@nums) {
   		my $deviation = $_ - $mean;
	 	$average_deviation += abs($deviation);
    		$variance += $deviation**2;
    		$skew += $deviation**3;
    		$kurtosis += $deviation**4;
	}
	$average_deviation /= $n;
	$variance /= ($n - 1);
	$standard_deviation = sqrt($variance);
	#print "variance:$variance ";
	if ($variance) {
		$skew /= ($n * $variance * $standard_deviation);
	    	$kurtosis = $kurtosis/($n * $variance * $variance) - 3.0;
		$kurtosis = substr($kurtosis, 0,5);
		$skew = substr($skew, 0,5);
	}
	#print "kurtosis:$kurtosis\n";
	return ($skew, $kurtosis);
}
sub fusionScore {
	my @params=@_; #0:split reads, 1:topsidesplit, 2:bottomsidesplit 3:spanreads 4;topspan 5;bottomspan 6:skew 7:chr1 8;loc1 9;strand1 10;chr2; 11;loc2; 12;strand2
	my $splitscore;
	# Score from Split Reads
	if ($params[1] >= $params[2]) {
		$splitscore = $params[0]/(($params[1]+1)/($params[2]+1));
	}
	else {$splitscore = $params[0]/(($params[2]+1)/($params[1]+1)); }
	my $spanscore;
	#Score from non-Split Reads
	if ($params[4] >= $params[5]) {
		$spanscore = $params[3]/(($params[4]+1)/($params[5]+1));
	}
	else { $spanscore = $params[3]/(($params[5]+1)/($params[4]+1)); }
	my $basescore = $splitscore/$splitscoremod + $spanscore/$spanscoremod ; 
	#Skew Penalty
	if ($params[6] >=0.25 ) { $basescore = $basescore*$skewpenalty; }	
	#Read-through penalty 
	if ($params[7] eq $params[10] && $params[9] eq $params[12] ) {
		my $dist = abs($params[8] - $params[11]);
		my $penalty = 20000/$dist;
		$basescore = $basescore - $penalty*$basescore ; 
	}
	$basescore = sprintf "%.1f", $basescore; 
	return($basescore); 
}
