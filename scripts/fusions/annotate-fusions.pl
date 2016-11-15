#! /usr/bin/perl
use warnings;
use strict;
use Cwd 'abs_path';
use Getopt::Long;
use File::Basename;
## usage: fusions-from-star.pl  outputname Chimeric.out.junction  

if (scalar(@ARGV) != 3 ) { die "Wrong number of inputs. Usage: fusions-from-star.pl output_seed input_chimeric.out.junction params.txt \n Be sure you have samtools, bedtools, and mafft available.\n";}

##Read in User Parameters (taken from Perl Cookbook "8.16. Reading Configuration Files")
my %Configs = ();
open CONFIG, "<$ARGV[2]" or die $!;
while (<CONFIG>) {
    chomp;                  # no newline
    s/#.*//;                # no comments
    s/^\s+//;               # no leading white
    s/\s+$//;               # no trailing white
    next unless length;     # anything left?
    my ($var, $value) = split(/\s*=\s*/, $_, 2);
    $Configs{$var} = $value;
}

my $numbcolumns=14; #need this one in case you junciton.out file/input changes.  this should be a non-existant final column (ie Chimeric.junciton.out has 0-13 columns)
#should fix this to be more robust...

 
#set variables
my $script_dir=abs_path($0);
$script_dir =~ s/annotate-fusions.pl//;
my $consensusloc= $script_dir . 'consensus.sh';
my $annotateloc= $script_dir . 'coordinates2genes.sh';
my $blastscript = $script_dir . 'check-pseudogenes.sh'; 
my $smithwaterman = $script_dir . 'smith_waterman.pl';
my $sw_matrix = $script_dir . 'data/' . 'sw_scoremat.txt';
my $readlength ; 
my $col_chrA=0; 
my $col_chrB=3;
my $col_cigarA=11;
my $col_cigarB=13;

#file management
my $outbase = $ARGV[0];
my $outsumm = $outbase . ".summary";
my $outsummtemp = $outsumm . ".temp";
my $outannotemp = $outsummtemp . ".annotated" ;
my $outanno = $outsumm . ".annotated";
my $junction = $ARGV[1];
my $sam = $junction;
$sam =~ s/junction$/sam/ ;
my ($trash, $stardir) = fileparse($junction) ;   
my $data_dir = $script_dir ;
$data_dir =~ s/scripts\/fusions//;


my $troublemakers = $data_dir .  $Configs{falsepositives} ;
my $familyfile = $data_dir . $Configs{familyfile} ;
my $cnvfile = $data_dir . $Configs{cnvs} ;

unless (-e $troublemakers ) { #if the file isn't in starchip/
	#print "$troublemakers\n"; 
        $troublemakers = $Configs{falsepositives} ; #check the absolute path
        unless (-e $troublemakers ) {
                print "Warning: Can not find your False Positives File: $Configs{falsepositives}\n";
        }
}
unless (-e $familyfile ) { #if the file isn't in starchip/
        $familyfile = $Configs{familyfile} ; #check the absolute path
        unless (-e $familyfile ) {
                print "Warning: Can not find your Gene Families File: $Configs{familyfile}\n";
        }
}
unless (-e $cnvfile ) { #if the file isn't in starchip/
        $cnvfile = $Configs{cnvs} ; #check the absolute path
        unless (-e $cnvfile ) {
                print "Warning: Can not find your Copy Number Variants File: $Configs{cnvs}\n";
        }
}

#create a hash of ensembl family members, these are common FP
my %famhash;
if (-e $familyfile) {
	open FAMILY, "<$familyfile" or die "opening $familyfile: $!"; 
	while (<FAMILY>) {
		chomp;
		my ($family, $gene) = split /\s+/;
		if (exists $famhash{$gene}) { $famhash{$gene} .= ",$family"; }
		else { $famhash{$gene} = $family; }
	}
}
else { print "Gene Family file $familyfile not found.  Will continue without it\n"; } 

#create a hash of known pseudogenes---loci that consistantly show FP fusions. 
my %trouble_spots;
if (-e $troublemakers) { 
	open TROUBLE, "<$troublemakers" or die "opening $troublemakers: $!";
	while (<TROUBLE>) {
	        chomp;
	        my ($gene1, $gene2) = split /\s+/;
	        if (exists $trouble_spots{$gene1}) { $trouble_spots{$gene1} .= ",$gene2";}
		else { $trouble_spots{$gene1} = "$gene2";}
	}
}
else { print "Known pseudogenes/paralogs/False-Positives file $troublemakers not found.  Will continue without it\n"; }
  
#CNVs
my @cnvlocs;
my $cnvindex = 0;
if (-e $cnvfile) {
	open CNVS, "<$cnvfile" or die $!; 
	while (my $x = <CNVS>) {
		$x =~ s/^chr// ; #strip out 'chr' if it exists (makes comparisons more robust) 
		my @line = split(/\s+/, $x);
		#cnv lines are formatted: chrm pos1 pos2
        	$cnvlocs[$cnvindex][0]{$line[0]} = $line[1]; ##data format: $AbParts[index][0=lower/1=upper position]{chrom}
        	$cnvlocs[$cnvindex][1]{$line[0]} = $line[2]; ##data format: $AbParts[index][0=lower/1=upper position]{chrom}
        	$cnvindex++;
	}
}
else { print "Copy number variants file $cnvfile not found.  Will Continue without it\n"; }

#Get read length
open JUNCTION, "<$junction" or die $!;
while (my $x = <JUNCTION>) {
	my @line = split(/\s+/, $x);
	my $cigarA=$line[$col_cigarA];
        my $cigarB=$line[$col_cigarB];
        my $lengthA = &splitCigar($cigarA);
        my $lengthB = &splitCigar($cigarB);
        if ($lengthA == $lengthB) {
                $readlength=$lengthA;
                print "Read length appears to be $readlength\n";
        }	
        elsif ($lengthA == (2*$lengthB)) {
                $readlength=$lengthB;
                print "Read length appears to be $readlength\n";
        }
        elsif ($lengthB == (2*$lengthA)) {
                $readlength=$lengthA;
                print "Read length appears to be $readlength\n";
        }
        else { print "read length error, please check input.  Using 100bp as a guess.\n";
                $readlength=100;
        }
	last ; #exit after the first line. 
}
close (JUNCTION) ;

##Post filter gene-annotation here:
#Run an External Script to Attach Gene Information to my Fusions (using bedtools intersect)
my $annotate_command = $annotateloc . " " . $outsummtemp . " " . $Configs{refbed} . " " . $Configs{repeatbed}; 
#print "running $annotate_command\n"; 
system($annotate_command); 
#print "now annotated\n"; 
#should create a file $outsumtemp.annotated with gene annotations.

#remove same gene fusions, simplify output
open ANNOTEMP, "<$outannotemp" or die $!;
open SUMM, ">$outsumm" or die $!;
open ANNOTATION, ">$outanno" or die $!; 
print SUMM "Partner1\tPartner2\tScore\tDiscordantReads\tSplitReads\tAvgAS\tNearGene1\tDistance1\tNearGene2\tDistance2\tConsensusSeq\n";
print ANNOTATION "Partner1\tPartner2\tScore\tSpanningReads\tSplitReads\tTopsideCrossing\tBottomsideCrossing\tChromAAnchors\tChromBAnchors\tUniqueSupportLeft\tUniqueSupportRight\tKurtosis\tSkew\tLeftAnchor\tRightAnchor\tTopsideSpanning\tBottomsideSpanning\tRepeats\tID1\tGeneInfo1\tDistance1\tID2\tGeneInfo2\tDistance2\tReferenceSequence\tConsensusSequence\tAvgAlignScore\n";
my $maxAS= $readlength/2;
my $minAS= 15; #this 15 is a changeable star parameter (--chimSegmentMin 15) 
my $avgAStarget= ($maxAS+$minAS)/2;  #a great fusion would have this average overhang.  (~32 for 100bp reads)
EXIT_ANNO_FILTER: while (my $x = <ANNOTEMP>) {
	chomp $x; 
        my @line = split(/\s+/, $x);
	my $cols = scalar(@line);
        my @geneannot = grep(/gene_id/, @line);
	#indices stores the positions of gene_id in .annotated.  # repeats is 1 less (than indices[0], #alignments score is 3 less.  For perl indexes, these become -2 and -4.
	my @indices=grep{ $line[$_] =~ /gene_id/ } 0..$#line ; 
	if (@geneannot) {
		my @anno1=split(/;/, $geneannot[0]);
       		my @anno2=split(/;/, $geneannot[1]);
		my @gene1name=grep(/gene_name/, @anno1);
        	my @gene2name=grep(/gene_name/, @anno2);
		#skip same-gene 'fusions'
		if ($gene1name[0] eq $gene2name[0]) {
			print "Skipping same-gene 'fusion' within $gene1name[0]\n"; 
        	    	next;  
		}
		#skip CNVs and circular-looking
		my ($chrom1, $position1, $strand1 ) = split ':', $line[0] ; $chrom1 =~ s/:.*//;  
		my ($chrom2, $position2, $strand2 ) = split ':', $line[1] ; $chrom2 =~ s/:.*//;  
		if ($chrom1 eq $chrom2) { #same chrom
			if ($strand1 eq $strand2) { #same strand
				#remove circular looking RNA.  If not, widen the fusion positions for the CNV hunt.
				if ($strand1 eq "-") { 
					if (($position2 - $position1) < $Configs{circlesize} && ($position2 - $position1) >0) { #ie a back splice
						print "skipping circle $line[0] $line[1]\n";
						next EXIT_ANNO_FILTER ;
					}#if not circular, adjust 
					my $temp = $position1 - $Configs{cnvwiggle} ; $position1 = $position2 + $Configs{cnvwiggle} ; $position2 = $temp;
				}
				else   { ## ie plus strand
					if (($position1 - $position2) < $Configs{circlesize} && ($position1 - $position2) > 0) {
						print "skipping circle $line[0] $line[1]\n";
						next EXIT_ANNO_FILTER ;
					} #if not circular, adjust 
					$position1 = $position1 + $Configs{cnvwiggle}; $position2 = $position2 - $Configs{cnvwiggle};
				}
				#check the positions for CNV status
				$chrom1 =~ s/^chr// ; ##CNVs chr values are stripped on read-in.  
				foreach my $z (0..$cnvindex) { 
	                 		no warnings 'uninitialized' ;# if positions 1 and 2 are within the CNV. 
		                	if ($cnvlocs[$z][0]{$chrom1} <= $position1 && $position1 <= $cnvlocs[$z][1]{$line[$col_chrA]}) {
                                		if ($cnvlocs[$z][0]{$chrom1} <= $position2 && $position2 <= $cnvlocs[$z][1]{$chrom1}) {
		                                        print "skipping CNV $line[0] $line[1]\n";
							next EXIT_ANNO_FILTER;
	                		        }
		                	}
                		}
			}
		}
		#define gene names/distances
		my $dist2 = $line[($cols-1)]; my $dist1 = $line[($cols-4)];
		$gene1name[0] =~ s/gene_name:// ; $gene1name[0] =~ s/"//g ;  
		$gene2name[0] =~ s/gene_name:// ; $gene2name[0] =~ s/"//g ;
		# check that the partners aren't family members:
		my @intersection;
		no warnings 'uninitialized'; 
		my @gene1fams = split (/,/, $famhash{$gene1name[0]}); my @gene2fams = split (/,/, $famhash{$gene2name[0]});
		my %count = ();
	        foreach my $element (@gene1fams, @gene2fams) { $count{$element}++ }
		        foreach my $element (keys %count) {
       			if ($count{$element} > 1) {
				print "skipping known family members $gene1name[0] and $gene2name[0] in family $element with $count{$element} count\n"; next EXIT_ANNO_FILTER;
			}
		}
		#use a hash to eliminate known pseudogene 'fusions' 
		no warnings 'uninitialized';
                if (($trouble_spots{$gene1name[0]} =~ m/$gene2name[0]/ )|| ($trouble_spots{$gene2name[0]} =~ m/$gene1name[0]/)) { print "skipping known false positive pair $gene1name[0] and $gene2name[0]\n" ; next; }
		#pull out consensus sequence and alignment score
		#calculate fusion score
		my ($refseq, $consSeq, $alignScore) = &extractSequence($line[0], $line[1]);
		my $score = $line[2];
		$score = $score*($Configs{repeatpenalty}**$line[($indices[0]-2)]) ;#penalize repeats : $line[($indices[0]-2)] is the # of repeats in this fusion
		#we want to penalize fusions with overhangs less than avgAStarget.  
		my $OverhangScoreMod;
		my $ASdifference = $alignScore - $minAS; 
		if ( $alignScore - $avgAStarget > 0) { #ie if our overhang is above where we expect it.
			$OverhangScoreMod = 1 ; 
		}
		else { # if our alignment score is above the min, below expected average
			my $ASrange = $avgAStarget - $minAS ;
			$OverhangScoreMod = $ASdifference/$ASrange ; 
		}
		$score = $score*($OverhangScoreMod);  #fairly heuristic.  Example : average overhang = 20 and read length 100: ScoreMod = 0.286 = (20-15)/(32.5-15)      
		print SUMM "$line[0]\t$line[1]\t$score\t$line[3]\t$line[4]\t$alignScore\t$gene1name[0]\t$dist1\t$gene2name[0]\t$dist2\t$consSeq";
		print ANNOTATION "$x\t$refseq\t$consSeq\t$alignScore";
		#get consensus sequence
		if ($Configs{simscore} eq 'TRUE') {
			&SIMSCORE($line[0], $line[1]);
		}
	print SUMM "\n";
	print ANNOTATION "\n"; 
	}
}
print "For a circos plot, run: Rscript ${script_dir}circos.plot.R $Configs{refbed} $outsumm myCircosPlotName \n"; 


### BEGIN SUBROUTINES ###
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

sub extractSequence {
	#input is format chr:pos:+ chr:pos:-
	my ($chrA, $posA, $strandA) = split(/:/, $_[0]); # (/__/, $_[0]);#0:chrm1 1:pos 2:chr2 3:pos2 4:strandA 5:strandB
	my ($chrB, $posB, $strandB) = split(/:/, $_[1]);
	#print "chrA $chrA posA $posA strandA $strandA chrb $chrB posb $posB strandb $strandB \n";
###obtain the REFERENCE sequence based on positions.
##Left side sequence
	my $seqpos1; my $seqpos2;
  ##if +
	if ($strandA eq "+") {
		$seqpos2 = $posA ;
		$seqpos1 = $seqpos2 - $readlength;
	}
  ##if - 
	if ($strandA eq "-") {
		$seqpos1 = $posA;
		$seqpos2 = $seqpos1 +$readlength;
	}
	my $cmdA = "samtools faidx $Configs{refFasta} '$chrA:$seqpos1-$seqpos2'";
##Right side sequence
	my $seqBpos1 ; my $seqBpos2; 
  ##if + 
	if ($strandB eq "+") {
		$seqBpos1 = $posB;
		$seqBpos2 = $seqBpos1 + $readlength ; 
	}
  ##if -
	if ($strandB eq "-") {
		$seqBpos2 = $posB;
		$seqBpos1 = $seqBpos2 - $readlength;
	}
	my $cmdB = "samtools faidx $Configs{refFasta} '$chrB:$seqBpos1-$seqBpos2'";
	#print "$cmdA\n$cmdB\n";
##Collect Reference fasta with samtools
	my @fastaA=`$cmdA`;
	my @fastaB=`$cmdB`;
	my $sequenceA;
	my $sequenceB;
  ## FASTA-->string of sequence
	while (@fastaA) {
		my $x = shift(@fastaA);
		chomp $x ;
		if ($x =~ m/^[actgnACTGN]/) {
			$sequenceA .=$x ;
		}
	}
	while (@fastaB) {
		my $x = shift(@fastaB);
		chomp $x ;
		if ($x =~ m/^[actgnACTGN]/) {
                        $sequenceB .=$x ;
                }
	}
	my $refseq; 
	if ($strandA eq "-") {
		my @reverseA=&revcompl($sequenceA);
		$refseq = join("", @reverseA);
	}
	else { 
		$refseq = $sequenceA; 
	}
	if ($strandB eq "-") {
                my @reverseB=&revcompl($sequenceB);
		$refseq = join("", $refseq, @reverseB); 
        }
	else {
		$refseq = join("", $refseq, $sequenceB); 
	}
##obtain the consensus sequence from the reads themselves
	if ($Configs{consensus} eq "TRUE"){
		my $tempID = int(rand(10000000)); 
		my ($unadjposA, $unadjposB) = &unadjustposition($posA, $strandA, $posB, $strandB); 
		#consensus command is : consensus.sh chrom1 pos1 chrom2 pos2 junctionfile samfile fusionID reference_sequence
		my $consensuscmd = "$consensusloc '$chrA' $unadjposA '$chrB' $unadjposB $junction $sam $tempID $refseq $script_dir";
		#print "$consensuscmd\n";
		my @consResults=`$consensuscmd`; 
		my $consensusSeq=$consResults[0];
		chomp $consensusSeq;
		my $avgAS=$consResults[1];
		chomp $avgAS;
		$avgAS = sprintf "%.1f", $avgAS ; 
		if ($consensusSeq ne ""){
			return ($refseq, $consensusSeq, $avgAS);		
		}
		else { 
			return ($refseq, ".", "0"); 
		}
	}
}
sub SIMSCORE { #input is format chr:pos:+ chr:pos:-
        my ($chrA, $posA, $strandA) = split(/:/, $_[0]); 
        my ($chrB, $posB, $strandB) = split(/:/, $_[1]); 
	$chrA =~ s/^chr//;
	$chrB =~ s/^chr//;
	my $seqpos1; my $seqpos2;
	my $seqBpos1; my $seqBpos2;
#check fusion site sequence similarity
	my $tempID = int(rand(10000000));
	my $seqfasta = $tempID . ".fa";
	open (FASTA, ">$seqfasta") or die $!; 
	#This time I want the 10bp on either side of the fusion site from both sites.
        ##if +
        if ($strandA eq "+") {
                $seqpos2 = ($posA+10);
                $seqpos1 = $seqpos2 - 20;
        }
        ##if - 
        if ($strandA eq "-") {
                $seqpos1 = ($posA-10);
                $seqpos2 = $seqpos1 +20;
        }
	my $cmdA = "samtools faidx $Configs{refFasta} '$chrA:$seqpos1-$seqpos2'";
     ##right side
        ##if + 
        if ($strandB eq "+") {
                $seqBpos1 = ($posB-10);
                $seqBpos2 = $seqBpos1 + 20 ;
        }
        ##if -
        if ($strandB eq "-") {
                $seqBpos2 = ($posB+10);
                $seqBpos1 = $seqBpos2 - 20;
        }
        my $cmdB = "samtools faidx $Configs{refFasta} '$chrB:$seqBpos1-$seqBpos2'";
	#run the commands
	#print "$cmdA\n$cmdB\n";
	my @fastaA=`$cmdA`;
        my @fastaB=`$cmdB`;
	#if neg strand, reverse comp.  otherwise print to fa file. 
	my $sequenceA = "";
	my $sequenceB = "";
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
        if ($strandA eq "-") {
                my @reverseA=&revcompl($sequenceA);
		print FASTA ">seq1\n";
		print FASTA "@reverseA\n";
        }
        else { print FASTA ">seq1\n$sequenceA\n"; }
        if ($strandB eq "-") {
                my @reverseB=&revcompl($sequenceB);
		print FASTA ">seq2\n@reverseB\n";
        }
        else { print FASTA ">seq2\n$sequenceB\n"; }
	my $sw_command = "$smithwaterman $sw_matrix $seqfasta";
	my $alignscore = `$sw_command`;
	print ANNOTATION "\tSIMSCORE:$alignscore";
	my $rm_fasta_cmd = "rm $seqfasta";
	system($rm_fasta_cmd);
}
sub reversestrand {
	if ($_[0] eq "+"){
		return "-";
	}
	elsif ($_[0] eq "-"){
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
	if ($_[1] eq "+") {
		$position1 = $_[0] -1;
	}
	elsif ($_[1] eq "-") {
		$position1 = $_[0] +1;
	}
	if ($_[3] eq "+") {
		$position2 = $_[2] +1;
	}
	elsif ($_[3] eq "-") {
		$position2 = $_[2] -1;
	}
	return ($position1, $position2); 
}
sub unadjustposition { #reverse the effect of adjustposition
	#takes in 0pos, 1strand, 2pos 3strand
	my $position1;
        my $position2;
        if ($_[1] eq "+") {
                $position1 = $_[0] +1;
        }
        elsif ($_[1] eq "-") {
                $position1 = $_[0] -1;
        }
        if ($_[3] eq "+") {
                $position2 = $_[2] -1;
        }
        elsif ($_[3] eq "-") {
                $position2 = $_[2] +1;
        }
        return ($position1, $position2);
}
