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
$script_dir =~ s/fusions-from-star.pl//;
my $consensusloc= $script_dir . 'consensus.sh';
my $annotateloc= $script_dir . 'coordinates2genes.sh';

#file management
if (length(@ARGV) != 1) { die "Wrong number of arguments!\n";}
my $outbase = $ARGV[0];
my $outsumm = $outbase . ".summary";
print "your final outputs will be in $outsumm and $outsumm.annotated\n";
my $junction = $ARGV[1];
my $sam = $junction;
$sam =~ s/junction$/sam/ ; 

# primary data format; $chr1_pos1_chr2_pos2_strandA_strandB[0/1/2][0-RL]
	#where strand is + or -
	#and [0/1/2] indicates left-side/right-side of fusion/non-split reads
	#0-Read length starts at the fusion site =0 and expands outwards from there, to pos2+RL on right side, pos1-RL on left side.
	# in most cases, we will use + strand notation for positions.

open JUNCTION, "<$junction" or die $!; 
while (my $x = <JUNCTION>) {
	my @line = split(/\s+/, $x); 
##some filtering 
	next if ($line[$col_chrA] eq "chrM");
	next if ($line[$col_chrB] eq "chrM");
	next if ($line[$col_chrA] eq "MT");
	next if ($line[$col_chrB] eq "MT");
	next if ($line[$col_chrA] =~ m/GL*/);
	next if ($line[$col_chrB] =~ m/GL*/);
	next if ($line[$col_overlapL] > $overlapLimit);
	next if ($line[$col_overlapR] > $overlapLimit);
	#skipping problem regions (Ig antibody parts regions).  This could be improved for easily changing genomes...
	next if ($line[$col_chrA] eq "chr2" && $line[$col_chrB] eq "chr2" && $line[$col_FusionposA] >= 89156874 && $line[$col_FusionposA] <=90471176 && $line[$col_FusionposB] >= 89156874 && $line[$col_FusionposB]<=90471176);
	next if ($line[$col_chrA] eq "chr22" && $line[$col_chrB] eq "chr22" && $line[$col_FusionposA] >= 22385572 && $line[$col_FusionposA] <=23265082 && $line[$col_FusionposB] >= 22385572 && $line[$col_FusionposB]<=23265082);
	next if ($line[$col_chrA] eq "chr14" && $line[$col_chrB] eq "chr14" && $line[$col_FusionposA] >= 105994256 && $line[$col_FusionposA] <=107283085 && $line[$col_FusionposB] >= 105994256 && $line[$col_FusionposB]<=107283085);
	
	next if ($line[$col_chrA] eq "2" && $line[$col_chrB] eq "2" && $line[$col_FusionposA] >= 89156874 && $line[$col_FusionposA] <=90471176 && $line[$col_FusionposB] >= 89156874 && $line[$col_FusionposB]<=90471176);
	next if ($line[$col_chrA] eq "22" && $line[$col_chrB] eq "22" && $line[$col_FusionposA] >= 22385572 && $line[$col_FusionposA] <=23265082 && $line[$col_FusionposB] >= 22385572 && $line[$col_FusionposB]<=23265082);
	next if ($line[$col_chrA] eq "14" && $line[$col_chrB] eq "14" && $line[$col_FusionposA] >= 105994256 && $line[$col_FusionposA] <=107283085 && $line[$col_FusionposB] >= 105994256 && $line[$col_FusionposB]<=107283085);
	#skipping other long genes
	next if ($line[$col_chrA] eq "chr3" && $line[$col_chrB] eq "chr3" && $line[$col_FusionposA] >= 149530000 && $line[$col_FusionposA] <=149680000 && $line[$col_FusionposB] >= 149530000 && $line[$col_FusionposB]<=149680000);
	next if ($line[$col_chrA] eq "3" && $line[$col_chrB] eq "3" && $line[$col_FusionposA] >= 149530000 && $line[$col_FusionposA] <=149680000 && $line[$col_FusionposB] >= 149530000 && $line[$col_FusionposB]<=149680000);
	
#calculate read length (where read length == the length of one pair of the sequencing if paired end)
	if ($linecount < 1 ) {
		if ($line[$col_chrA] =~ /chr/ ) { $chrflag = 1; }
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
		else { print "read length error, please check input\n"; 
			$readlength=100;
			$linecount=-1;
		}
	}
#Create two 'names' and check if they exist.  Checking the written + reverse allows us to collapse reads on opposite strands  
	my $fusionname=$line[$col_chrA] . "_" . $line[$col_FusionposA] . "_" . $line[$col_chrB] . "_" . $line[$col_FusionposB] . "_" . $line[$col_strandA] . "_" . $line[$col_strandB] ; 
	my $invStrandA = &reversestrand($line[$col_strandA]);
	my $invStrandB = &reversestrand($line[$col_strandB]);
	my $fusionnameInv=$line[$col_chrB] . "_" . $line[$col_FusionposB] . "_" . $line[$col_chrA] . "_" . $line[$col_FusionposA] . "_" . $invStrandB . "_" . $invStrandA ;
	## chrA_pos1_+ fused to chrB_pos2_+ equals chrB_pos2_- fused to chrA_pos_- etc.  
	##check the existence
	if (exists $fusions{$fusionname}) {
		&supportCigar(@line, $fusionname, "1");

	}
	elsif (exists $fusions{$fusionnameInv}) {
	  #because this is the reverse compliment of the already indexed read, we'll feed in a rearranged line from Chimeric.out.junction, with the strands flipped.
		my $chrA=$line[$col_chrB];
		my $posA=$line[$col_FusionposB];
		my $strandA = &reversestrand($line[$col_strandB]); 
		my $chrB=$line[$col_chrA];
		my $posB=$line[$col_FusionposA];
		my $strandB = &reversestrand($line[$col_strandA]);
		#6 7,8,9 unchanged/unimportant here
		my $startposA = $line[$col_startposB];
		my $cigarA = $line[$col_cigarB];
		my $starposB = $line[$col_startposA];
		my $cigarB = $line[$col_cigarA];
		#print "@line\n";
		#print "$chrA $posA $strandA $chrB $posB $strandB $line[7] $line[8] $line[9] $startposA $cigarA $starposB $cigarB\n";
		&supportCigar($chrA, $posA, $strandA, $chrB, $posB, $strandB, $line[6], $line[7], $line[8], $line[9], $startposA, $cigarA, $starposB, $cigarB, $fusionnameInv, "2" );
	}
	else {
		#create the array in the fusions hash
		for my $x (0..($readlength-1)) {
			$fusions{$fusionname}[0][$x] =0; #jxn crossing read support (Left side)
			$fusions{$fusionname}[1][$x] =0; #jxn crossing read support (Right side)
		}
		$fusions{$fusionname}[2][0] = 0; #jxn spanning read support on the strands as named
		$fusions{$fusionname}[2][1] = 0; # strand distribution (jxn crossing reads in 'fusionname' orientaiton with first listed chrm on the left)
		$fusions{$fusionname}[2][2] = 0; # strand distribution (jxn crossing reads in 'fusionnameInv' orientation with first listed chrm on the right)
		$fusions{$fusionname}[2][3] = 0; # A chrom anchored (for split reads, the pair lies on chromosome A)
		$fusions{$fusionname}[2][4] = 0; # B chrom anchored (for split reads, the matepair lies on chromosome B)
		$fusions{$fusionname}[2][5] = 0; #jxn spanning read support inverse strands
		&supportCigar(@line, $fusionname, "1");
	}
	$linecount++;
	#print "$fusionname\t$fusions{$fusionname}[0][0]\t$fusions{$fusionname}[1][0]\n";
	#if ($linecount > 10 ) { die; }
}
my $numbkeys = scalar keys %fusions;
print "Finished catologing fusion reads, now processing over $numbkeys\n";
close(JUNCTION);

##Post Processing and Filtering
#SETUP
my $fusioncounter =0;
open (SUMM, ">$outsumm") or die; 
print SUMM "Partner1\tPartner2\tScore\tSpanningReads\tSplitReads\tTopsideCrossing\tBottomsideCrossing\tChromAAnchors\tChromBAnchors\tUniqueSupportLeft\tUniqueSupportRight\tKurtosis\tSkew\tLeftAnchor\tRightAnchor\tTopsideSpanning\tBottomsideSpanning\tReferenceSeq\tConsensusSeq\n";
my $keycount ;

#Join and Evaluate Fusions
for my $key (keys %fusions) {#go through all 'fusions'
	$keycount++ ;
	#first filter by read support
	if ($fusions{$key}[0][0] >=$cutoff && $fusions{$key}[1][0] >=$cutoff) {
                my @keyarray=split(/_/, $key); #0:chrm1 1:pos 2:chr2 3:pos2 4:strand 5:strand 
                 #skip same-chrom proximal fusions
		if ($keyarray[0] ~~$keyarray[2] && (abs($keyarray[1]-$keyarray[3])<=$samechrom_wiggle) ){
               		next;
                }
		#skip lopsided fusions (ie all the read support is on one side)
		if ( ($fusions{$key}[2][1]+0.1)/($fusions{$key}[2][2]+0.1) >= $lopsidedupper || ($fusions{$key}[2][1]+0.1)/($fusions{$key}[2][2]+0.1) <= $lopsidedlower) {
			next;
		}
		my $topspancount = $fusions{$key}[2][0];
		my $bottomspancount = $fusions{$key}[2][5];
		##now we need to do a broad check for spanning reads that didn't map to the exact fusion location
		# I do this by doing a rough check for close locations.  This mostly works, but gets messy with larger values of $wiggle, and when fusion sites are near each other on the same chromosome.  
		# in those cases, take results with a grain of salt.
		if ($pairedend == 1) {  
			for my $key2 (keys %fusions) { #cycle through all the keys again
				next if ($key ~~$key2); #skip itself
				if ($fusions{$key2}[2][0] >= 1 || $fusions{$key2}[2][5] >=1 ) { #we are only really interested in sites with spanning fusions
					my @key2array=split(/_/, $key2); #see above for indices
					next if ($key2array[0] ~~$key2array[2] && (abs($key2array[1]-$key2array[3])<=$samechrom_wiggle) ); #skip same chrom proximal
					#check if the fusion from keyarray (has jxn crossing) has the same coordinates (within $wiggle bp) as $keyarray2
					if ($keyarray[4]~~$key2array[4] && $keyarray[5]~~$key2array[5] && $keyarray[0]~~$key2array[0] && $keyarray[2]~~$key2array[2] && (abs($keyarray[1]-$key2array[1])<=$wiggle) && (abs($keyarray[3]-$key2array[3])<=$wiggle) ) {
						#print "1: $key2 spans $key reads: $fusions{$key2}[2][0] sum before: $fusions{$key}[2][0]\n";
						$topspancount += $fusions{$key2}[2][0];
						$bottomspancount += $fusions{$key2}[2][5]; 
					}
					#check the inverse fusion
					elsif (&reversestrand($keyarray[4])~~$key2array[5] && &reversestrand($keyarray[5])~~$key2array[4] && $keyarray[0]~~$key2array[2] && $keyarray[2]~~$key2array[0] && (abs($keyarray[1]-$key2array[3])<=$wiggle) && (abs($keyarray[3]-$key2array[1])<=$wiggle) ) {
						$topspancount += $fusions{$key2}[2][5];
						$bottomspancount += $fusions{$key2}[2][0];
						#print "2: $key2 spans $key reads: $fusions{$key2}[2][0] sum before: $fusions{$key}[2][0]\n";
					}
				}
			}
		}
		#Second, filter on spanning read pairs
		my $spancount = $topspancount + $bottomspancount ; 
		if ($spancount >= $spancutoff) { 
			#next get an estimate of the 'unique reads' mapped by looking at how many unique read support values there are
			my @array0; my @array1;
			my $leftanchor ; my $rightanchor; 
	                for my $x (0..($readlength-1)) { 
				#cycle across possible support positions.  create arrays of # of reads of support at each position.  (skip 0 read support) Also note largest overlap in read support
   	                	if ($fusions{$key}[0][$x] != 0) {
					push (@array0, $fusions{$key}[0][$x]);	$leftanchor = $x;}
				if ($fusions{$key}[1][$x] != 0) {
					push (@array1, $fusions{$key}[1][$x]); $rightanchor = $x;}
        	        }
                	my %count0; my %count1;
        	        @count0{@array0} =(); @count1{@array1} =(); #turn the arrays into hashes.
	                my $unique0=scalar keys %count0; my $unique1=scalar keys %count1; #count the unique hash indices. 
			my @kurtosisarray;
               		for my $x (reverse (15..($readlength-16))) { push (@kurtosisarray, $fusions{$key}[0][$x]);}
			for my $x (15..($readlength-16)) { push (@kurtosisarray, $fusions{$key}[1][$x]);}
			my ($skew, $kurtosis)=&kurtosis(@kurtosisarray); #still calculating kurtosis, but not sure how useful it is.  
			#filter by the number of these unique reads
                	if ($unique0 >= $cutoff2 && $unique1 >= $cutoff2) {
				$fusioncounter++;
				my $splitreads = $fusions{$key}[2][1] + $fusions{$key}[2][2] ;
				my ($position1, $position2) = &adjustposition($keyarray[1],$keyarray[4],$keyarray[3],$keyarray[5]); 
				#0:split reads, 1:topsidesplit, 2:bottomsidesplit 3:spanreads 4;topspan 5;bottomspan 6:skew 7:chr1 8;loc1 9;strand1 10;chr2; 11;loc2; 12;strand2
				#print "$splitreads,$fusions{$key}[2][1],$fusions{$key}[2][2],$spancount,$topspancount, $bottomspancount,$skew,$keyarray[0],$position1,$keyarray[4],,$keyarray[2],$position2,$keyarray[5])\n";
				my $score = &fusionScore($splitreads,$fusions{$key}[2][1],$fusions{$key}[2][2],$spancount,$topspancount, $bottomspancount,$skew,$keyarray[0],$position1,$keyarray[4],$keyarray[2],$position2,$keyarray[5]);
				
				##Output. This can be changed as needed, but the first two columns need to be chr1:pos:str.  They are fed into coordinates2genes.sh for gene annotation later.  
				if ($chrflag == 0 ) {
					print SUMM "chr$keyarray[0]:$position1:$keyarray[4]\tchr$keyarray[2]:$position2:$keyarray[5]\t$score\t$spancount\t$splitreads\t$fusions{$key}[2][1]\t$fusions{$key}[2][2]\t$fusions{$key}[2][3]\t$fusions{$key}[2][4]\t$unique0\t$unique1\t$kurtosis\t$skew\t$leftanchor\t$rightanchor\t$topspancount\t$bottomspancount\t";
				}
				else {
					print SUMM "$keyarray[0]:$position1:$keyarray[4]\t$keyarray[2]:$position2:$keyarray[5]\t$score\t$spancount\t$splitreads\t$fusions{$key}[2][1]\t$fusions{$key}[2][2]\t$fusions{$key}[2][3]\t$fusions{$key}[2][4]\t$unique0\t$unique1\t$kurtosis\t$skew\t$leftanchor\t$rightanchor\t$topspancount\t$bottomspancount\t";
				}
				&extractSequence($key);	
				print "Filtered fusions count:$fusioncounter, searched $keycount\n";
			}
		}
	}
}
print "Total fusions found: $fusioncounter\nNow Annotating these Fusions\n";
close (SUMM);

##Post filter gene-annotation here:
my $annotate_command = $annotateloc . " " . $outsumm . " " . $refbed . " " . $repeatbed; 
system($annotate_command); 
#should create a file $outsum.annotated with gene annotations.
#remove same gene fusions, simplify output
my $outanno = $outsumm . ".annotated" ; 
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
		if (($gene1name[0] eq "HLA-B" && $gene2name[0] eq "HLA-C") || ($gene1name[0] eq "HLA-C" && $gene2name[0] eq "HLA-B")) { next; }
		if (($gene1name[0] eq "FTL" && $gene2name[0] eq "FTLP3") || ($gene1name[0] eq "FTLP3" && $gene2name[0] eq "FTL")) { next; }
		if (($gene1name[0] eq "CES1" && $gene2name[0] eq "CES1P1") || ($gene1name[0] eq "CES1P1" && $gene2name[0] eq "CES1")) { next; }
		if (($gene1name[0] eq "CES1" && $gene2name[0] eq "CES1P1") || ($gene1name[0] eq "CES1P1" && $gene2name[0] eq "CES1")) { next; }
		my $score = $line[2]; 
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
