#! /usr/bin/perl

# Sequence Alignment Program by Kipp Akers
# This program will align 2 sequences based on a scoring matrix provided by the user.
#running this program from the command prompt should be performed as follows: align3.pl SCOREFILENAME SEQUENCEFILENAME
#Your score file should be in the following format: Top line: all possible letters to be aligned separated by whitespace (ie" A T C G ")
#Each line after that should provide the alignment score for the letter above and the corresponding letter on a transposed axis. Finally, the 'gap penalty should be on the last line.
# For example: 
#A T C G
#1 0 0 0
#0 1 0 0
#0 0 1 0
#0 0 0 1
#-4
#Your sequence file should be FASTA formatted with at least 2 sequences
#One output will be the scoring matrix from your two sequences.  I decided that the matrix is best viewed in a text editor rather than in the command prompt, so the similarity score
#matrix is output to a text file called 'matrix.txt'. I advise opening this in a text editor that does not have a word-wrap function enabled. It could also be viewed handily in open office calc. 
#The other form of output, the alignment, will show up in your command prompt, 80 characters at a time.  
#Enjoy!  



##########Reading in the scoring matrix##############
#########Featuring: lots of error checking###########
open (SCOREFILE, $ARGV[0]); 
	if (! open SCOREFILE, $ARGV[0]) {
		die "Cannot find or open your scoring matrix file: $! \n";
	}
	$s = 0;
	while ($line = <SCOREFILE>) { #read each line of the score file to an array
		chomp $line;
		$scorein[$s] = lc($line);
		$s += 1;
	}
	$possibles = $scorein[0]; #take the top line
	$possibles =~ s/\s//g; #remove the whitespace
	$gap = $scorein[$#scorein]; #assign the gap score as the last line of the array
	if ($gap =~ /^[+-]?\d+$/ ) {} # if the gap score is not a number (restricted to integers here), then die
	else {
			die "Program aborted because your gap penalty is not an integer";
	}
	if ($possibles =~ /[+-]?\d+$/) { #if the top line has numbers then die
		die "Program aborted because the top line of your score matrix contiains a numerical value\n";
	}


	foreach $k (0..$#scorein) {#iterate through the lines of the score matrix input array
		@scorearray = split (/\s/, $scorein[$k]); #split each line into an array matching on (removing) whitespace
		$matlength[$k] = scalar(@scorearray); #create a length array to check our matrix is even
		if ($k != $#scorein) { #while it's not the gap penalty line
			if ($matlength[$k] != $matlength[0]) {#check all lines are equal in length to the top line
				die "Program aborted:One of the lines of your score matrix is shorter/longer than the top line\n";#if not then die.
			}
		}	
		foreach $x (0..$#scorearray) { #then iterate through the freshly minted array
			$scoremat[$k][$x]= $scorearray[$x];#and assign the value to an array of arrays (our score matrix)
			
		} 
	}
	$length = length($possibles)-1; #this is the length of the score matrix, or the number of different possibilities. the -1 is because perl counts from 0. 
	foreach $x (0..$length) { #iterate through the top line
		foreach $y (0..$length) { #again iterate through the top line
			$alignpair = join ('', $scoremat[0][$x],$scoremat[0][$y]);#create pairs, ie AA or TG
			$scorehash {$alignpair} = $scoremat[$x+1][$y]; #assign the pairs to their score in the matrix via a hash
		}	
	}
if ($matlength[0] != ($#scorein-1)){ #if the length of your top line is not equal to the number of rows in the matrix (except the gap penalty row)
	die "Program aborted because your scoring matrix is not even\n"; #die and give an error
}	
close SCOREFILE;


###### Opening the file containing sequences ######################
open (ALIGNFILE, $ARGV[1]); #Open the file 
	if (! open ALIGNFILE, $ARGV[1]){
		die "Cannot find or cannot open your file containg sequences to be aligned: $! \n";
		}
	$inputcount = 0; #the count will increase with each sequence
	while ($line = <ALIGNFILE>) {
		chomp $line;
		if ($line =~ /^>/) {#if a line looks like a name
				$inputarray[$inputcount][0] = $name;#put the previous name into a name/sequence matrix
				$inputarray[$inputcount][1] = lc($sequence);#put the sequence in
				$inputcount++; #increase the count
				$name = $line; #assign the new name
				$sequence = '';#clear the sequence
			
		}
		else {#if it looks like a sequence
			$sequence .= $line;#concatenate
		}
		$inputarray[$inputcount][0] = $name; #assign the last name to the matrix
		$inputarray[$inputcount][1] = lc($sequence);#assign the last sequence to the matrix
	}
	if ($inputarray[1][1] =~ /[^$scorein[0]]/) {#if the first sequence has character NOT in the first line of the score file
		die "Program aborted because your FASTA sequence 1 contains characters not in the scoring matrix";
	}
	if ($inputarray[2][1] =~ /[^$scorein[0]]/) { #if the 2nd sequence has characters not in the first line of the score file
		die "Program aborted because your FASTA sequence 2 contains characters not in the scoring matrix";
	}	
close ALIGNFILE;


####Setting some variables/Book keeping### 

$ilength = length($inputarray[1][1])-1; #the length of the first sequence
$jlength = length($inputarray[2][1])-1; #length of the second sequence
$inputarray[1][1] = "*" . "$inputarray[1][1]"; # the * is a placeholder for the '0' spot
$inputarray[2][1] = "*" . "$inputarray[2][1]";
@iseq = split(//, $inputarray[1][1]); #split up seq1 characters into an array
@jseq = split(//, $inputarray[2][1]); #split up seq1 characters into an array

for ($i = 0, $tener = 0;$i <= $#iseq; $scorematrix[$i][0] = $tener, $tener += $gap, $i++){}; #set (1, $i) values to decreasing values of -10x.  
for ($j = 0, $tener = 0;$j <= $#jseq; $scorematrix[0][$j] = $tener, $tener += $gap, $j++){}; #set ($j, 1) values to decreasing values of -10x.


#####iterate through and create an alignment score matrix################
foreach $i (1..$#iseq) { 
	foreach $j (1..$#jseq) {
	&maxscore;
	}
}


##########Printout of the matrix ##################
#not necassary in starchimp
#open (MYFILE, ">matrix.txt"); #I'm printing the matrix to a tab delimited text file, I think that's more helpful for long sequences than output in the command prompt 80 characters at a time.
#print MYFILE"|*	";
#foreach $i (0..$#iseq) { #creating the top line of the score matrix printout
#	print MYFILE"|$iseq[$i]	";
#}
#	print MYFILE"|\n";
#foreach $j (0..$jlength){ #iterate through seq2 (the 'y axis')
#	print MYFILE"|$jseq[$j]	"; #print the row header (letters)
#	foreach $i (0..$ilength) { #then iterate through the 'x axis' (seq1)
#		print MYFILE"|$scorematrix[$i][$j]	"; #and print the score for that spot and a tab
#	}
#	print MYFILE"|\n";
#}
#close (MYFILE);
#print "Your max score is $scorematrix[$ilength][$jlength].  Your scoring matrix is called \"matrix.txt\"\n";
#print "ilength is $ilength, jlength is $jlength\n";
print "$scorematrix[$ilength][$jlength]"; #my output I'm after.  

#don't need to do any of the below.
exit(1);

#######Creating the alignment###############
$aligni = ''; #empty variables to be filled with alignment seqs
$alignj = '';
$itracker = $ilength;
$jtracker = $jlength;#starting from the bottom right
while ($itracker > 0 && $jtracker > 0) {#work until you hit the upper left
	##print "$itracker\t$jtracker\n";
	##sleep(0.25);
	if ($alignmatrix[$itracker][$jtracker] == 1) { #if the score came from diagonal
		$aligni = join ('', substr($inputarray[1][1], $itracker, 1),$aligni);#add a letter to seq1
		$alignj = join ('', substr($inputarray[2][1], $jtracker, 1),$alignj);#add a letter to seq2
		$itracker -= 1; #move left and
		$jtracker -= 1; #up in our alignmatrix

	}
	elsif ($alignmatrix[$itracker][$jtracker] == 3) {#if the score came from the left
		$aligni = join ('', '_',$aligni); #add a placeholder to seq1
		$alignj = join ('', substr($inputarray[2][1], $jtracker, 1),$alignj);#add a letter to seq2
		$jtracker -= 1;#move left in our alignmatrix
	}
	elsif ($alignmatrix[$itracker][$jtracker] == 2) {#if the score came from above
		$alignj = join ('', '_',$alignj);#add a spacer to sequence 2
		$aligni = join ('', substr($inputarray[1][1], $itracker, 1),$aligni);#add a letter to seq1
		$itracker -= 1; #move up in the alignmatrix
	}
}

##########Output of aligned Sequences##################
#print "Seq1 name: $inputarray[1][0]\n";
#print "Seq2 name: $inputarray[2][0]\n"; 
$lengthi = length($aligni);
$lengthj = length($alignj);
$i = 0;
#while ($i < $lengthi && $i < $lengthj) {
#	if ($i == 0) {
#		print "Seq1:";
#		print substr($aligni, $i, 75);
#		print "\n";
#		print "Seq2:";
#		print substr($alignj, $i, 75);
#		print "\n";
#		$i += 75
#	}	
#	else {
#		print substr($aligni, $i, 80);
#		print "\n";
#		print substr($alignj, $i, 80);
#		print "\n";
#		$i += 80;
#	}
#}



###Begin Subroutines###

sub maxscore {
	&alignscore; #start by figuring out the alignment score for our position
	$diaganol = $scorematrix[$i-1][$j-1] + $align; #create a diaganol score
	$up = $scorematrix[$i-1][$j] + $gap; #an up score
	$left = $scorematrix[$i][$j-1] + $gap; #and a left score
	@scorearray = ("$diaganol", "$up", "$left"); #put the scores in an array
	@sortedscorearray = sort { $a <=>$b } @scorearray; #sort by numerical value
	$scorematrix[$i][$j] = $sortedscorearray[2]; #whichever value was highest, throw it into our scorematrix
		if ($scorematrix[$i][$j] == $diaganol) {#this if/elsif area defines where a cell in scorematrix got its value from. 
			$alignmatrix[$i][$j] = 1; #1 for diagonal (up and left)
		}
		elsif ($scorematrix[$i][$j] == $up) {
			$alignmatrix[$i][$j] = 2; # 2 for above (space in sequence 2)
		}
		elsif ($scorematrix[$i][$j] == $left) { #3 for from the left (space in sequence 1)
			$alignmatrix[$i][$j] = 3; #alignmatrix will later be used to iterate back
		}
}
sub alignscore { #this subroutine creates a key value for our alignment score hash then returns the hash key value.  
	$iletter = substr($inputarray[1][1], $i, 1);
	$jletter = substr($inputarray[2][1], $j, 1);
	$pair = $iletter;
	$pair .= $jletter;
	$align = $scorehash{$pair};
}

