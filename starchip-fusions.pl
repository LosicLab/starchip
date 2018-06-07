#! /usr/bin/perl
use warnings;
use strict;
use Cwd 'abs_path';
use Getopt::Long;
use File::Basename;
## usage: fusions-from-star.pl  outputname Chimeric.out.junction  

##IMPORTANT NOTE!  this calls 'samtools' and 'bedtools' (>= v.2.24) and 'mafft'  please have these installed and in your path under those aliases.  
	#It shouldn't crash without them, but you will get error messages rather than some of the outputs.
	#$ On MSSM/minerva you can just do: module load starchip

# to do:    
	# blat/blast output.  blat is REALLY slow, and because of mem loading, probably faster if run in batches.  
		# check-pseudogenes.sh is a work in progress, but should work eventually, once I get gtf-->fasta parsing correct. 
		# run all blat at once.
		# $ blat t=dna q=rna /path/to/database /path/to/query 
		# the other issue is that blat/blast don't lend themselves to computational filtering easily.  taking the top hit works ok though.   

if (scalar(@ARGV) != 3 ) { die "Wrong number of inputs. Usage: starchip-fusions.pl output_seed input_chimeric.out.junction params.txt \n Be sure you have samtools, bedtools, and mafft available.\n";}

#Output the Verion
my $CHANGES=abs_path($0) ;
$CHANGES =~ s/starchip-fusions.pl/CHANGES.md/;
open CHANGESFILE, "<$CHANGES" or die $! ;
my $version = <CHANGESFILE> ;
close CHANGESFILE ;
print "STARChip $version";

my $numbcolumns=14; #need this one in case you junction.out file/input changes.  this should be a non-existant final column (ie Chimeric.junciton.out has 0-13 columns)
#should fix this to be more robust...

#file management
my $script_dir=abs_path($0);
$script_dir =~ s/starchip-fusions.pl/scripts\/fusions\//;
my $outbase = $ARGV[0];
my $outsumm = $outbase . ".summary";
my $outsummtemp = $outsumm . ".temp";
my $outannotemp = $outsummtemp . ".annotated" ;
my $outanno = $outsumm . ".annotated";
print "your final outputs will be in $outsumm and $outsumm.annotated\n";

#RUN FUSIONS FROM STAR
	#needs most values.  
	my $fusions_from_star_cmd = $script_dir . "fusions-from-star.pl" . " $ARGV[0] " . $ARGV[1] . " " . $ARGV[2] ; 
	print "$fusions_from_star_cmd\n";
	system("$fusions_from_star_cmd"); 
#RUN ANNOTATE-FUSIONS
	#check that there are more than just the header line
	open (OUTSUMTEMP, "<$outsummtemp") or die "$! cannot open $outsummtemp\n"; 
	while (<OUTSUMTEMP>) {}; #stores # of lines in $. 
	if ($. > 0) {
		my $annotate_fusions_cmd = $script_dir . "annotate-fusions.pl" . " $ARGV[0] " . $ARGV[1] . " " . $ARGV[2] ;
		print "$annotate_fusions_cmd\n"; 
		system("$annotate_fusions_cmd"); 	
	}
#cleanup
	my $cleanupcmd = "rm -f $outsummtemp $outannotemp";
	system($cleanupcmd); 
