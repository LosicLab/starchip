#! /usr/bin/perl
use warnings;
use strict;
use Cwd 'abs_path';
use Getopt::Long;
use File::Basename;

#This script will read in user parameters and feed them on to circle-star.sh.  
#Usage: starchimp-circles.pl list_of_star_dirs.txt params.txt

if (scalar(@ARGV) != 2 ) { die "Wrong number of inputs. Usage: starchimp-circles.pl star_dirs.txt params.txt \n Be sure you have R (Rscript) available.\n";}
#file management
my $script_dir=abs_path($0);
$script_dir =~ s/\/starchimp-circles.pl//;

##Read in User Parameters (taken from Perl Cookbook "8.16. Reading Configuration Files")
my %Configs = ();
my $configfile = $ARGV[1];
open CONFIG, "<$configfile" or die $!;
while (<CONFIG>) {
    chomp;                  # no newline
    s/#.*//;                # no comments
    s/^\s+//;               # no leading white
    s/\s+$//;               # no trailing white
    next unless length;     # anything left?
    my ($var, $value) = split(/\s*=\s*/, $_, 2);
    $Configs{$var} = $value;
}

#check if the cpus is unlisted/automatic:
if ($Configs{cpus} =~ /^\d+$/ ) {#this should see if we have a positive integer
	#do nothing.  
}
else { 
	print "Non-numeric value $Configs{cpus} for cpus.  Automatically assessing CPUs available: ";
	my $foundcpus=`nproc`;
	chomp($foundcpus);  
	$Configs{cpus} = $foundcpus - 1; 
	print "$foundcpus found, will use $Configs{cpus}\n"; 
}
if (lc $Configs{do_splice} eq "true") { $Configs{do_splice}="true"; }
elsif (lc $Configs{do_splice} eq "false") { $Configs{do_splice}="false"; }
else { 
	print "do_splice parameter should be true or false\n";
	die;
}
if (lc $Configs{annotate} eq "true") { 
	$Configs{annotate}="true"; 
	unless (-e $Configs{refbed}) { #unless the refbed file exists
 		die "Could not find file $Configs{refbed} needed for cRNA gene annotation\n";
 	} 
}
elsif (lc $Configs{annotate} eq "false") { 
	$Configs{annotate}="false"; 
	$Configs{refbed}="dummy_not_annotating"; 
}
else { 
	print "annotate parameter should be true or false\n";
	die;
}
if ($Configs{starprefix} eq "" ) { 
	$Configs{starprefix}="NoPrefix123456789"
}

my $run_cmd = "$script_dir/circle_star.sh $Configs{readsCutoff} $Configs{minSubjectLimit} $ARGV[0] $Configs{do_splice} $Configs{cpus} $Configs{cpmCutoff} $Configs{subjectCPMcutoff} $Configs{annotate} $Configs{refbed} $Configs{starprefix} $Configs{IDstepsback}" ; 
print "$run_cmd\n"; 
system($run_cmd); 

