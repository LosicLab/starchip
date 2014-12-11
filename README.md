starchimp
==========
starchimp is short for Star Chimeric Post, written by Kipp Akers.  

This software is designed to take the chimeric output from the STAR alignment tool and discover high confidence fusions and circular RNA in the data. 
Before running, you must have used a recent version of STAR with chimeric output turned on, to align your RNA-Seq data.

Currently, there are two main modules, which need to be run separately. 

###############
##  Fusions  ##
###############

Usage:

/path/to/StarChimPo/fusions/fusions-from-star.pl output_seed Chimeric.junction.out


This will output to the directory from which the script is from. 

Softare Dependencies:

	samtools
	bedtools
	razerS3

Files Needed:

	reference_fasta_file (indexed with samtools)
	reference gtf file in .bed format
	repetative elements in .bed format

##############
## Circles ##
##############
Circular RNA is not poly-A modified, so generally data from poly-T amplified RNA is not appropriate for this.    

Usage:

	/path/to/StarChimPo/circles/circle_star.sh [Reads_threshold] [Subjects_threshold] star_dirs.txt

		[Reads_threshold] is the minimum read support required for a circleRNA.  I.e. 5
		[Subjects_threshold] is the minimum subjects support required I.e. 2
		Where star_dirs.txt is a file with the full path to a star output directory on each line.  These should be distinctive directories.  Ie /path/to/star_subject1_run1/  
		This will output to the directory from which the script is run. 

Software Dependencies:
	
	R with the following packages: limma, edgeR

Files Needed:

	-

Resources:

	It will use 1 cpu per sample.  Run times vary but an hour is typical. 

Output:

	Count matrix: circRNA.cutoff[readthreshold]reads.[subjectthreshold]ind.countmatrix
	Info about each circRNA:  Circs[reads].[subjects].spliced.consensus
	There are other files as well, but they are mostly useful only if you need to really dig into the data. 

Known issues:

	Haven't really dug into how well this utilizes paired-end data
