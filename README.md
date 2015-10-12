starchimp
==========
starchimp is short for Star Chimeric Post, written by Kipp Akers as part of his work in Bojan Losic's group at the Icahn Institute of Genomics and Multiscale Biology of Mount Sinai.   

This software is designed to take the chimeric output from the STAR alignment tool and discover high confidence fusions and circular RNA in the data. 
Before running, you must have used a recent version of STAR with chimeric output turned on, to align your RNA-Seq data.

Currently, there are two main modules, which need to be run separately. 

###############
##  Fusions  ##
###############

Usage:

	/path/to/StarChimPo/fusions/fusions-from-star.pl output_seed Chimeric.junction.out Parameters.txt
	
		This will output to the directory from which the script is run. 

Softare Dependencies:

	samtools
	bedtools
	mafft

Files Needed:
	
	parameters file (examples included in /starchimp/fusions/paramfiles )
	reference_fasta_file (indexed with samtools)
	reference gtf file in .bed format
	repetative elements in .bed format

Output:

	you'll get a *.annotated and a *.annotated.summary file
	The summary file is further filtered, and is more compact, while the annotated file has more information.  

Known Issues:

	Fusions are really tricky.  This still has a high false positive rate, due to short overhangs and  paralogues/pseudogenes. 
	A good second step is to BLAST/BLAT your fusion consensus sequence, especially to the newest genome assembly available. 
	In the future I'll implement BLAST on-the-fly to check that fusion partners don't share sequence similarity. 


##############
## Circles ##
##############
Circular RNA is not poly-A modified, so generally data from poly-T amplified RNA is not appropriate for this.    

Usage:

	/path/to/StarChimPo/circles/circle_star.sh [Reads_threshold] [Subjects_threshold] star_dirs.txt [splice/nosplice]

		[Reads_threshold] is the minimum read support required for a circleRNA.  I.e. 5
		[Subjects_threshold] is the minimum subjects support required I.e. 2
		Where star_dirs.txt is a file with the full path to a star output directory on each line.  These should be distinctive directories.  Ie /path/to/star_subject1_run1/  
		the 4th argument should be the word 'splice' or 'nosplice' .  splice will cause the software to run splice analysis, generating information about the likely splices within each cRNA. This can take significant time. 
		This will output to the directory from which the script is run. 

Software Dependencies:
	
	R with the following packages: limma, edgeR

Files Needed:

	-

Resources:

	The default is to use 4 CPUs.  This can be modified easily in the circle_star.sh script.  Run times vary and scale up with more samples.  Expect at least 1 minute per sample. 

Output:

	Count matrix: circRNA.cutoff[readthreshold]reads.[subjectthreshold]ind.countmatrix
	Info about each circRNA:  Circs[reads].[subjects].spliced.consensus
	There are other files as well, but they are mostly useful only if you need to really dig into the data. 

Known issues:

	Haven't really dug into how well this utilizes paired-end data
