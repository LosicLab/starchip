STARChip
==========
STARChip is short for Star Chimeric Post, written by Dr. Nicholas Kipp Akers as part of his work in Bojan Losic's group at the Icahn Institute of Genomics and Multiscale Biology at Mount Sinai School of Medicine

This software is designed to take the chimeric output from the STAR alignment tool and discover high confidence fusions and circular RNA in the data. 
Before running, you must have used a recent version of STAR with chimeric output turned on, to align your RNA-Seq data.

Please see http://starchip.readthedocs.io/en/latest/ for complete STARChip documentation

Currently, there are two main modules, which need to be run separately. 

## Setup ##

Both Fusions and Circles work from a completed STAR run with chimeric output turned on.  I STRONGLY recommend using the same fasta and gtf files to build your STAR index and run STARChip.    
You will need to run /path/to/starchip/setup.sh to create gtf files in bed format as well as .genome files:
	
		Usage: setup.sh /path/to/my.gtf  /path/to/my.fasta /path/to/desired/output_directory/

This will create a a my.gtf.bed file and a my.gtf.exons.bed file.  The only difference is that the 2nd file is restricted to gtf lines that have exon information, excluding annotations like full transcripts.  This will only affect your annotations, but I reccomend the exons file for the best annotations.   

Addionally, annotated repeats should be downloaded from the UCSC table browser.  See the manual for details.  

You will also need to create parameters files for your run.  See examples in /starchip/paramfiles/ or the manual for complete details.

##  Fusions  ##

Usage:

	/path/to/starchip/starchip-fusions.pl output_seed Chimeric.junction.out Parameters.txt
	
		This will output to the directory from which the script is run. 

Softare Dependencies:

	samtools
	bedtools (>= 2.24.0)
	mafft
	star (starchip doesn't call star, but you do need star output)

Files Needed:
	
	parameters file (examples included in /starchip/paramfiles )
	Fasta File and GTF File used to generate STAR index 
	repetative elements/antibody parts regions files in .bed format (included for hg19,hg38)
	Optional: known gene families and paralogs (included for human)

Output:

	you'll get a [output_seed].annotated and a [output_seed].annotated.summary file
	The summary file is more compact, while the annotated file has more information.  

Known Issues:

	Fusions are tricky.  This still has a high false positive rate, due to short overhangs and paralogues/pseudogenes. 
	A good second step is to BLAST/BLAT your fusion consensus sequence, especially to the newest genome assembly available. 


## Circles ##

Circular RNA is not poly-A modified, so generally data from poly-T amplified RNA is not appropriate for this.    

Usage:

	/path/to/starchip/starchip-circles.pl star_dirs.txt parameters.txt 

	This command will generate scripts Step1.sh, Step2.sh, and optionally Step3.sh and Step4.sh  This is to allow users to optimize computing resources for a variety of computing environments.  
	starchip-circles is designed to work on a set of multiple samples (though it works fine on 1 sample). 
	star_dirs.txt should have 1 full pathway to star output per line: ie 
		/path/to/file1/star/
		/path/to/file2/star/
	Alternatively, if you want STARChip to perform STAR realignment, star_dirs.txt should have fastq files:
		/path/to/file1.1_R1.fastq,/path/to/file1.2_R1.fastq  /path/to/file1.1_R2.fastq,/path/to/file1.2_R2.fastq	
		/path/to/file2.1_R1.fastq,/path/to/file2.2_R1.fastq  /path/to/file2.1_R2.fastq,/path/to/file2.2_R2.fastq	


Run Modes: 

	STARChip circRNA can be run on just STAR outputs, or it can work directly from fastq.  If running from fastq files, set runSTAR=true in your parameters file.  
	If runSTAR=true, STARChip will create scripts to align (Step1.sh), discover circRNA (Step2.sh), re-align (Step3.sh), and quantify/annotate (Step4.sh).
	If runSTAR=false, STARChip will create scrips to discover (Step1.sh), and quantify/annotate (Step2.sh) 

Software Dependencies:
	
	R with the following packages: limma, edgeR
	bedtools (>=2.24.0)
	star

Files Needed:

	GTF file with the same chromosome naming conventions as your star alignments (in bed format via setup.sh).  
	
Resources:

	starchip-circles can make use of multiple processors, modifiable in the paramters file.  Run times vary and scale up with more samples.  Expect at least 1 minute per sample. 

Output:

	Count matrixes : 
		raw cRNA backsplice counts: circRNA.cutoff[readthreshold]reads.[subjectthreshold]ind.countmatrix
		log2CPM of above: norm_log2_counts_circRNA.[readthreshold]reads.[subjectthreshold]ind.0cpm_0samples.txt 
		Maximum Linear Splices at Circular Loci: rawdata/linear.[readthreshold]reads.[subjectthreshold]ind.sjmax
	Info about each circRNA:  
		Consensus Information about Internal Splicing: Circs[reads].[subjects].spliced.consensus
		Complete Gene Annotation:  circRNA.[readthreshold]reads.[subjectthreshold]ind.annotated
		Consise Gene Annotation + Splice Type:  circRNA.[readthreshold]reads.[subjectthreshold]ind.genes
	Images:
		PCA plots: circRNA.[readthreshold]reads.[subjectthreshold]ind.0cpm_0samples_variance_PCA.pdf
		Heatmap: circRNA.[readthreshold]reads.[subjectthreshold]ind.heatmap.pdf 
	There are other files as well stored in rawdata/ but they are mostly useful only if you need to really dig into the data. 



