STARChip
==========
STARChip manuscript now published in Bioinformatics! Open Access [Here](https://academic.oup.com/bioinformatics/advance-article/doi/10.1093/bioinformatics/bty091/4883488?guestAccessKey=9f40eec1-96cc-4b0f-b0b5-bb1eaa2e20eb)

STARChip is short for Star Chimeric Post, written by Dr. Nicholas Kipp Akers as part of his work in Bojan Losic's group at the Icahn Institute of Genomics and Multiscale Biology at Mount Sinai School of Medicine

This software is designed to take the chimeric output from the STAR alignment tool and discover high confidence fusions and circular RNA in the data. 
Before running, you must have used a recent version of STAR with chimeric output turned on, to align your RNA-Seq data.

Please see http://starchip.readthedocs.io/en/latest/ for complete STARChip documentation

Currently, there are two main modules, which need to be run separately. 

## Setup ##

STARChip works hand-in-hand with the STAR aligner.  I STRONGLY recommend using the same fasta and gtf files to build your STAR index and run STARChip.    
You will need to run /path/to/starchip/setup.sh to create gtf files in bed format as well as .genome files:
	
		Usage: setup.sh /path/to/my.gtf  /path/to/my.fasta /path/to/desired/output_directory/

This will create a a my.gtf.bed file and a my.gtf.exons.bed file.  The only difference is that the 2nd file is restricted to gtf lines that have exon information, excluding annotations like full transcripts.  This will only affect your annotations, but I reccomend the exons file for the best annotations.   

For fusion analysis, annotated repeats should be downloaded from the UCSC table browser.  See the manual for details.

You will also need to create parameters files for your run.  See examples in /starchip/paramfiles/ or the manual for complete details.

## Circles ##

Circular RNA is not poly-A modified, so generally data from poly-T amplified RNA is not appropriate for this.    

Usage:

	/path/to/starchip/starchip-circles.pl star_dirs.txt parameters.txt 

	or

	/path/to/starchip/starchip-circles.pl fastq_files.txt parameters.txt


	This command will generate scripts Step1.sh, Step2.sh, and optionally Step3.sh and Step4.sh  This is to allow users to optimize computing resources for a variety of computing environments.  
	starchip-circles is designed to work on a set of multiple samples (though it works fine on 1 sample). 
	star_dirs.txt should have 1 full pathway to star output per line: ie 
		/path/to/file1/star/
		/path/to/file2/star/
	
	Alternatively, if you want STARChip to perform STAR realignment, fastq_files.txt should have fastq files:
		/path/to/file1.1_R1.fastq,/path/to/file1.2_R1.fastq  /path/to/file1.1_R2.fastq,/path/to/file1.2_R2.fastq	
		/path/to/file2.1_R1.fastq,/path/to/file2.2_R1.fastq  /path/to/file2.1_R2.fastq,/path/to/file2.2_R2.fastq	

	An example parameter file is located at /path/to/starchip/paramfiles/starchip-circles.params and full details can be found in the manual.

Run Modes: 

	STARChip circRNA can be run on just STAR outputs, or it can work directly from fastq.  If running from fastq files, set runSTAR=true in your parameters file.  
	If runSTAR=true, STARChip will create scripts to align (Step1.sh), discover circRNA (Step2.sh), re-align (Step3.sh), and quantify/annotate (Step4.sh).
	If runSTAR=false, STARChip will create scripts to discover (Step1.sh), and quantify/annotate (Step2.sh).  If runSTAR=false, you must have run STAR to 
	create both Chimeric.out.junction and Chimeric.out.sam files.  Prior to STAR v2.6, these files are generated anytime Chimeric output is on. From STAR v2.6.0c
	forward, you must have run STAR with chimeric output on and the flag --chimOutType Junctions SeparateSAMold .  

Software Dependencies:
	
	R with the following packages: limma, edgeR, gplots
	bedtools (>=2.24.0)
	STAR (NOT version 2.6.0a or 2.6.0b)

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



##  Fusions  ##

Usage:

	/path/to/starchip/starchip-fusions.pl output_seed Chimeric.junction.out Parameters.txt
	
		This will create output_seed.summary and output_seed.summary.annotated. 

Softare Dependencies:

	SAMtools
	BEDtools (>= 2.24.0)
	mafft
	STAR (starchip doesn't call star, but you do need star output)

Files Needed:
	
	parameters file (examples included in /starchip/paramfiles )
	Fasta File and GTF File used to generate STAR index 
	repetative elements/antibody parts regions files in .bed format (included for hg19,hg38)
	Optional: known gene families and paralogs (included for human)
	You must have run STAR to create both Chimeric.out.junction and Chimeric.out.sam files.  Prior to STAR v2.6, these files are generated anytime Chimeric output is on. From STAR v2.6.0c forward, you must have run STAR with chimeric output on and the flag --chimOutType Junctions SeparateSAMold .

Output:

	you'll get a [output_seed].annotated and a [output_seed].annotated.summary file
	The summary file is more compact, while the annotated file has more information.  

Known Issues:

	Fusions are tricky.  This still has a high false positive rate, due to short overhangs and paralogues/pseudogenes. 
	A good second step is to BLAST/BLAT your fusion consensus sequence, especially to the newest genome assembly available. 


