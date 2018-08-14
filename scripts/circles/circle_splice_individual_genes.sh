#!/bin/bash

#usage: ./circle_splice_individual.sh star_alignment_dir/ readscutoff star_out_prefix uniqID  >output_file

#the circular candidates should be in the format output by filter_circs.pl
# ie  reads_support chrom finish strand chrom start etc etc. 

#echo -e "subject\ttotalreads\t+strand\t-strand\tchrm\tpos\tstrand\tchrm\tpos\tstrand\tjxntype\toverlapL\toverlapR\talignmentscoreMedian\tscoreL\tscoreR\tGenomicSize\tSplicedSize\tLeftBorder\tRightBorder\tLeftEnvelope\tRightEnvelope\tEnvREads\tEnvelopingStrand\tEnvelopingJxnType\tExons\tExonStarts\tExonsizes"
star_aligns=$1
cutoff=$2
prefix=$3
IDstepsback=$4
star_aligns=${star_aligns}/
star_aligns=`echo $star_aligns | sed 's/\/\/$/\//'`
array=(${star_aligns//\// })
length=${#array[@]}
uniqIDindex=$((length-$IDstepsback))
uniqID=${array[${uniqIDindex}]}
joinstrands=rawdata/backsplices.${uniqID}
outfile=`echo $joinstrands |sed 's/backsplices.//'`
outfile=${outfile}.${cutoff}.spliced
if [ $prefix == "NoPrefix123456789" ] ; then
        prefix="";
fi
sj_out=${star_aligns}${prefix}SJ.out.tab
joinstrands=rawdata/circRNA.gene.coords
# line is a circular RNA candidate
if [ -f $joinstrands ] ; then
	while read line ; do
   	linearray=(${line// / })
  	#echo ${linearray[0]}
   	#if  (( ${linearray[0]} >= $cutoff )) ; then
		rm -f ${uniqID}.tempfile
		rm -f ${uniqID}.border-jxns-right
		rm -f ${uniqID}.border-jxns-left
		rm -f ${uniqID}.enveloping-reads
		rm -f ${uniqID}.sjmax 
		#put all splice sites within the circle into a tempfile (SJ.out.tab must be in the parent directory)
		#also put any border splices into two files
		awk -v id=$uniqID -v chrom="${linearray[3]}" -v start="${linearray[4]}" -v finish="${linearray[7]}" -v reads="${linearray[0]}" ' {
 			if ($1 == chrom && $2 < start && $3 > finish && $7>=(reads*0.6) ) print $0,start,finish >>id".enveloping-reads"
			if ($1 == chrom && ( (($2 - finish) <=5 ) && (($2 - finish) >=-5 )) ) print $0,start,finish >> id".border-jxns-right"
			if ($1 == chrom && ( (($3 - start) <=5 ) && (($3 - start) >=-5 )) ) print $0,start,finish >> id".border-jxns-left"
			if ($1 == chrom && $2 > start && $3 < finish && $7>=(reads*0.6) ) print $0,start,finish >> id".tempfile"
			if ($1 == chrom && $2 >= start && $3 <= finish ) print $0,start,finish,"internal" >>id".sjmax"
			if ($1 == chrom && $2 <= start && $3 >= start && $3 <= finish ) print $0,start,finish,"splice-left" >>id".sjmax"
			if ($1 == chrom && $3 >= finish && $2  >= start && $2 <= finish ) print $0,start,finish,"splice-right" >>id".sjmax"
			} ' ${sj_out}
		spliceend=1
		spliceOutSize=0
		blockCount=1
		unset blockStart
		unset blockSize
		blockStart[1]=0
		let "blockSize[1] = ${linearray[7]} - ${linearray[4]}"
		#if a splice in the circle was found
		if [ -f ${uniqID}.tempfile ] ; then
			while read templine; do
				templinearray=(${templine// / }) ##templinearray: 0:chromosome 1:start  2:stop  3:strand 4:jxn type 5:annotation(0=un 1=anno) 6:unique spanning 7:multimapping spanning  8:max overhang
				#if we have overlapping splice sites, select the one with read count closest to crna read count
				if (( ${templinearray[1]} <= $spliceend )) ; then
					#compare the two lines
					if (( $topreadcount < ${templinearray[6]} )) ; then
        					topreadcount=${templinearray[6]}
						#remove the old intron, insert the new intron
		                        	let "spliceOutSize -= intron"
						let "intron = ${templinearray[2]} - ${templinearray[1]}"
						let "spliceOutSize += intron"
						spliceend=${templinearray[2]}
						#update the exon information
						let "blockStart[$blockCount] = ${templinearray[2]} - ${linearray[4]}"
						let "blockCountless1 = $blockCount - 1"
						let "blockSize[$blockCountless1] = ${templinearray[1]} - ${blockStart[$blockCountless1]} - ${linearray[4]}"
						#QC check
						#echo "circreads:"${linearray[0]} "intronreads:"${templinearray[6]} "start:"${templinearray[1]} "stop:"${templinearray[2]} "isize:"$intron "splicesize:"$spliceOutSize 
					fi
				#if it's not an overlapping splice site, add up the intron lengths
				else
					#add intron information
					spliceend=${templinearray[2]}
					let "intron = ${templinearray[2]} - ${templinearray[1]}"
					let "spliceOutSize += intron"
					#add exon information
					let "blockCount = $blockCount + 1"
					let "blockStart[$blockCount] = ${templinearray[2]} - ${linearray[4]}"
					let "blockCountless1 = $blockCount - 1"
                                	let "blockSize[$blockCountless1] = ${templinearray[1]} - ${blockStart[$blockCountless1]} - ${linearray[4]}"
					topreadcount=${templinearray[6]}
				fi
			done < ${uniqID}.tempfile
			let "blockSize[$blockCount] = ${linearray[7]} - ${linearray[4]} - ${blockStart[$blockCount]}"
		fi
		#above I created border-jxns-right which contains lines from SJ.out.tab that contain splices from the right side of my circle to upstream.  here I find the smallest of these introns
		declare -i smallestRight
		smallestRight=999999999
		rightIntron=.
		if [ -f ${uniqID}.border-jxns-right ] ; then
			while read rightjunc ; do
				rightjuncarray=(${rightjunc// / }) #0:chromosome 1:start  2:stop  3:strand 4:jxn type 5:annotation(0=un 1=anno) 6:unique spanning 7:multimapping spanning  8:max overhang
				#echo "${rightjuncarray[2]} ${linearray[4]} $smallestRight"
				if (( "${rightjuncarray[2]}" < "$smallestRight" )) ; then
					smallestRight="${rightjuncarray[2]}"
					let "rightIntron = $smallestRight - ${linearray[7]}"
					#echo "$smallestRight - ${linearray[4]} eq $rightIntron"
					#echo "$rightIntron"
				fi
			done < ${uniqID}.border-jxns-right
		fi
		#figure out the smallest left border intron
		declare -i largestLeft
		largestLeft=0
		leftIntron=.
		if [ -f ${uniqID}.border-jxns-left ] ; then
			while read leftjunc ; do
				leftjuncarray=(${leftjunc// / }) #0:chromosome 1:start  2:stop  3:strand 4:jxn type 5:annotation(0=un 1=anno) 6:unique spanning 7:multimapping spanning  8:max overhang
				#echo "${leftjuncarray[1]} ${linearray[7]} $largestLeft"
				if (( "${leftjuncarray[1]}" > "$largestLeft" )) ; then
					largestLeft="${leftjuncarray[1]}"
					let "leftIntron = ${linearray[4]} - largestLeft"
					#echo "$largestLeft  - ${linearray[7]} eq $leftIntron"
					#echo "$leftIntron"
				fi
			done < ${uniqID}.border-jxns-left
		fi
		#find the enveloping read with the most read support
		leftEnvelope=.
		rightEnvelope=.
		declare -i envSupport
		envSupport=0
		envStrand=.
		envType=.
		if [ -f ${uniqID}.enveloping-reads ] ; then
			while read envreads ; do
				envarray=(${envreads// / }) ##envarray: 0:chromosome 1:start  2:stop  3:strand 4:jxn type 5:annotation(0=un 1=anno) 6:unique spanning 7:multimapping spanning  8:max overhang
				if  (( ${envarray[6]} > $envSupport )) ; then
					leftEnvelope=${envarray[1]}
					rightEnvelope=${envarray[2]}
					envSupport=${envarray[6]}
					envStrand=${envarray[3]}
					envType=${envarray[4]}
				fi
			done < ${uniqID}.enveloping-reads
		fi
		#find the SJmax
		splicetype=.
		declare -i sjmax
		sjmax=0;
		if [ -f ${uniqID}.sjmax ] ; then
			while read sjreads ; do
				sjarray=(${sjreads// / }) ##sjarray: 0:chromosome 1:start  2:stop  3:strand 4:jxn type 5:annotation(0=un 1=anno) 6:unique spanning 7:multimapping spanning  8:max overhang
				if (( ${sjarray[6]} > $sjmax )) ; then
					sjmax=${sjarray[6]}
					splicetype=${sjarray[11]}
				fi
			done < ${uniqID}.sjmax
		fi

		#calculate the total size of the crna
		let "circleSizeGenomic = ${linearray[7]} - ${linearray[4]}"
		let "circleSizeSpliced = $circleSizeGenomic - $spliceOutSize"
		echo -ne "$uniqID\t$line\t$circleSizeGenomic\t$circleSizeSpliced\t$leftIntron\t$rightIntron\t$leftEnvelope\t$rightEnvelope\t$envSupport\t$envStrand\t$envType\t$sjmax\t$splicetype\t"
		(IFS=,; echo -e "${blockCount}\t${blockStart[*]}\t${blockSize[*]}")
		#echo -ne "$uniqID\n$line\n$circleSizeGenomic\n$circleSizeSpliced\n$leftIntron\n$rightIntron\n$leftEnvelope\n$rightEnvelope\n$envSupport\n$envStrand\n$envType\n$sjmax\n$splicetype\n"

		rm -f ${uniqID}.tempfile
		rm -f ${uniqID}.border-jxns-right
		rm -f ${uniqID}.border-jxns-left
		rm -f ${uniqID}.enveloping-reads
		rm -f ${uniqID}.sjmax
  	#fi
	done < $joinstrands > $outfile
	sed -i '1isubject\ttotalreads\t+strand\t-strand\tchrm\tpos\tstrand\tchrm\tpos\tstrand\tjxntype\toverlapL\toverlapR\talignmentscoreMedian\tscoreL\tscoreR\tID\tGenomicSize\tSplicedSize\tLeftBorder\tRightBorder\tLeftEnvelope\tRightEnvelope\tEnvREads\tEnvelopingStrand\tEnvelopingJxnType\tLinearMaxSplice\tLinearSpliceMaxType\tExons\tExonStarts\tExonsizes' $outfile
fi


