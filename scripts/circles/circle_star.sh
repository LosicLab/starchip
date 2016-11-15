#!/bin/bash

# calls Rscript (R required)

##list of arguments required:
#  1: readsCutoff  2: minSubjectsLimit  3:star_dir_list.txt  4:do_splice(true/false)  5:cpus  6: cpmCutoff  7: subjectCPMcutoff
#  8: annotate(true/false)  9:refbed(if annotate==true) 10: prefix on your star output 11: steps back to your data ID. 

if [ "$#" -ne 11 ]; then
    echo "Illegal number of parameters, please check your paramters file. Note that this script should be run by parent script starchip-circles.pl"
    exit 1; 
fi


#make a clean exit, if something goes wrong
#trap "exit" INT TERM
#trap "kill 0" EXIT

#How many subjects must express a cRNA species?
minSubjLimit=${2}
#How many reads for a subject to express a cRNA species? To look at multiple thresholds, parameter should look like: cutofflist=( 1 3 5 10 100)
cutofflist=${1}
#How many CPUs do we have?
cpus=${5}
#CPM thresholds
cpmcutoff=${6}
subjectcpm=${7}
annotate=${8}
refbed=${9}
starprefix=${10}
IDstepsback=${11}

echo "Circular RNA species must have at least "${cutofflist}" reads in at least "${minSubjLimit}" subjects/output files.  Using $cpus CPUs."
echo "For CPMs: Rscript must be callable.  Must be ${subjectcpm} subjects/outputs with ${cpmcutoff} Counts per million circular reads to count a given circular RNA"

# get the full path to this script
DIR=$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )

# create a junk directory 
mkdir -p rawdata

##cleanup any old files
for cutoff in "${cutofflist[@]}" ; do
	rm -f rawdata/cRNA.cutoff.${cutoff}
done
rm -f rawdata/temp_subject_IDS.txt

#The while read loop takes in a file as input, with a star output directory on each line.  ie /path/to/star_align_file1/
#create backsplices files.  these have the same cRNA concatenated to one strand. 
while read f ; do
	#define subject names/paths
	f=${f}/
	f=`echo $f | sed 's/\/\/$/\//'`
	array=(${f//\// })
	length=${#array[@]}
	uniqIDindex=$((length-$IDstepsback))
	uniqID=${array[${uniqIDindex}]}
	Parent="$(dirname "$f")"
	#echo $uniqID >> rawdata/temp_subject_IDS.txt
	IDarray=("${IDarray[@]}" "${f}")
done < ${3}

#filter chimera to circular, sort, create backsplice file. #chimera-score.pl checks the directory $f for "Chimeric.out.sam"  and "Chimeric.out.junction"
echo "Filtering out chimeric reads that appear circular"
cat $3 | sed 's/$/\//' |sed 's/\/\/$/\//' | xargs --max-procs=${cpus} -I {} ${DIR}/chimera-score.pl {} $starprefix $IDstepsback ${DIR}/filter_circsv3.pl
#done < ${3}
wait  #wait for all backsplices files to finish.
echo "Filtering circular reads based on assigned thesholds of ${cutofflist[@]} reads in ${minSubjLimit} individuals"
##Generate a cutoff file with all cRNA above the reads cutoff
for cutoff in "${cutofflist[@]}" ; do
        ##now toss all candidates that pass a read support cutoff into a file named for the cutoff# (ie 10,20,50)
        ls rawdata/backsplices* | xargs --max-procs=${cpus} -I {} awk -v var="$cutoff" -v var2="rawdata/cRNA.cutoff.${cutoff}" '{ if ($1 >= var) print $0 >> var2 }' {}
done


##  create the list of interesting cRNA to investigate (those cRNA with > reads cutoff and > indiv. cutoff)
for cutoff in "${cutofflist[@]}" ; do
	cut -f4,5,8 rawdata/cRNA.cutoff.${cutoff} | sort -k1,1 -k2,2n -k3,3n | uniq -c | awk -v var=${minSubjLimit} '{ if ($1 >= var) print $2":"$3"-"$4,$1 }' OFS="\t" \
	|sort -k2,2nr | ${DIR}/merge_close_crna_1file.pl | sed 's/[: \t-]/ /g' > rawdata/circs"${cutoff}"."${minSubjLimit}".investigate
	# sed -i '1itotalreads\t+strand\t-strand\tchrm\tpos\tstrand\tchrm\tpos\tstrand\tjxntype\toverlapL\toverlapR\tAvgScore\tLeftScore\tRightScore' ${cutoff}
done 
#generates circs"${cutoff}"."${minSubjLimit}".investigate

##
## Create a count matrix using backsplice files + circs"${cutoff}"."${minSubjLimit}".investigate
##
echo "Creating a count matrix.  This may take a long time."
for cutoff in ${cutofflist[@]} ; do
	#clean up any potential old files
	rm -f row.names.${cutoff}
	#create the rown.names files
	echo "cRNA" >> row.names.${cutoff}
	while read circ ; do
		circarray=(${circ// / })
		echo ${circarray[0]}":"${circarray[1]}"-"${circarray[2]} >> row.names.${cutoff}
	done < rawdata/circs${cutoff}.${minSubjLimit}.investigate
	#cycle throgh all subject files, create columns
	rm -f rawdata/backsplices.*grepcircles.temp
	ls rawdata/backsplices.* | xargs --max-procs=${cpus} -I {} ${DIR}/grep-circles.pl rawdata/circs${cutoff}.${minSubjLimit}.investigate {} count
	#assemble the matrix by pasting columns
	paste row.names.${cutoff} rawdata/backsplices.*grepcircles.temp > circRNA.${cutoff}reads.${minSubjLimit}ind.countmatrix
	#clean up: remove columns
	rm -f row.names.${cutoff} rawdata/backsplices.*grepcircles.temp
done
wait


#Annotate Genes
if [ ${annotate} == 'true' ] ; then
	echo "annotating your cRNA"
	for f in ${cutofflist[@]} ; do
		cut -f1 circRNA.${cutoff}reads.${minSubjLimit}ind.countmatrix > circRNA.${cutoff}reads.${minSubjLimit}ind
		${DIR}/cRNA_coordinates2genes.sh circRNA.${cutoff}reads.${minSubjLimit}ind $refbed
		rm circRNA.${cutoff}reads.${minSubjLimit}ind
		#summarize the data
		tail -n +2 circRNA.${cutoff}reads.${minSubjLimit}ind.annotated | ${DIR}/summarize_geneinfo.pl > circRNA.${cutoff}reads.${minSubjLimit}ind.genes
	done
fi

#generate counts per million. 
echo "generating counts per million, transforming your data with voom, and generating plots"
for f in ${cutofflist[@]} ; do
	Rscript ${DIR}/cpm2.R circRNA.${cutoff}reads.${minSubjLimit}ind.countmatrix ${cpmcutoff} ${subjectcpm} circRNA.${cutoff}reads.${minSubjLimit}ind.genes &
done


#splicing
if [ ${4} == 'true' ] ; then
	echo "Now splicing your cRNA "
	for cutoff in "${cutofflist[@]}" ; do
		#splice the candidates using the SJ.out.tab file for the individual. 
		cat ${3} | xargs --max-procs=${cpus} -I {} ${DIR}/circle_splice_individual.sh {} ${cutoff} ${starprefix} ${IDstepsback}
	done
	wait

	##Calculate Statistics on the candidates
	for f in ${cutofflist[@]} ; do 
		${DIR}/candidate_stats.pl rawdata/circs${f}.${minSubjLimit}.investigate ${f}.spliced &
	done
	##generate sj max matrix
	echo "creating matrix of maximum linear splices"
	for cutoff in ${cutofflist[@]} ; do
		#create the rownames
        	cut -f1 circRNA.${cutoff}reads.${minSubjLimit}ind.countmatrix > row.names.${cutoff}
		#create all of the columns as rawdata/*spliced.temp and rawdata/*spliced.temp2 
		ls rawdata/*spliced | xargs --max-procs=${cpus} -I {} ${DIR}/grep-circles.pl rawdata/circs${cutoff}.${minSubjLimit}.investigate {} linear
		# paste all the columns together
		paste row.names.${cutoff} rawdata/*spliced.grepcircles.temp > rawdata/linear.${cutoff}reads.${minSubjLimit}ind.sjmax
		paste row.names.${cutoff} rawdata/*spliced.grepcircles.temp2 > rawdata/linear.${cutoff}reads.${minSubjLimit}ind.splicetype
		#remove all the column files
        	rm rawdata/*spliced.grepcircles.temp rawdata/*spliced.grepcircles.temp2 row.names.${cutoff}
	done
	wait
fi
#generates .consensus and .allvariants files, as well as .sjmax and .splicetype matrixes.

