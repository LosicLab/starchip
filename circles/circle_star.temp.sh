##usage: ./circle_star.sh readsMinimum minSubjects stardirs.txt splice
#for example with the structure: /data/star_alignments/star1, star2, star3, star4
# you would run ./circle_star.sh 5 5 /data/star_alignments/ splice

# calls Rscript

##Check input arguments:
if [ $# -ne 4 ]
  then
    echo "Wrong number of arguments supplied.  usage: ./circle_star.sh readsMinimum minSubjects stardirs.txt splice"
    echo "where stardirs.txt has one star output directory per line, splice if you want to perform splicing (takes time) or nosplice if not."
    exit 1
fi

#make a clean exit, if something goes wrong
trap "exit" INT TERM
trap "kill 0" EXIT

#How many subjects must express a cRNA species?
minSubjLimit=${2}
#How many reads for a subject to express a cRNA species? To look at multiple thresholds, modify the line below to something like: cutofflist=( 1 3 5 10 100)
cutofflist=${1}
#How many CPUs do we have?
cpus=8
#CPM thresholds
cpmcutoff=0
subjectcpm=0


echo "Circular RNA species must have at least "${cutofflist}" reads in at least "${minSubjLimit}" subjects/output files.  Using $cpus CPUs."
echo "For CPMs: Rscript must be callable.  Must be ${subjectcpm} subjects/outputs with ${cpmcutoff} Counts per million circular reads to count a given circular RNA"

# get the full path to this script
DIR=$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )

##cleanup any old files
for cutoff in "${cutofflist[@]}" ; do
	rm -f ${cutoff}
done
rm -f temp_subject_IDS.txt

#The while read loop takes in a file as input, with a star output directory on each line.  ie /path/to/star_align_file1/
#create joinstrands files.  these have the same cRNA concatenated to one strand. 
while read f ; do
	#define subject names/paths
	array=(${f//\// })
	length=${#array[@]}
	uniqIDindex=$((length-1))
	uniqID=${array[${uniqIDindex}]}
	Parent="$(dirname "$f")"
	echo $uniqID >> temp_subject_IDS.txt
	IDarray=("${IDarray[@]}" "${f}")
done < ${3}

#filter chimera to circular, sort, create joinstrands file. #chimera-score.pl checks the directory $f for "Chimeric.out.sam"  and "Chimeric.out.junction"
cat $3 | xargs --max-procs=${cpus} -I {} ${DIR}/chimera-score.pl {}/ ${DIR}/filter_circsv3.pl
#done < ${3}
wait  #wait for all joinstrands files to finish.

##Generate a cutoff file with all cRNA above the reads cutoff
for cutoff in "${cutofflist[@]}" ; do
	##now toss all candidates that pass a read support cutoff into a file named for the cutoff# (ie 10,20,50)
	cat temp_subject_IDS.txt | xargs --max-procs=${cpus} -I {} awk -v var="$cutoff" '{ if ($1 >= var) print $0 >> var }' joinstrands.{}
done

##  create the list of interesting cRNA to investigate (those cRNA with > reads cutoff and > indiv. cutoff)
for cutoff in "${cutofflist[@]}" ; do
	cut -f4,5,8 ${cutoff} | sort -k1,1 -k2,2n -k3,3n | uniq -c | awk -v var=${minSubjLimit} '{ if ($1 >= var) print $2,$3,$4,$1 }' > circs"${cutoff}"."${minSubjLimit}".investigate
	# sed -i '1itotalreads\t+strand\t-strand\tchrm\tpos\tstrand\tchrm\tpos\tstrand\tjxntype\toverlapL\toverlapR\tAvgScore\tLeftScore\tRightScore' ${cutoff}
done 
#generates circs"${cutoff}"."${minSubjLimit}".investigate

## Create a count matrix using joinstrands files + circs"${cutoff}"."${minSubjLimit}".investigate (separate file is called make_circ_table.sh)
echo "Creating a count matrix.  This may take a long time."
rm -f circ.headers
rm -f  row.names
for f in ${cutofflist[@]} ; do
        rm -f circ.headers.cutoff${f}.${minSubjLimit}
        rm -f circRNA.cutoff${f}.${minSubjLimit}
done
#create the row names.
while read circ ; do 
	circarray=(${circ// / })
	echo "chr"${circarray[0]}":"${circarray[1]}"-"${circarray[2]} >>row.names
done < circs${cutoff}.${minSubjLimit}.investigate
#cycle through the cutoffs
for cutoff in ${cutofflist[@]} ; do
	cat row.names > circRNA.${cutoff}reads.${minSubjLimit}ind.countmatrix
	#cycle throgh all subject files
        for f in ./joinstrands.* ; do
		#cycle through circles
		rm -f ${f}.temp
		while read circ ; do
			circarray=(${circ// / })
			ucscname="chr"${circarray[0]}":"${circarray[1]}"-"${circarray[2]}
			#searchstring="${circarray[0]}\t${circarray[1]}\t.\t${circarray[0]}\t${circarray[2]}"
		        #grephit=$(grep -P $searchstring $f )
                	grephit=$(fgrep $ucscname $f )
			if  [[ $? -eq  0 ]] ; then 
        		        echo $grephit |cut -f1 -d" " >>${f}.temp
                        else
                                echo "0" >> ${f}.temp
                        fi
                done <circs${cutoff}.${minSubjLimit}.investigate
		#add this subject's data to the countmatrix
		paste circRNA.${cutoff}reads.${minSubjLimit}ind.countmatrix  ${f}.temp > tempmatrix
        	mv tempmatrix circRNA.${cutoff}reads.${minSubjLimit}ind.countmatrix
		#add to the header file. 
		echo -n -e ${f}"\t" | sed 's/\.\/joinstrands\.//' >>circ.headers.cutoff${cutoff}.${minSubjLimit}
		rm ${f}.temp
	done
	#clean up the header, put it on top of the count matrix. 
	sed -i 's/\t$/\n/' circ.headers.cutoff${cutoff}.${minSubjLimit}
	cat circ.headers.cutoff${cutoff}.${minSubjLimit} circRNA.${cutoff}reads.${minSubjLimit}ind.countmatrix >tempmatrix
	mv tempmatrix circRNA.${cutoff}reads.${minSubjLimit}ind.countmatrix 
	rm circ.headers.cutoff${cutoff}.${minSubjLimit}
done 
wait ${!}
rm row.names

#generate counts per million. 
for f in ${cutofflist[@]} ; do
	Rscript ${DIR}/cpm.R circRNA.${cutoff}reads.${minSubjLimit}ind.countmatrix ${cpmcutoff} ${subjectcpm} &
done 
wait 



#splicing
if [ $4 == 'splice' ] ; then
	echo "splicing... "
	for cutoff in "${cutofflist[@]}" ; do
		#splice the candidates using the SJ.out.tab file for the individual. 
		cat temp_subject_IDS.txt | xargs --max-procs=${cpus} -I {} ${DIR}/circle_splice_individual.sh joinstrands.{} ${Parent}/{} ${cutoff}
	done
	wait

	##Calculate Statistics on the candidates
	for f in ${cutofflist[@]} ; do 
		${DIR}/candidate_stats.pl circs${f}.${minSubjLimit}.investigate ${f}.spliced &
	done
	wait
fi
#generates .consensus and .allvariants files.  



