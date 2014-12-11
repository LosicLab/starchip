##usage: ./circle_star.sh readsMinimum minSubjects stardirs.txt
#for example with the structure: /data/star_alignments/star1, star2, star3, star4
# you would run ./circle_star.sh /data/star_alignments/

##Check input arguments:
if [ $# -ne 3 ]
  then
    echo "Wrong number of arguments supplied.  usage: ./circle_star.sh readsMinimum minSubjects stardirs.txt "
    exit 1
fi

#make a clean exit, if something goes wrong
trap "exit" INT TERM
trap "kill 0" EXIT

#How many subjects must express a cRNA species?
minSubjLimit=${2}
#How many reads for a subject to express a cRNA species? To look at multiple thresholds, modify the line below to something like: cutofflist=( 1 3 5 10 100)
cutofflist=${1}
echo "Circular RNA species must have at least "${cutofflist}" reads in at least "${minSubjLimit}" subjects/output files"

# get the full path to this script
DIR=$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )

##cleanup any old files
for cutoff in "${cutofflist[@]}" ; do
	rm -f ${cutoff}
done

#Two ways to run.  The for loop takes as input just the parent directory of your star directories.  
#	The while read loop takes in a file as input, with a star output directory on each line.  ie /path/to/star_align_file1/
#this for loops through all the subjects by taking the directory you give it, then looking in all subdirectories.  

#create joinstrands
#for f in ${1}*/ ; do
while read f ; do
	#define subject names/paths
	array=(${f//\// })
	length=${#array[@]}
	uniqIDindex=$((length-1))
	uniqID=${array[${uniqIDindex}]}
	#filter chimera to circular, sort, create a bed file for each person, create joinstrands file #chimera-score.pl checks the directory $f for "Chimeric.out.sam"  and "Chimeric.out.junction"
	${DIR}/chimera-score.pl $f/ |	awk '$1!="chrM" && $1!="MT" && $8<=5 && $9<=5 { print $1,$2,$3,$4,$5,$6,$7,$8,$9,$15}' |
	awk '{ if ($3 =="+") print $1,$2,$3,$4,$5,$6,$7,$8,$9,$10 ;else print $1,$5,$3,$4,$2,$6,$7,$8,$9,$10 }' | sort -k1,1 -k2,2n -k5,5n -k6,6 |
	${DIR}/filter_circsv3.pl ${uniqID} &
	IDarray=("${IDarray[@]}" "${f}")
done < ${3}
#done
wait  #wait for all joinstrands files to finish.

##create  .spliced files, $cutoff files.  the spliced process runs in bg while the rest of the script carries on. 
for f in "${IDarray[@]}" ; do 
	array=(${f//\// })
        length=${#array[@]}
        uniqIDindex=$((length-1))
        uniqID=${array[${uniqIDindex}]}
	for cutoff in "${cutofflist[@]}" ; do
		##now toss all candidates that pass a read support cutoff into a file named for the cutoff# (ie 10,20,50)
		awk -v var="$cutoff" '{	if ($1 >= var) print $0 >> var }' joinstrands.${uniqID}
		#splice the candidates using the SJ.out.tab file for the individual.  
		awk -v var="$cutoff" '{if ($1 >= var) print $0 }' joinstrands.${uniqID} | ${DIR}/circle_splice_individual.sh ${uniqID} $f > ${uniqID}.${cutoff}.spliced &
	done
done 

##  create the list of interesting cRNA to investigate.  
for cutoff in "${cutofflist[@]}" ; do
	cut -f4,5,8 ${cutoff} | sort -k1,1 -k2,2n -k3,3n | uniq -c | awk -v var=${minSubjLimit} '{ if ($1 >= var) print $2,$3,$4,$1 }' > circs"${cutoff}"."${minSubjLimit}".investigate
	# sed -i '1itotalreads\t+strand\t-strand\tchrm\tpos\tstrand\tchrm\tpos\tstrand\tjxntype\toverlapL\toverlapR\tAvgScore\tLeftScore\tRightScore' ${cutoff}
done 
#generates circs"${cutoff}"."${minSubjLimit}".investigate


## Create a count matrix using joinstrands files + circs"${cutoff}"."${minSubjLimit}".investigate (separate file is called make_circ_table.sh)
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
			#ucscname="chr"${circarray[0]}":"${circarray[1]}"-"${circarray[2]}
			searchstring="${circarray[0]}\t${circarray[1]}\t.\t${circarray[0]}\t${circarray[2]}"
		        grephit=$(grep -P $searchstring $f )
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
done &
wait ${!}
for f in ${cutofflist[@]} ; do
	Rscript ${DIR}/cpm.R circRNA.${cutoff}reads.${minSubjLimit}ind.countmatrix 0 0 &
done 

wait 
rm row.names
##Calculate Statistics on the candidates
for f in ${cutofflist[@]} ; do 
	${DIR}/candidate_stats.pl circs${f}.${minSubjLimit}.investigate ${f}.spliced &
done
wait
#generates .consensus and .allvariants files.  



