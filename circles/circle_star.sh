##usage: ./circle_star.sh readsMinimum minSubjects stardirs.txt splice
#for example with the structure: /data/star_alignments/star1, star2, star3, star4
# you would run ./circle_star.sh 5 5 /data/star_alignments/ splice

# calls Rscript

##list of arguments required:
#  1: readsCutoff  2: minSubjectsLimit  3:star_dir_list.txt  4:do_splice(true/false)  5:cpus  6: cpmCutoff  7: subjectCPMcutoff
#  8: annotate(true/false)  9:refbed(if annotate==true) 10: prefix on your star output 11: steps back to your data ID. 

if [ "$#" -ne 11 ]; then
    echo "Illegal number of parameters, please check your paramters file. Note that this script should be run by parent script starchimp-circles.pl"
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
	echo $uniqID >> rawdata/temp_subject_IDS.txt
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
        cat rawdata/temp_subject_IDS.txt | xargs --max-procs=${cpus} -I {} awk -v var="$cutoff" -v var2="rawdata/cRNA.cutoff.${cutoff}" '{ if ($1 >= var) print $0 >> var2 }' rawdata/backsplices.{}
done


##  create the list of interesting cRNA to investigate (those cRNA with > reads cutoff and > indiv. cutoff)
for cutoff in "${cutofflist[@]}" ; do
	cut -f4,5,8 rawdata/cRNA.cutoff.${cutoff} | sort -k1,1 -k2,2n -k3,3n | uniq -c | awk -v var=${minSubjLimit} '{ if ($1 >= var) print $2,$3,$4,$1 }' > rawdata/circs"${cutoff}"."${minSubjLimit}".investigate
	# sed -i '1itotalreads\t+strand\t-strand\tchrm\tpos\tstrand\tchrm\tpos\tstrand\tjxntype\toverlapL\toverlapR\tAvgScore\tLeftScore\tRightScore' ${cutoff}
done 
#generates circs"${cutoff}"."${minSubjLimit}".investigate

## Create a count matrix using backsplice files + circs"${cutoff}"."${minSubjLimit}".investigate (separate file is called make_circ_table.sh)
echo "Creating a count matrix.  This may take a long time."
for f in ${cutofflist[@]} ; do
        rm -f circ.headers.cutoff${f}.${minSubjLimit}
        rm -f circRNA.cutoff${f}.${minSubjLimit}
	rm -f row.names.${cutoff}
done
#create the row names.
for cutoff in ${cutofflist[@]} ; do
	while read circ ; do 
		circarray=(${circ// / })
		echo ${circarray[0]}":"${circarray[1]}"-"${circarray[2]} >> row.names.${cutoff}
	done < rawdata/circs${cutoff}.${minSubjLimit}.investigate   
done
#cycle through the cutoffs
for cutoff in ${cutofflist[@]} ; do
	cat row.names.${cutoff} > circRNA.${cutoff}reads.${minSubjLimit}ind.countmatrix
	#cycle throgh all subject files
        for f in rawdata/backsplices.* ; do
		#cycle through circles
		rm -f ${f}.temp
		while read circ ; do
			circarray=(${circ// / })
			ucscname=${circarray[0]}":"${circarray[1]}"-"${circarray[2]}
			#searchstring="${circarray[0]}\t${circarray[1]}\t.\t${circarray[0]}\t${circarray[2]}"
		        #grephit=$(grep -P $searchstring $f )
                	grephit=$(fgrep $ucscname $f )
			if  [[ $? -eq  0 ]] ; then 
        		        echo $grephit |cut -f1 -d" " >>${f}.temp
                        else
                                echo "0" >> ${f}.temp
                        fi
                done <rawdata/circs${cutoff}.${minSubjLimit}.investigate
		#add this subject's data to the countmatrix
		paste circRNA.${cutoff}reads.${minSubjLimit}ind.countmatrix  ${f}.temp > tempmatrix
        	mv tempmatrix circRNA.${cutoff}reads.${minSubjLimit}ind.countmatrix
		#add to the header file. 
		echo -n -e ${f}"\t" | sed 's/rawdata\/backsplices\.//' >>circ.headers.cutoff${cutoff}.${minSubjLimit}
		rm ${f}.temp
	done
	#clean up the header, put it on top of the count matrix. 
	sed -i 's/\t$/\n/' circ.headers.cutoff${cutoff}.${minSubjLimit}
	cat circ.headers.cutoff${cutoff}.${minSubjLimit} circRNA.${cutoff}reads.${minSubjLimit}ind.countmatrix >tempmatrix
	mv tempmatrix circRNA.${cutoff}reads.${minSubjLimit}ind.countmatrix 
	rm -f circ.headers.cutoff${cutoff}.${minSubjLimit} row.names.${cutoff}
done
##wait ${!}

#generate counts per million. 
echo "generating counts per million, and transforming your data with voom"
for f in ${cutofflist[@]} ; do
	Rscript ${DIR}/cpm.R circRNA.${cutoff}reads.${minSubjLimit}ind.countmatrix ${cpmcutoff} ${subjectcpm} &
done 
wait 



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
	wait
	##generate sj max matrix
	echo "creating matrix of maximum linear splices"
	for cutoff in ${cutofflist[@]} ; do
        	rm -f circ.headers.cutoff${cutoff}.${minSubjLimit}
        	cut -f1 circRNA.${cutoff}reads.${minSubjLimit}ind.countmatrix |tail -n +2 > rawdata/linear.${cutoff}reads.${minSubjLimit}ind.sjmax
        	cp rawdata/linear.${cutoff}reads.${minSubjLimit}ind.sjmax rawdata/linear.${cutoff}reads.${minSubjLimit}ind.splicetype
        	for f in rawdata/*spliced ; do 
                	rm -f ${f}.temp ${f}.temp2
                	while read circ ; do 
                    	circarray=(${circ// / })
                    	ucscname=${circarray[0]}":"${circarray[1]}"-"${circarray[2]}
                    	grephit=$(fgrep $ucscname $f )
                    	if  [[ $? -eq  0 ]] ; then 
                            	echo $grephit |cut -f27 -d" " >>${f}.temp
                            	echo $grephit |cut -f28 -d" " >>${f}.temp2
                 	else
                            	echo "0" >> ${f}.temp
                            	echo "." >> ${f}.temp2
                    	fi
                	done <rawdata/circs${cutoff}.${minSubjLimit}.investigate
                	#paste columns to full matrix
                	paste rawdata/linear.${cutoff}reads.${minSubjLimit}ind.sjmax ${f}.temp > tempmatrix
                	paste rawdata/linear.${cutoff}reads.${minSubjLimit}ind.splicetype ${f}.temp2 >tempmatrix2
                	mv tempmatrix rawdata/linear.${cutoff}reads.${minSubjLimit}ind.sjmax
                	mv tempmatrix2 rawdata/linear.${cutoff}reads.${minSubjLimit}ind.splicetype
                	#add to the header file.
                	echo -n -e ${f}"\t" | sed 's/rawdata\///' | sed "s/\.${cutoff}.spliced//"  >>circ.headers.cutoff${cutoff}.${minSubjLimit}
                	rm ${f}.temp ${f}.temp2
        	done
        	#clean up the header, put it on top of the count matrix. 
        	sed -i 's/\t$/\n/' circ.headers.cutoff${cutoff}.${minSubjLimit}
        	cat circ.headers.cutoff${cutoff}.${minSubjLimit} rawdata/linear.${cutoff}reads.${minSubjLimit}ind.sjmax >tempmatrix
        	cat circ.headers.cutoff${cutoff}.${minSubjLimit} rawdata/linear.${cutoff}reads.${minSubjLimit}ind.splicetype >tempmatrix2
        	mv tempmatrix rawdata/linear.${cutoff}reads.${minSubjLimit}ind.sjmax
        	mv tempmatrix2 rawdata/linear.${cutoff}reads.${minSubjLimit}ind.splicetype
        	rm circ.headers.cutoff${cutoff}.${minSubjLimit} 
	done
fi
#generates .consensus and .allvariants files, as well as .sjmax and .splicetype matrixes.  

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
