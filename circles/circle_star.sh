##usage: ./circle_star.sh /path/to/star_outputs/
#for example with the structure: /data/star_alignments/star1, star2, star3, star4
# you would run ./circle_star.sh /data/star_alignments/

#How many subjects must express a cRNA species?
minSubjLimit=1
#How many reads for a subject to express a cRNA species? (multiple values separated by a space are OK)
cutofflist=(5)

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
	/hpc/users/akersn01/scripts/circRNA/chimera-score.pl $f |	awk '$1!="chrM" && $1!="MT" && $8<=5 && $9<=5 { print $1,$2,$3,$4,$5,$6,$7,$8,$9,$15}' |
	awk '{ if ($3 =="+") print $1,$2,$3,$4,$5,$6,$7,$8,$9,$10 ;else print $1,$5,$3,$4,$2,$6,$7,$8,$9,$10 }' | sort -k1,1 -k2,2n -k5,5n -k6,6 |
	/hpc/users/akersn01/scripts/circRNA/filter_circsv3.pl ${uniqID} &
	IDarray=("${IDarray[@]}" "${f}")
done < ${1}
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
		awk -v var="$cutoff" '{if ($1 >= var) print $0 }' joinstrands.${uniqID} | /hpc/users/akersn01/scripts/circRNA/circle_splice_individual.sh ${uniqID} $f > ${uniqID}.${cutoff}.spliced &
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
#create the row names column
for f in ./joinstrands.* ; do
        echo $f | sed 's/\.\/joinstrands.//' >> row.names
done
#cycle through the cutoffs
for cutoff in ${cutofflist[@]} ; do
	cat row.names > circRNA.${cutoff}reads.${minSubjLimit}ind.countmatrix
	#cycle through all circRNAs
	while read circ ; do
        	circarray=(${circ// / })
	        newfile=(${circ// /-})
        	rm -f $newfile
       		#newfile is a shitty name, but it means a given circRNA, ie 1-3000-3500
        	#searchstring: chrm position  strand chrm position 
        	searchstring="${circarray[0]}\t${circarray[1]}\t.\t${circarray[0]}\t${circarray[2]}"
        	#searchstring="${circarray[0]}\t${circarray[2]}\t.\t${circarray[0]}\t${circarray[1]}"
        	#cycle throgh all subject files
        	for f in ./joinstrands.* ; do
        	        grephit=$(grep -P $searchstring $f )
                	        if  [[ $? -eq  0 ]] ; then 
                        	        echo $grephit |cut -f1 -d" " >>$newfile
                        	else
                            		echo "0" >> $newfile
                        	fi
       		done
        	#now $newfile (named by the position of hte circRNA) has the spanning reads for each subject
        	paste circRNA.${cutoff}reads.${minSubjLimit}ind.countmatrix  $newfile > tempmatrix
        	mv tempmatrix circRNA.${cutoff}reads.${minSubjLimit}ind.countmatrix
        	#create a line file of headers
        	echo -n $newfile" " >>circ.headers.cutoff${cutoff}.${minSubjLimit}
        	rm -f $newfile
	done < circs${cutoff}.${minSubjLimit}.investigate
	sed -i 's/ $/\n/' circ.headers.cutoff${cutoff}.${minSubjLimit}
	cat circ.headers.cutoff${cutoff}.${minSubjLimit} circRNA.${cutoff}reads.${minSubjLimit}ind.countmatrix >tempmatrix
	mv tempmatrix circRNA.${cutoff}reads.${minSubjLimit}ind.countmatrix 
	rm circ.headers.cutoff${cutoff}.${minSubjLimit}
done &
rm row.names
wait 
##Calculate Statistics on the candidates
for f in ${cutofflist[@]} ; do 
	/hpc/users/akersn01/scripts/circRNA/candidate_stats.pl circs${f}.${minSubjLimit}.investigate ${f}.spliced &
done
wait
#generates .consensus and .allvariants files.  



