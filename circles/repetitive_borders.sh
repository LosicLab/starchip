##NOTE--this is untested.  I think it will work but I haven't actually run it yet.
#usage: repetitive_borders.sh A.bed B.bed
##Starting with 2 bed files: 
	#A: a bed file of candidate splices from your study 
	#B: a bed file of reference repeats
##create a new file, similar to A.bed, with columns for classes of reptetive sequences found. 
#to create A.bed form candidates.investigate: awk '{print "chr"$1,$2,$3}' candidates.investigate | tr ' ' "\t" >Afile.bed

#adjustable variables:
windowsize=500

#use bedtools window to find all lines in b which overlap in a window around each line in A.

module load BEDTools/2.19.1
bedtools window -w ${windowsize} -a $1 -b $2 | awk '{ if ($6 <= $2) print $0 ; else if ($5 >=$3) print $0 }' > ${1}.window${windowsize}

#use matched_comp_repsv2.pl to join multiple hits from B in a single A line and look for complimentary Alu 
~/scripts/circRNA/matched_comp_repsv2.pl ${1}.window${windowsize} > ${1}.window${windowsize}.complimentary

echo ${1}.window${windowsize}.complimentary
for f in {4..9} ; do
	header=`head -n 1 ${1}.window${windowsize}.complimentary `
	headerarray=(${header//;/ })
	echo -n ${headerarray[(${f}-1)]} " "
	tail -n +2 ${1}.window${windowsize}.complimentary |cut -f${f} |sed "s/[1-9][0-9]*/1/g" |paste -sd+ | bc  
done


#the problem here is that the above file *.complimentary has missing lines if the candidate has no hits.  
# so create all_sites.txt with all your candidates "chr pos pos", then run the below
#while read line ; do fgrep "$line" SJ.out.tab.5.5.subset.bed.window500.complimentary ; if  [[ $? -eq  1 ]] ; then echo $line "0 0 0 0 0 0"; fi ; done <all_sites.txt |sed "s/chr//g" |tr " " "\t" >controls_comp_repeats.out
