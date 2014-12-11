ONCOFUSE_JAR=~/software/oncofuse-v1.0.6/Oncofuse.jar
#The first 9 columns give information about the chimeric junction:
#1: chromosome of donor
#2: first base of the intron of the donor
#3: strand of the donor
#4: chromosome of the acceptor
#5: first base of the intron of the acceptor
#6: strand of the acceptor
#7: junction type: -1=junction is between the mates, 1=GT/AG, 2=CT/AC
#8: repeat length to the left of the junction
#9: repeat length to the right of the junction

#Columns 10-14 describe the alignments of the two chimeric segments, it is SAM like. Alignments are given with respect to the + strand
#10: read name
#11: first base of the first segment (on the + strand)
#12: CIGAR of the first segment
#13: first base of the second segment
#14: CIGAR of the second segment

mkdir -p chimera

#for i in `find . -name "Chimeric.out.junction"`; do nom=`echo "${i}" | sed 's/.\/star_//g' | sed 's/_L00.*//g' `; cp "${i}" chimera/"${nom}"_chimeric.out.junction; done
for i in `find . -name "Chimeric.out.junction"`; do nom=`echo "${i}" | sed 's/\.\///' | sed 's/\/.*//' ` ; cp "${i}" chimera/"${nom}"_chimeric.out.junction; done

#collapse readbased chimeric calls from STAR into summary information

for i in `find . -name *_chimeric.out.junction`; do awk '$1!="chrM" && $1!="MT" && $4!="chrM" && $4!="MT" && $7>0 && $8+$9<=5 {print "chr"$1,$2,$3,"chr"$4,$5,$6,$7,$8,$9}' "${i}" | sed 's/\s/\t/g' | sort | uniq -c | sort -k1,1rn | sed -e 's/^[ \t]*//' | sed 's/\s/\t/g' | awk '$1 > 10' > "${i}"_filtered; done

#this keeps canonical junctions with repeat length less than 5 and removes mitchondrial chimera: I keep only those junctions with greater than 10 reads supporting it

for i in `find . -name *filtered`; do awk '{print "'${i}'" "\t" $0}' "${i}"; done | sort -k3,3 -k4,4n > star_chimera_summary

module load oncofuse

for i in `find . -name *filtered`; do awk '{print $2,$3,$5,$6, "HEM"}' "${i}" | sed 's/\s/\t/g' > "${i}"_finput; java -Xmx1G -jar $ONCOFUSE_JAR "${i}"_finput coord - "${i}"_foutput; done

for i in *_filtered  ; do
	while read line ; do 
		coord=`echo $line | awk '{ print $2":"$3">"$5":"$6 }'`
		#echo $coord
		grephit="null"
		grephit=$(grep $coord ${i}_foutput)
		if  [[ $? -eq  0 ]] ; then
                	hit=`echo $grephit |awk '{print $4,$5,$12,$20,$21,$22}'`
                else
                      	hit="0" 
                fi
		echo -n $line" "$hit | tr " " "\t" 
		echo ""
	done < $i > ${i}_with_oncofuse
	sed -i '1ireads\tchr\tpos\tstrand\tchr\tpos\tstrand\tjxntype\tleftrepeat\trightrepeat\tcoord\tgene1\tgene2\tpval\tdriverprob\tpred_gene_expr_change' ${i}_with_oncofuse
done 
