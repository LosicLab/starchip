#!/bin/bash

#usage coordinates2genes.sh file.summary gtf.bed repeats.bed
##inputs: 1: file.summary
#gtfbed=/sc/orga/work/akersn01/ref/Homo_sapiens.GRCh37.74.chr.gtf.bed
gtfbed=$2
repeatbed=$3
genomefile=`echo $gtfbed |sed 's/.bed$/.genome/'`
# this was written for fusion files that start like this: chr2:135252043:-      chr10:70243177:+
# turn coordinates into a sorted bed file.
tail -n +2 $1 | cut -f1,2 |tr "\t" "\n" |sed 's/:/\t/g' |awk '{ print $1,($2-1),($2+1),$1":"$2":"$3 }' |tr " " "\t" | sort -k1,1 -k2,2n -k3,3n >${1}.p1.bed
#gene annotation
bedtools closest -g $genomefile -t first -D b -a ${1}.p1.bed -b $gtfbed |cut -f4,8,13> ${1}.p1.anno
#repeats
bedtools intersect -sorted -a ${1}.p1.bed -b $repeatbed -c |cut -f4,5 > ${1}.p1.rep

#header management:
tail -n +2 $1 > ${1}.headerless
head -n 1 $1 > ${1}.annotated

#now we have gene annotations in the .anno file, and repeat counts in the .rep file. 
# use grep to join fusion partners back to the same line.
while read line ; do
	linearray=(${line/\t/ / })
	coord1=${linearray[0]}
	coord2=${linearray[1]}
	#all grep calls have -m 1 which tells grep to stop after one match.
	#find the number of repeats for each partner (0 or 1), sum them. 
	#echo "running fgrep -m 1 $coord1 ${1}.p1.rep "
	coord1repgrep=$(fgrep -m 1 $coord1 ${1}.p1.rep |awk '{print $2}' )
	#echo "coord1repgrep is "$coord1repgrep
	coord2repgrep=$(fgrep -m 1 $coord2 ${1}.p1.rep |awk '{print $2}' )
	#let "repeatssum = $coord1repgrep + $coord2repgrep"
	#((repeatssum = $coord1repgrep + $coord2repgrep))
	#thanks to Tianshi Lu for providing the fix here:
	repeatssum=$(($coord1repgrep + $coord2repgrep))
	#assign gene annotations to vars. 
	coord1genegrep=$(fgrep -m 1 $coord1 ${1}.p1.anno )
        if  [[ $? -eq  0 ]] ; then 
       	        gene1anno=$coord1genegrep
        else
        	gene1anno="${coord1}\tgene_name:\"NoGene\"\tNA"
        fi
	coord2genegrep=$(fgrep -m 1 $coord2 ${1}.p1.anno )
        if  [[ $? -eq  0 ]] ; then 
                gene2anno=$coord2genegrep
        else
                gene2anno="${coord2}\tgene_name:\"NoGene\"\tNA"
        fi
	gene1anno2=`echo $gene1anno | sed 's/\s\.\s-1/\tgene_name:NoGene\tNA/' `
	gene2anno2=`echo $gene2anno | sed 's/\s\.\s-1/\tgene_name:NoGene\tNA/' `
	echo -e $line"\t"$repeatssum"\t"$gene1anno2"\t"$gene2anno2
done < ${1}.headerless  | tr " " "\t" >> ${1}.annotated

#clean up

sed -i -e '1 s/$/\tRepeats\tPartner1\tP1GeneInfo\tP1GeneDistance\tPartner2\tP2GeneInfo\tP2GeneDistance/' ${1}.annotated
rm ${1}.p1.bed ${1}.p1.anno ${1}.p1.rep ${1}.headerless


#output looks like:
# chr2:135252043:-      chr10:70243177:+  ANYTHINGELSE	numbRepeats	p1	p1GeneInfo	p1genedist	p2	p2GeneInfo	p2genedist
