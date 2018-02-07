#!/bin/bash

#inputs: 1: cRNA.consensus (or anything with the cRNA in 1st column.
gtfbed=$2
genomefile=`echo $gtfbed |sed 's/.bed$/.genome/'`

#cRNA format 1-12345-6789-12 or 1:2345-6789
tail -n +2 $1 | cut -f1 |sed 's/[-:]/\t/g' |awk '{ print $1,$2,$2,$1":"$2"-"$3":p1\n"$1,$3,$3,$1":"$2"-"$3":p2"}' OFS="\t"| sort -k1,1 -k2,2n -k3,3n  >${1}.p1.bed

bedtools closest -g $genomefile -nonamecheck -t first -D b -a ${1}.p1.bed -b $gtfbed |cut -f4,8,10,13 > ${1}.p1.anno
sed -i -e "1i\GeneInfo\tStrand\tDistance" ${1}.p1.anno
#paste $1 ${1}.p1.anno > ${1}.temp
#paste ${1}.temp ${1}.p2.anno > ${1}.annotated
#rm ${1}.p1.bed ${1}.p1.anno ${1}.temp ${1}.p2.anno ${1}.p2.bed

#header management:
tail -n +2 $1 > ${1}.headerless
head -n 1 $1 > ${1}.annotated

#now we have gene annotations in the .anno file
# use grep to join fusion partners back to the same line.
while read line ; do
        linearray=(${line/\t/ / })
        cRNA=${linearray[0]}
	search1=${cRNA}":p1"
	search2=${cRNA}":p2"
        #all grep calls have -m 1 which tells grep to stop after one match.
        #assign gene annotations to vars.
        coord1genegrep=$(fgrep -m 1 $search1 ${1}.p1.anno )
        if  [[ $? -eq  0 ]] ; then 
                gene1anno=$coord1genegrep
        else
            	gene1anno=".\t.\t-1"
        fi
	coord2genegrep=$(fgrep -m 1 $search2 ${1}.p1.anno )
        if  [[ $? -eq  0 ]] ; then 
                gene2anno=$coord2genegrep
        else
                gene2anno=".\t.\t-1"
        fi
        gene1anno2=`echo $gene1anno | sed 's/\.\s\.\s-1/gene_name:NoGene\tNA\tNA/' | awk '{print $2,$3,$4}' OFS="\t"`
        gene2anno2=`echo $gene2anno | sed 's/\.\s\.\s-1/gene_name:NoGene\tNA\tNA/' | awk '{print $2,$3,$4}' OFS="\t"`
        echo -e $line"\t"$gene1anno2"\t"$gene2anno2
done < ${1}.headerless  | tr " " "\t" >> ${1}.annotated

#clean up
sed -i -e '1 s/$/\tP1GeneInfo\tP1Strand\tP1GeneDistance\tP2GeneInfo\tP2Strand\tP2GeneDistance/' ${1}.annotated
rm ${1}.p1.bed ${1}.p1.anno ${1}.headerless


##gene name + exon: which columns will depend on your initial input.
#just gene/exone
# $ cut -f8,10 tissue_top_cRNA.annotated |sed 's/exon_number:/\t/g' | sed 's/gene_name:/\t/g' | sed 's/gene_biotype:/\t/g' |cut -f2,3,6,7 |sed 's/"//g' | sed 's/;//g' >tissue_anno_slim
#gene/exon + distance
# $ cut -f8,9,10,11 tissue_top_cRNA.annotated |sed 's/exon_number:/\t/g' | sed 's/gene_name:/\t/g' | sed 's/gene_biotype:/\t/g' | cut -f2,3,5,7,8,10  |sed 's/"//g' | sed 's/;//g' >tissu_anno_slim

