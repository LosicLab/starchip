#usage coordinates2genes.sh file.summary gtf.bed repeats.bed
##inputs: 1: file.summary
#gtfbed=/sc/orga/work/akersn01/ref/Homo_sapiens.GRCh37.74.chr.gtf.bed
gtfbed=$2
repeatbed=$3

# this was written for fusion files that start like this: chr2:135252043:-	chr10:70243177:+
# modify as needed, so long as the end result is from the first two lines here is a bed format file.
tail -n +2 $1 | cut -f1 |sed 's/:/\t/g' |awk '{ print $1,($2-1),($2+1),$1":"$2":"$3 }' |tr " " "\t" >${1}.p1.bed
tail -n +2 $1 | cut -f2 |sed 's/:/\t/g' |awk '{ print $1,($2-1),($2+1),$1":"$2":"$3 }' |tr " " "\t" >${1}.p2.bed
#gene annotation
bedtools closest -t first -D b -a ${1}.p1.bed -b $gtfbed |cut -f4,8,13> ${1}.p1.anno
bedtools closest -t first -D b -a ${1}.p2.bed -b $gtfbed |cut -f4,8,13> ${1}.p2.anno
#repeats
bedtools intersect -a ${1}.p1.bed -b $repeatbed -c |cut -f5 > ${1}.p1.rep
bedtools intersect -a ${1}.p2.bed -b $repeatbed -c |cut -f5 > ${1}.p2.rep
paste ${1}.p1.rep ${1}.p2.rep |awk '{ if ($1>0 && $2>0) print "2" ; else if ($1 >0 || $2 >0) print "1" ; else print "0"}' >${1}.repeats
#clean up
sed -i -e "1i\Repeats" ${1}.repeats
sed -i -e "1i\ID\tGeneInfo\tDistance" ${1}.p1.anno
sed -i -e "1i\ID\tGeneInfo\tDistance" ${1}.p2.anno
paste $1 ${1}.repeats >${1}.temp0
paste ${1}.temp0 ${1}.p1.anno > ${1}.temp
paste ${1}.temp ${1}.p2.anno > ${1}.annotated
rm ${1}.p1.bed ${1}.p2.bed ${1}.p1.anno ${1}.p2.anno ${1}.temp ${1}.p1.rep ${1}.p2.rep ${1}.repeats ${1}.temp0
