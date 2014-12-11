##inputs: 1: file.summary
gtfbed=/sc/orga/work/akersn01/ref/Homo_sapiens.GRCh37.74.chr.gtf.bed
#gtfbed=$2

# this was written for circRNA files formatted like this: X2.135252043.70243177.89
# modify as needed, so long as the end result is from the first two lines here is a bed format file. 
#tail -n +2 $1 | cut -f1 |sed 's/X\./chrX./g' |sed 's/Y\./chrY./' |sed 's/^X/chr/g'|sed 's/\./\t/g' |awk '{ print $1,$2,($2+1)}' |tr " " "\t" >${1}.p1.bed
#tail -n +2 $1 | cut -f1 |sed 's/X\./chrX./g' |sed 's/Y\./chrY./' |sed 's/^X/chr/g'|sed 's/\./\t/g' |awk '{ print $1,$3,($3+1)}' |tr " " "\t" >${1}.p2.bed

#alternate start for format 1-12345-6789-12
tail -n +2 $1 | cut -f1 |sed 's/^/chr/g'|sed 's/-/\t/g' |awk '{ print $1,$2,($2+1)}' |tr " " "\t" >${1}.p1.bed
tail -n +2 $1 | cut -f1 |sed 's/^/chr/g'|sed 's/-/\t/g' |awk '{ print $1,$3,($3+1)}' |tr " " "\t" >${1}.p2.bed

bedtools closest -t first -D b -a ${1}.p1.bed -b $gtfbed |cut -f7,12 > ${1}.p1.anno 
bedtools closest -t first -D b -a ${1}.p2.bed -b $gtfbed |cut -f7,12 > ${1}.p2.anno 
sed -i -e "1i\GeneInfo\tDistance" ${1}.p1.anno 
sed -i -e "1i\GeneInfo\tDistance" ${1}.p2.anno 
paste $1 ${1}.p1.anno > ${1}.temp
paste ${1}.temp ${1}.p2.anno > ${1}.annotated
rm ${1}.p1.bed ${1}.p1.anno ${1}.temp ${1}.p2.anno ${1}.p2.bed


##gene name + exon: which columns will depend on your initial input.
#just gene/exone
# $ cut -f8,10 tissue_top_cRNA.annotated |sed 's/exon_number:/\t/g' | sed 's/gene_name:/\t/g' | sed 's/gene_biotype:/\t/g' |cut -f2,3,6,7 |sed 's/"//g' | sed 's/;//g' >tissue_anno_slim
#gene/exon + distance
# $ cut -f8,9,10,11 tissue_top_cRNA.annotated |sed 's/exon_number:/\t/g' | sed 's/gene_name:/\t/g' | sed 's/gene_biotype:/\t/g' | cut -f2,3,5,7,8,10  |sed 's/"//g' | sed 's/;//g' >tissu_anno_slim
