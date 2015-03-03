#! /bin/bash
module load samtools

chrom1=$1
pos1=$2
chrom2=$3
pos2=$4
junction=$5
sam=$6
fusionname=$7
refseq=$8
script_dir=$9

#create fusion_ref.fa, index it
#echo ">"${fusionname} >${fusionname}_ref.fa
#echo ${refseq} >>${fusionname}_ref.fa
#samtools faidx ${fusionname}_ref.fa
#pull IDs
grep ${chrom1} ${junction} |grep ${pos1} |grep ${chrom2} | grep ${pos2} | cut -f10 >${fusionname}.ID 
#pull sequences from sam file
fgrep -f ${fusionname}.ID ${sam} |cut -f1,10 |sort -u |sed 's/^/>/' |tr "\t" "\n" >${fusionname}.fasta
avgAS=`fgrep -f ${fusionname}.ID ${sam} |sed 's/.*AS:i://' |sed 's/\s.*//' |awk '{ total += $1; count++ } END { print total/count }' `
	#align sequences to ref
	#razers3 -i 90 -ng -o ${fusionname}.sam ${fusionname}_ref.fa ${fusionname}.fasta
mafft --reorder --adjustdirection --auto --quiet ${fusionname}.fasta >${fusionname}.msa
${script_dir}consensus.pl ${fusionname}.msa
	#sort aligned reads
	#samtools view -bS ${fusionname}.sam | samtools sort - ${fusionname}.sorted
	#generate consensus
	#samtools mpileup -uf ${fusionname}_ref.fa ${fusionname}.sorted.bam |bcftools view -cg - | vcfutils.pl vcf2fq #> ${fusionname}.consensus
echo $avgAS
rm ${fusionname}.ID ${fusionname}.fasta ${fusionname}.msa
