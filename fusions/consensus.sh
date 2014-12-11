#! /bin/bash
module load optitype samtools

chrom1=$1
pos1=$2
chrom2=$3
pos2=$4
junction=$5
sam=$6
fusionname=$7
refseq=$8

#create fusion_ref.fa, index it
echo ">"${fusionname} >${fusionname}_ref.fa
echo ${refseq} >>${fusionname}_ref.fa
samtools faidx ${fusionname}_ref.fa
#pull IDs
grep ${chrom1} ${junction} |grep ${pos1} |grep ${chrom2} | grep ${pos2} | cut -f10 >${fusionname}.ID 
#pull sequences from sam file
fgrep -f ${fusionname}.ID ${sam} |cut -f1,10 |sort -u |sed 's/^/>/' |tr "\t" "\n" >${fusionname}.fasta
#align sequences to ref
razers3 -i 90 -ng -o ${fusionname}.sam ${fusionname}_ref.fa ${fusionname}.fasta
#sort aligned reads
samtools view -bS ${fusionname}.sam | samtools sort - ${fusionname}.sorted
#generate consensus
samtools mpileup -uf ${fusionname}_ref.fa ${fusionname}.sorted.bam |bcftools view -cg - | vcfutils.pl vcf2fq #> ${fusionname}.consensus
rm ${fusionname}.ID ${fusionname}.fasta ${fusionname}.sam ${fusionname}_ref.fa ${fusionname}_ref.fa.fai ${fusionname}.sorted.bam
