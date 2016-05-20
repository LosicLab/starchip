#!/bin/bash

if [ "$#" -ne 3 ]; then
    echo "Illegal number of parameters.  usage: gtf2bed.sh /path/to/gtf.file /path/to/genome.fasta.file /path/to/desired/output/directory/"
    exit 1
fi
filename=`basename $1`
#create a bed format version of the supplied GTF.
echo "creating ${3}${filename}.bed"
#sed 's/ /:/g' $1 |sed 's/;:/;/g' | awk 'BEGIN { FS = "\t" }{ print $1,$4,$5,$9,$6,$7,$2,$3 }' OFS="\t" |sort -k1,1 -k2,2n -k3,3n > ${3}${filename}.bed
#create a 2nd bed format GTF, with only features that have exon_numbers, for better cRNA annotation
echo "creating ${3}${filename}.exons.bed"
#sed 's/ /:/g' $1 |sed 's/;:/;/g' | awk 'BEGIN { FS = "\t" }{ print $1,$4,$5,$9,$6,$7,$2,$3 }' OFS="\t" |fgrep exon_number | sort -k1,1 -k2,2n -k3,3n > ${3}${filename}.exons.bed
#create a .genome file for use with bedtools
 #get dir of this script
echo "creating ${3}${filename}.genome"
DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
#${DIR}/fasta_getlength.pl $2 | sort -k1,1 > ${3}${filename}.genome 
ln -s $(cd $(dirname "${3}${filename}.genome") && pwd -P)/$(basename "${3}${filename}.genome") ${3}${filename}.exons.genome
#cut -f1,3 ${3}${filename}.bed |awk '{ if ($1 != chrm) print chrm,pos ; chrm=$1 ; pos=$2 } END {print chrm,pos}' OFS="\t" |grep -v -P "^#" > ${2}${filename}.genome
