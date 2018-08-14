#!/bin/bash

if [ "$#" -ne 3 ]; then
    echo "Spliced circRNA Sequence usage: circRNA_sequence.sh /path/to/myCircs.consensus /path/to/genome.fasta.file /path/to/desired/output_dir/"
    exit 1
fi
myname=$(basename $1)
cut -f1,17-19 $1 | tail -n +2 | sed 's/[:-]/\t/g' | awk '{ print $1,$2,$3,$1":"$2"-"$3,".",".",".",".",".",$4,$6,$5 }' OFS="\t" > ${3}/${myname}.bed

bedtools getfasta -fi $2 -bed ${3}/${myname}.bed -fo ${3}/${myname}.fa -name -split
