#!/bin/bash

#usage: genes_to_coords.sh   genes_file   gtf.bed

while read gene ; do 
	fgrep -w $gene $2 |cut -f1-3 > temp.${gene}.txt ; 
	if [[ $(wc -l <temp.${gene}.txt) -ge 1 ]] 
	then
		chromosome=`head -1 temp.${gene}.txt |cut -f1 `;
		pos1=`sort -k2,2n temp.${gene}.txt |head -1 |cut -f2 `; 
		pos2=`sort -k3,3nr temp.${gene}.txt |head -1 |cut -f3 `; 
		echo -e "0\t.\t.\t"${chromosome}"\t"${pos1}"\t.\t"${chromosome}"\t"${pos2}"\t.\t.\t.\t.\t0\t.\t.\t"${gene}
	else
		>&2 echo "Cannot find gene ${gene}"
	fi
	rm temp.${gene}.txt
done < $1
