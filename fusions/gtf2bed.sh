sed 's/ /:/g' $1 |sed 's/;:/;/g' | awk 'BEGIN { FS = "\t" }{ print $1,$4,$5,$9,$6,$7,$2,$3 }' OFS="\t" |sort -k1,1 -k2,2n -k3,3n > ${1}.bed
cut -f1,3 ${1}.bed |awk '{ if ($1 != chrm) print chrm,pos ; chrm=$1 ; pos=$2 } END {print chrm,pos}' OFS="\t" |grep -v -P "^#" > ${1}.genome
