if [ "$#" -ne 2 ]; then
    echo "Illegal number of parameters.  usage: gtf2bed.sh /path/to/gtf.file /path/to/desired/output/directory/"
    exit 1
fi
filename=`basename $1`
#create a bed format version of the supplied GTF.
sed 's/ /:/g' $1 |sed 's/;:/;/g' | awk 'BEGIN { FS = "\t" }{ print $1,$4,$5,$9,$6,$7,$2,$3 }' OFS="\t" |sort -k1,1 -k2,2n -k3,3n > ${2}${filename}.bed
#create a .genome file for use with bedtools
cut -f1,3 ${2}${filename}.bed |awk '{ if ($1 != chrm) print chrm,pos ; chrm=$1 ; pos=$2 } END {print chrm,pos}' OFS="\t" |grep -v -P "^#" > ${2}${filename}.genome
