gtf=$2 
fasta=$3
cut -f7,9 $1 |tail -n +2 |sort |uniq > ${1}.genes
while read line ; do 
	gene1=`echo $line |sed 's/ .*//'`
	gene2=`echo $line |sed 's/.* //'`
	fgrep "gene_name \"${gene1}\"" $gtf |sed 's/.* transcript_id//' |sed 's/; gene_type.*//' |sort |uniq |sed 's/"//g' >${1}.${gene1}.ensembl
	fgrep "gene_name \"${gene2}\"" $gtf |sed 's/.* transcript_id//' |sed 's/; gene_type.*//' |sort |uniq |sed 's/"//g' >${1}.${gene2}.ensembl
	while read ID ;do
		samtools faidx $fasta ${ID} 
	done < ${1}.${gene1}.ensembl >${1}.${gene1}.fasta
	while read ID ;do
		samtools faidx $fasta ${ID} 
	done < ${1}.${gene2}.ensembl >${1}.${gene2}.fasta
	tophit=`blastn -query ${1}.${gene1}.fasta -subject ${1}.${gene2}.fasta -outfmt "6 qseqid sseqid evalue pident" |sort -k4,4g |head -n 1`
	echo $gene1 $gene2 $tophit
	rm ${1}.${gene1}.ensembl ${1}.${gene1}.fasta ${1}.${gene2}.ensembl ${1}.${gene2}.fasta
done < ${1}.genes > ${1}.blastout
