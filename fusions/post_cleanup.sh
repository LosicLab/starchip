#usage: post_cleanup.sh ListOfFiles output_prefix  Go

if [ -z ${1+x} ] || [ $1 == "-h" ] || [ $1 == "--help" ] || [ $1 == "-help" ]; then
	echo "usage: $ post_cleanup.sh ListOfFiles output_prefix  Go"
	echo "use the 3rd argument if you want gene level summaries performed"
	exit 1
fi

mylist=$1
myprefix=$2
while read eachfile ; do
	while read eachline ; do
		echo -e $eachfile"\t"$eachline ;
	done < ${eachfile} ;
done < $mylist |grep -v Partner |tr " " "\t" > $myprefix

if [ -z ${3+x} ]; then
	echo "Not performing Gene summaries, add any 3rd argument to perform gene summaries";
else
	cut -f1,8,10 $myprefix | awk '{ if ($2 >= $3) print $1,$2,$3 ; else print $1,$3,$2 } ' OFS="\t" |sort |uniq |cut -f2,3 |sort | uniq -c |sort -gr > ${myprefix}.genepairs
	cut -f1,8,10 $myprefix | awk '{ if ($2 >= $3) print $1,$2,$3 ; else print $1,$3,$2 } ' OFS="\t" |sort |uniq |awk '{print $1,$2"\n"$1,$3}' OFS="\t" |cut -f2 | sort |uniq -c |sort -gr > ${myprefix}.genecount_uniqpairs
	cut -f1,8,10 $myprefix | awk '{ if ($2 >= $3) print $1,$2,$3 ; else print $1,$3,$2 } ' OFS="\t" |sort |uniq |awk '{print $1,$2"\n"$1,$3}' OFS="\t" |sort | uniq | cut -f2 | sort |uniq -c |sort -gr > ${myprefix}.genecount_uniqinds
fi
