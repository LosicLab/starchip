#!/bin/bash
##requires file ./junctionfiles.txt  (see last line)
##this file should include full paths to Chimeric.out.junction file, 1 per line.  ie /user/data/star_runs/sample_2/Chimeric.out.junction
##note that the script will also look for Chimeric.out.sam in the same folder.

#get dir of this script
DIR=$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )
#get directory 2 level up for starchimp-fusions.pl
#DIR=`echo $DIR |sed 's/\/[a-zA-Z0-9\._-]*$//' `;
#DIR=`echo $DIR |sed 's/\/[a-zA-Z0-9\._-]*$//' `;
DIR=`echo $DIR |sed 's/\/[^/]*$//' `;
DIR=`echo $DIR |sed 's/\/[^/]*$//' `;

mkdir -p log
mkdir -p runfiles
mkdir -p data
params=$1
junctionFiles=$2
while read line ; do
ID=` echo $line | sed "s/\/Chimeric.out.junction//g" | sed 's/.*\///'`
echo "
#! /bin/bash
#BSUB -J fusionfilter
#BSUB -e log/"${ID}".e
#BSUB -o log/"${ID}".o
#BSUB -q alloc
#BSUB -P acc_PBG
#BSUB -W 1:00
#BSUB -n 2
###BSUB -P acc_PBG
###BSUB -m manda
#BSUB -R "rusage[mem=5000]"
#BSUB -R "span[hosts=1]"
module load starchimp
${DIR}/starchimp-fusions.pl data/"${ID}" "${line}" "${params}"
" > runfiles/${ID}.lsf
done < $junctionFiles
