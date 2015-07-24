#!/bin/bash
##requires file ./junctionfiles.txt  (see last line)
##this file should include full paths to Chimeric.out.junction file, 1 per line.  ie /user/data/star_runs/sample_2/Chimeric.out.junction
##note that the script will also look for Chimeric.out.sam in the same folder.

mkdir -p log
mkdir -p runfiles
mkdir -p data
params=$1
while read line ; do
ID=` echo $line | sed "s/\/Chimeric.out.junction//g" | sed 's/.*\///'`
echo "
#! /bin/bash
#BSUB -J fusionfilter
#BSUB -e log/"${ID}".e
#BSUB -o log/"${ID}".o
#BSUB -q low
#BSUB -P acc_PBG
#BSUB -W 1:00
#BSUB -n 2
###BSUB -P acc_PBG
###BSUB -m manda
#BSUB -R "rusage[mem=5000]"
#BSUB -R "span[hosts=1]"
module load bedtools samtools
/hpc/users/akersn01/scripts/starchimp/fusions/fusions-from-star.pl data/"${ID}" "${line}" "${params}"
" > runfiles/${ID}.lsf
done < $2
