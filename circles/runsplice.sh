#usage: runsplice.sh cutoff cpus script_dir data_parent_dir min_subjects
# for use if you want to splice samples after running circle_star with splice turned off. 

cutoff=$1
cpus=$2
DIR=$3
Parent=$4
minSubjLimit=$5

#splice the candidates using the SJ.out.tab file for the individual.  
##temp_subject_IDS.txt was created by circle_star.sh it should just be the individual IDs.  
cat temp_subject_IDS.txt | xargs --max-procs=${cpus} -I {} ${DIR}/circle_splice_individual.sh joinstrands.{} ${Parent}/{} ${cutoff}
##Calculate Statistics on the candidates
${DIR}/candidate_stats.pl circs${cutoff}.${minSubjLimit}.investigate ${cutoff}.spliced
#generates .consensus and .allvariants files.  


