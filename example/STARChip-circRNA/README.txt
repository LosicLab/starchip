#Example data + commands for runnign STARChip-circRNA

#Modify starchip-circles.params to include the path to your gtf.bed and fasta files
# (Use your favorite text editor)

#Run the following command from this directory
/path/to/starchip/starchip-circles.pl stardirs.txt starchip-circles.params

#This should create Step1.sh and Step2.sh run these scripts in successsion:
./Step1.sh
./Step2.sh

#The output should match what is found in ExpectedOutput/
