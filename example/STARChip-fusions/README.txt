# STARChip-fusions Example Run:

# Modify hg38.params.txt to point to your generated reference files
# (use your favorite text editor to change the values for refbed, refFasta, and repeatbed ) 

#Run STARChip Fusions
/path/to/starchimp/starchip-fusions.pl Fusions STARout/sample1/Chimeric.out.junction hg38.params.txt

#The output should match what is found in ExpectedOutput/
