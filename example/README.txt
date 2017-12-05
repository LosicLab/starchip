# Example data + commands for running STARChip  

# You will need to run the following for fusion or circRNA examples:  


# Download matching reference fasta and gtf files for Hg38 
# These files will total ~ 7GB, use a directory with enough space
wget ftp://ftp.sanger.ac.uk/pub/gencode/Gencode_human/release_26/gencode.v26.primary_assembly.annotation.gtf.gz
wget ftp://ftp.sanger.ac.uk/pub/gencode/Gencode_human/release_26/GRCh38.primary_assembly.genome.fa.gz
# Fusions only: download repeat regions:
wget http://hgdownload.soe.ucsc.edu/goldenPath/hg38/database/rmsk.txt.gz

# Decompress the references 
gunzip gencode.v26.primary_assembly.annotation.gtf.gz
gunzip GRCh38.primary_assembly.genome.fa.gz
# Fusions only:
gunzip rmsk.txt.gz


# Set up the references 
/path/to/starchip/setup.sh /path/to/gencode.v26.primary_assembly.annotation.gtf /path/to/GRCh38.primary_assembly.genome.fa /path/to/desired/output/dir/
# Fusions only:
cut -f6-8 rmsk.txt > hg38.repeats.bed

#Now you're set up to run STARChip.  Move to examples/STARChip-circRNA or examples/STARChip-fusions for usage examples. 
