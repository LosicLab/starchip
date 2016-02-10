#!/usr/bin/perl
use strict;
use warnings;

print "ID\tGene(s)\tExon1\tExon2\tStrand\tSpliceType\n";
while (my $x = <>) { #read in argv[0] or stdin:  
	#format should be 0:cRNA 1:p1GeneInfo 2:p1Strand 3:p1Distance 4:p2GeneInfo 5:p2Strand 6:p2Distance
	my ($id, $p1Info, $p1Strand, $p1Distance, $p2Info, $p2Strand, $p2Distance) = split(/\s+/, $x); 
	my $strand ; my $GeneID ; my $p1Exon ; my $p2Exon ; my $circleType ; 
	##Get the Gene Names.
	my @p1geneinfo=split(/\;/, $p1Info);
	my @p2geneinfo=split(/\;/, $p2Info);
	my @gene1name=grep(/gene_name/, @p1geneinfo);
        my @gene2name=grep(/gene_name/, @p2geneinfo);
	## If common gene name not available, use gene_id
	if (! @gene2name || ! @gene1name) { 
		@gene1name=grep(/gene_id/, @p1geneinfo);
		@gene2name=grep(/gene_id/, @p2geneinfo);
	}
	if ($gene1name[0] eq $gene2name[0]) {
		$strand=$p1Strand;
		$GeneID=$gene1name[0];
		#determine type of splice
		if ($p1Distance == 0 && $p2Distance == 0 ) {
			$circleType="Exon-Exon";
		}
		elsif ($p1Distance == -1 && $p2Distance == 1) {
			$circleType="Intron-Lariat"; 
		}
		elsif (($p1Distance == 0 && $p2Distance != 0 ) || ($p1Distance != 0 && $p2Distance == 0 )) {
			$circleType="Exon-Intron"; 
		}
		else { $circleType="Unannotated"; }
	}
	#if we have different gene names
	elsif ($gene1name[0] ne $gene2name[0] ) {
		#first check if it is in 1 gene, not another
		if ($p1Distance == 0 && $p2Distance != 0 ) {
			$GeneID=$gene1name[0]; 
			$strand=$p1Strand;
			$circleType="Exon-Intron"; 
		}
		elsif ($p2Distance == 0 && $p1Distance != 0 ) {
                        $GeneID=$gene2name[0]; 
                        $strand=$p2Strand;
                        $circleType="Exon-Intron"; 
                }
		## If it is within two different genes
		elsif ($p1Distance == 0 && $p2Distance == 0 ) {
			$GeneID=$gene1name[0] . "/" . $gene2name[0]; 
			$circleType="Exon-Exon"; 
			if ($p1Strand eq $p2Strand ) {
				$strand=$p1Strand; 
			}
			else { $strand=$p1Strand . "/" . $p2Strand ; 
			}
		}
		# If it isn't directly in any gene 
		elsif ($p1Distance != 0 && $p2Distance != 0 ) {
			if (abs( $p1Distance ) > abs ($p2Distance) ){
				$GeneID=$gene1name[0]; 
				$strand=$p1Strand;
			}
			else {
				$GeneID = $gene2name[0];
				$strand=$p2Strand;
			}
			$circleType="Unannotated";
		}
	}
	$GeneID =~ s/gene_name:"//g;
	$GeneID =~ s/gene_id:"//g;
	$GeneID =~ s/"//g;
	##Get Exon Information
	my @p1exon=grep(/exon_number/, @p1geneinfo);
	my @p2exon=grep(/exon_number/, @p2geneinfo);
	$p1Exon = $p1exon[0];
	$p2Exon = $p2exon[0];
	$p1Exon =~ s/exon_number://g;
	$p2Exon =~ s/exon_number://g;
	$p1Exon =~ s/"//g;
	$p2Exon =~ s/"//g;
	if ($strand eq "-" ) { 
		#swap the exons
		my $temp = $p1Exon; 
		$p1Exon = $p2Exon; 
		$p2Exon = $temp; 
	}
	print "$id\t$GeneID\t$p1Exon\t$p2Exon\t$strand\t$circleType\n"; 
} 
