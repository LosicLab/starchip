args <- commandArgs(TRUE);
library(limma)
library(edgeR)
#library(rgl)

counts <- read.table( args[1], row.names = 1, header=T)
#dim(counts)

pdf( paste(args[1], "_cpm_", args[2], "_samples_", args[3], "_variance_PCA.pdf", sep="") )
#voom with diagonal design matrix (all replicates of one pretend condition, i.e. variance model including all samples)
v <- voom(counts[rowSums(cpm(counts) >= as.numeric(args[2]) ) >= as.numeric(args[3]), ], plot = TRUE )
#dim(v$E)

#scaled PCA
pca.voom <- as.data.frame( prcomp(t(na.omit(v$E)) , scale = TRUE, center = TRUE)$x)
pca.voom.2 <- prcomp(t(na.omit(v$E)), scale = TRUE, center = TRUE)
#
plot(PC2 ~ PC1, data = pca.voom, pch=21, main = "(SVD, scaled, centered) PCA 1,2 on variance-stabilized counts")
text( pca.voom.2$x[,1], pca.voom.2$x[,2], strtrim( colnames(counts), 32), pos = 3, cex=0.5)
#
plot(PC3 ~ PC2, data = pca.voom, pch=21, main = "(SVD, scaled, centered) PCA 2,3 on variance-stabilized counts")
text( pca.voom.2$x[,2], pca.voom.2$x[,3], strtrim( colnames(counts), 32), pos = 3, cex=0.5)
#
dev.off()


write.table( v$E, paste( "norm_log2_counts", args[1],  "cpm", args[2], "samples", args[3], sep = "_"), row.names = TRUE, col.names = TRUE, sep = "\t", quote = FALSE)

#note that cor( 2^(v$E), normalizedCounts) should have diagonal entries very close to +1, where v$E here is calculated by setting <sample_threshold> above to zero
#(i.e. no count filtering whatsoever), and normalizedCounts is
#
#count_table <- raw counts Hardik gives you
#expt_design <- data.frame(row.names=colnames(count_table),condition=c(colnames(count_table))
#conditions=expt_design$condition
#library(DESeq)
#data <- newCountDataSet(count_table,conditions)
#data <- estimateSizeFactors(data)
#normalizedCounts <-  t( t(counts(data)) / sizeFactors(data) )
#
#
#

