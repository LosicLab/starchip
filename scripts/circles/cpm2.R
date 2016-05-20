args <- commandArgs(TRUE);
#arg1=counts file, 2=cpm cutoff, 3=samples cutoff, 4=gene annotation.
library(limma)
library(edgeR)
#library(rgl)

ccounts <- read.table( args[1], row.names = 1, header=T)
ann_c <- read.table( args[4], row.names=1 , header=T)
#dim(counts)

pdf( paste(args[1], "_cpm_", args[2], "_samples_", args[3], "_variance_PCA.pdf", sep="") )
#voom with diagonal design matrix (all replicates of one pretend condition, i.e. variance model including all samples)
dge_cRNA <- DGEList( counts = ccounts, genes = ann_c)
dge_cRNA <- dge_cRNA[ , colSums(dge_cRNA$counts) > 1000, keep.lib.sizes = F]
dge_cRNA <- calcNormFactors(dge_cRNA)
v_cRNA <- voom(dge_cRNA)

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


