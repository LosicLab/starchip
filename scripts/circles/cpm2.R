args <- commandArgs(TRUE);
#arg1=counts file, 2=cpm cutoff, 3=samples cutoff, 4=gene annotation.
library(limma)
library(edgeR)
library(gplots)
library(methods)
#library(rgl)

readscutoff <- 1000 # minimum number of circRNA reads.  Less and the sample will be eliminated.  
ccounts <- read.table( args[1], row.names = 1, header=T)
ann_c <- read.table( args[4], row.names=1 , header=T)
samples<-ncol(ccounts)
#voom with diagonal design matrix (all replicates of one pretend condition, i.e. variance model including all samples)
dge_cRNA <- DGEList( counts = ccounts, genes = ann_c)
	#keep only samples with > 1000 backsplice (circular) reads
badsamples<-colnames(dge_cRNA$counts)[which(colSums(dge_cRNA$counts) < readscutoff)]
if (length(badsamples) > 0) {
	print(paste("The following samples were eliminated from PCA/heatmap analysis because they have less than", readscutoff, "circRNA reads", sep=" "))
	print(badsamples)
}
if (samples - length(badsamples) > 1) {
#dim(ccounts)
pdf( paste(gsub("countmatrix", "", args[1]), args[2], "cpm_", args[3], "samples", "_variance_PCA.pdf", sep="") )
dge_cRNA <- dge_cRNA[ , colSums(dge_cRNA$counts) > readscutoff, keep.lib.sizes = F]
dge_cRNA <- calcNormFactors(dge_cRNA)
v_cRNA <- voom(dge_cRNA)

write.table( v_cRNA$E, paste( "norm_log2_counts_", gsub("countmatrix", "", args[1]),  args[2], "cpm_", args[3], "samples.txt", sep = ""), row.names = TRUE, col.names = TRUE, sep = "\t", quote = FALSE)

#scaled PCA
pca.voom <- as.data.frame( prcomp(t(na.omit(v_cRNA$E)) , scale = TRUE, center = TRUE)$x)
pca.voom.2 <- prcomp(t(na.omit(v_cRNA$E)), scale = TRUE, center = TRUE)
#
if (ncol(pca.voom) > 1)	{
plot(PC2 ~ PC1, data = pca.voom, pch=21, main = "(SVD, scaled, centered) PCA 1,2 on variance-stabilized counts")
text( pca.voom.2$x[,1], pca.voom.2$x[,2], strtrim( colnames(ccounts), 32), pos = 3, cex=0.5)
}
#
if (ncol(pca.voom) > 2) {
	plot(PC3 ~ PC2, data = pca.voom, pch=21, main = "(SVD, scaled, centered) PCA 2,3 on variance-stabilized counts")
	text( pca.voom.2$x[,2], pca.voom.2$x[,3], strtrim( colnames(ccounts), 32), pos = 3, cex=0.5)
}
#
dev.off()

#clunky: would be better to display the actual table with nice spacing on the heatmap.2, but couldn't find better way
cfeature <- as.matrix(paste(rownames(ann_c), ann_c$Gene.s. ,sep = "_"))

#dim(cfeature)

#upper quantile of circles
#as.numeric(as.matrix(quantile(rowMeans(v_cRNA$E), probs = .75 )))

idx<-order(apply(v_cRNA$E, 1, mad), decreasing=T)[1:100]
limit<-quantile(v_cRNA$E[idx,], probs=.98)
#make heatmap: don't know if the margins are going to be good with very long sample names in the circ count matrix

pdf(paste(gsub("countmatrix", "", args[1]), 'heatmap.pdf', sep=""), height = 10, width = 10)
my_palette <- colorRampPalette(c("red", "black", "green"))(n = 20)
heatmap.2(pmin(v_cRNA$E[idx,], limit), labRow = cfeature[idx], Rowv = T, scale = 'none', Colv = T, trace = 'none', srtCol = 90, cexRow = 0.5, cexCol=0.7, margins = c(7, 7), col = my_palette)

dev.off()

}

