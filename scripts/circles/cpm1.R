args <- commandArgs(TRUE);
#arg1=star dir, 2=output file name 3=strand imbalance details output

#gather up the read count files
flist <- list.files( path = args[1], pattern = "*ReadsPerGene.out.tab$", full.names = T, recursive = T)
#remove empty files, convert to count matrix
#gcounts_list <- lapply( flist, read.table, skip = 4 )
gcounts_list <- lapply(flist, function(x) { 
	tmp <- try(read.table(x, skip = 4))
	if (!inherits(tmp, 'try-error')) tmp
})
empties <-unlist(lapply(gcounts_list, is.null))
if(sum(empties) > 0 ) { 
	print("The following count files were found to be empty, they will be skipped, you may need to rerun STAR for these files")
	print(flist[empties])
}
gcounts_list <- gcounts_list[! empties]	
gcounts <- as.data.frame( sapply( gcounts_list, function(x) x[,2] ) )
# calculate the ratio of strand read support
gcounts_ratio <- (as.data.frame( sapply( gcounts_list, function(x) x[,3] ) )+0.5)/( as.data.frame( sapply( gcounts_list, function(x) x[,4] ) )+0.5)
# turn ratios < 1 to > 1 
gcounts_ratio_adjusted <- gcounts_ratio
gcounts_ratio_adjusted[gcounts_ratio_adjusted < 1] <- 1/gcounts_ratio_adjusted[gcounts_ratio_adjusted < 1]
#add colnames + rownames
mycolnames <- gsub("STARout\\/", "", gsub("\\/ReadsPerGene.out.tab", "", flist[! empties]))
colnames(gcounts_ratio_adjusted) <- mycolnames
colnames(gcounts) <- mycolnames
rownames(gcounts) <- gcounts_list[[1]]$V1
rownames(gcounts_ratio_adjusted) <- gcounts_list[[1]]$V1

#create output tables
output_all<-data.frame( medianStrandRatio=apply(gcounts_ratio_adjusted, 1, median), samplesWithStrandImbalanceOver2=rowSums(gcounts_ratio_adjusted > 2)  )
output_all$Kept<-output_all$samplesWithStrandImbalanceOver2 <= (ncol(gcounts_ratio_adjusted)/2) 
countMatrix<-gcounts[rownames(output_all)[which(output_all$Kept)],,drop=F]
output_all2<-cbind(output_all, gcounts_ratio_adjusted)

#write output
write.table(countMatrix, file=args[2], sep="\t", quote=F)
write.table(gcounts, file=paste(args[2], "NotStrandImbalanceReduced", sep="."), sep="\t", quote=F)
write.table(output_all2, file=args[3], sep="\t", quote=F)
