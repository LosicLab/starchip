#### i.e. grep chr8 data/*.summary |sed 's/:/\t/g'> test.chr8 > args[2]
#### args[1] = gtf 
### args[3] = outputFile ; 
options(stringsAsFactors=FALSE)
args <- commandArgs(TRUE);
library(OmicCircos)

#Read in a GTF/BED to get a good apporoximation of geneome info (args[1])
seg = read.table(args[1],header=F, row.names=NULL,sep="\t")
if (grepl(".bed$", args[1], perl=T, ignore.case=T)) {
	seg = subset(seg, seg$V8=="gene")
	segs = data.frame(cbind(seg$V1,seg$V2,seg$V3,seg$V7,seg$V8))
} else { 
	seg = subset(seg,seg$V3=="gene")
	segs = data.frame(cbind(seg$V1,seg$V4,seg$V5,seg$V7,seg$V8))
}
colnames(segs) = c("chrom","chromStart","chromEnd","strand","name")
segs = subset(segs,substr(segs$chrom,1,3) %in% "chr")

#create f, which stores chromosome info
seg.name= unique(segs$chrom)
f = segs[1:length(seg.name),]
f$chrom=unique(segs$chrom)

for (i in 1:nrow(f)){
f[i,"chromStart"]=min(as.numeric(subset(segs,segs$chrom==f[i,"chrom"])[,"chromStart"]))
f[i,"chromEnd"]=max(as.numeric(subset(segs,segs$chrom==f[i,"chrom"])[,"chromEnd"]))
}
chrOrder<-c(paste("chr",1:22,sep=""),"chrX","chrY","chrM")
rownames(f)<-f$chrom
f<-f[chrOrder,]

db = segAnglePo(f,seg=f$chrom)
colors = rainbow(length(unique(seg.name)),alpha=0.5)

#Read in the Fusion Data (args[2])
fusions =read.table(args[2],header=F,row.names=NULL,sep="\t")
fusions.1<-do.call(rbind.data.frame, strsplit(fusions[,1], split=":"))
fusions.2<-do.call(rbind.data.frame, strsplit(fusions[,2], split=":"))
fusions<-cbind(fusions.1, fusions.2, fusions[,c(-1,-2)])
fusions<-merge(fusions, cbind(chrOrder, colors), by.x=1, by.y="chrOrder")
fusions[,1]<-factor(fusions[,1], levels=chrOrder)
fusions[,4]<-factor(fusions[,4], levels=chrOrder)

#Restrict to just Interchromosomal Fusions
#fusions.inter<-fusions[which(fusions[,1] != fusions[,4]),]
fusions.inter<-fusions


#Genes to Label
#seg$GeneID<-gsub(".*gene_name:", "", seg$V4) 
#seg$GeneID<-gsub(";.*", "", seg$GeneID)
#Most commonly found genes
#common_genes<-names(sort(table(c(fusions[,11], fusions[,13])), decreasing=T))[1:40]
#common_gene_locations<-(seg[which(seg$GeneID %in% common_genes),c("V1", "V2", "GeneID")])
#fusions2<-fusions[which(fusions[,11] %in% common_genes | fusions[,13] %in% common_genes),]


par ( mar=c ( 2 , 2 , 2 , 2 ) ) ;
pdf(paste(args[3], ".pdf", sep=""))
plot(c(1,800),c(1,800),type="n",axes=F,xlab="",ylab="",main="")
#circos(R=310,cir=db, W=20, side="out", cex=0.25, type="label2", mapping=labelGenes_dat_noRP, col=c("black"))
circos(R=300,cir=db,scale=F, type="chr",side="in", col=rainbow(nrow(f),alpha=1.0),print.chr.lab=T, W=4)
circos(R=280,cir=db, mapping=fusions.inter,type="link",lwd=0.8, col=fusions.inter$colors)
dev.off()

