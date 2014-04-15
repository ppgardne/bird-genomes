#R CMD BATCH --no-save scripts/heatmaps.R

library(gplots)
breaks<-c(-1000,0,1,5,10,50,100,500,1000)
lb<-length(breaks)-2

lmat = rbind(c(1,2),c(3,4))
lwid = rbind(c(4,.5))
lhei = rbind(c(4,.25))

key<-c(0,1,5,10,50,100,500,1000)
key<-matrix(c(key,key),ncol=2,nrow=8,dimnames = list(key,c("","")), byrow = FALSE)
pdf(file="paper/figures/key.pdf",width=9, height=9)
heatmap.2(key,trace="none",col=gray(lb:0/lb),
          dendrogram = "none", scale = "none",
          breaks=breaks,
          sepcolor="black", Rowv=F, Colv=F, key = F,
          lmat=lmat,
          cexCol=1.1, cexRow=1.1
          )
dev.off()

rna<-read.table("data/R/RNA.dat",header = T, sep = "\t",row.names=52)
pdf(file="paper/figures/RNA.pdf",width=8, height=16)
heatmap.2(as.matrix(rna),trace="none",col=gray(lb:0/lb),
          dendrogram = "none", scale = "none",
          breaks=breaks,
          sepcolor="black", Rowv=F, Colv=F, key = F,
          lmat=lmat, lwid=lwid, lhei=lhei,
          cexCol=1.1, cexRow=1.1
          )
dev.off()



mirna<-read.table("data/R/miRNA.dat",header = T, sep = "\t",row.names=52)
pdf(file="paper/figures/miRNA.pdf",width=8, height=50)
heatmap.2(as.matrix(mirna),trace="none",col=gray(lb:0/lb),
          dendrogram = "none", scale = "none",
          breaks=breaks,
          sepcolor="black", Rowv=F, Colv=F, key = F,
          lmat=lmat, lwid=lwid, lhei=lhei,
          cexCol=1.1, cexRow=1.1
          )
dev.off()


snorna<-read.table("data/R/snoRNA.dat",header = T, sep = "\t",row.names=52)
pdf(file="paper/figures/snoRNA.pdf",width=8, height=25)
heatmap.2(as.matrix(snorna),trace="none",col=gray(lb:0/lb),
          dendrogram = "none", scale = "none",
          breaks=breaks,
          sepcolor="black", Rowv=F, Colv=F, key = F,
          lmat=lmat, lwid=lwid, lhei=lhei,
          cexCol=1.1, cexRow=1.1
          )
dev.off()

snorna<-read.table("data/R/snoRNA-human-yeast-correspondences.dat",header = T, sep = "\t",row.names=52)
pdf(file="paper/figures/snoRNA-human-yeast-correspondences.pdf",width=8, height=16)
heatmap.2(as.matrix(snorna),trace="none",col=gray(lb:0/lb),
          dendrogram = "none", scale = "none",
          breaks=breaks,
          sepcolor="black", Rowv=F, Colv=F, key = F,
          lmat=lmat, lwid=lwid, lhei=lhei,
          cexCol=1.1, cexRow=1.1
          )
dev.off()

trna<-read.table("data/R/tRNA.dat",header = T, sep = "\t",row.names=52)
pdf(file="paper/figures/tRNA.pdf",width=8, height=8)
heatmap.2(as.matrix(trna),trace="none",col=gray(lb:0/lb),
          dendrogram = "none", scale = "none",
          breaks=breaks,
          sepcolor="black", Rowv=F, Colv=F, key = F,
          lmat=lmat, lwid=lwid, lhei=lhei,
          cexCol=0.8, cexRow=0.8
          )
dev.off()

lncrna<-read.table("data/R/lncRNA.dat",header = T, sep = "\t",row.names=52)
pdf(file="paper/figures/lncRNA.pdf",width=8, height=8)
heatmap.2(as.matrix(lncrna),trace="none",col=gray(lb:0/lb),
          dendrogram = "none", scale = "none",
          breaks=breaks,
          sepcolor="black", Rowv=F, Colv=F, key = F,
          lmat=lmat, lwid=lwid, lhei=lhei,
          cexCol=0.8, cexRow=0.8
          )
dev.off()

allrna<-read.table("data/R/allRNA.dat",header = T, sep = "\t",row.names=52)
pdf(file="paper/figures/allRNA.pdf",width=8, height=80)
heatmap.2(as.matrix(allrna),trace="none",col=gray(lb:0/lb),
          dendrogram = "none", scale = "none",
          breaks=breaks,
          sepcolor="black", Rowv=F, Colv=F, key = F,
          lmat=lmat, lwid=lwid, lhei=lhei,
          cexCol=1.1, cexRow=1.1
          )
dev.off()

allrna<-read.table("data/R/diverged.dat",header = T, sep = "\t",row.names=52)
pdf(file="paper/figures/diverged.pdf",width=25, height=30)
heatmap.2(as.matrix(allrna),trace="none",col=gray(lb:0/lb),
          dendrogram = "none", scale = "none",
          breaks=breaks,
          sepcolor="black", Rowv=F, Colv=F, key = F,
          lmat=lmat, lwid=lwid, lhei=lhei,
          cexCol=3.0, cexRow=3.0, margins=c(15, 5)
          )
par(mar=c(21, 0, 0, 18) + 0.1, oma=c(0,0,0,0))
#xvals<-c(8.9,10.6,  58.5,60.2,  65.4,67.1)
xvals<-c(8.9,10.6,  65.4,67.1)
for (x in xvals){
    lines(c(x,x)/100, c(0,2), lwd=5)
}
box("plot", col="black",lwd=3)  


dev.off()

allrna<-read.table("data/R/unusual-conserved.dat",header = T, sep = "\t",row.names=52)
pdf(file="paper/figures/unusual-conserved.pdf",width=25, height=30)

heatmap.2(as.matrix(allrna),trace="none",col=gray(lb:0/lb),
          dendrogram = "none", scale = "none",
          breaks=breaks,
          sepcolor="black", Rowv=F, Colv=F, key = F,
          lmat=lmat, lwid=lwid, lhei=lhei,
          cexCol=3.0, cexRow=3.0, margins=c(15, 5)
          )
dev.off()

######################################################################

allrna<-read.table("data/R/high-copy-numbers.dat",header = T, sep = "\t",row.names=52)
pdf(file="paper/figures/high-copy-numbers.pdf",width=8, height=10)
op<-par(las=2,mar=c(7, 4, 4, 2) + 0.1,cex=1.5); #bottom, left, top, right
plot(log10(allrna$human), pch=1, yaxt="n", xaxt="n", ylab="Genomic copy number", xlab="", main="Genomic copy numbers for high-copy ncRNAs"); 
axis(2,0:3, 10^(0:3)); 
axis(1,1:(length(rownames(allrna))), rownames(allrna)); 
points(log10(allrna$chicken), pch=2, col="red"); 
points(log10(allrna$alligator), pch=3); 
points(log10(allrna$turtle), pch=4);  
points(log10(allrna$zebrafinch), pch=6, col="red"); 
points(log10(allrna$budgerigar), pch=7, col="red"); 
points(log10(rowMeans(allrna[,4:51])),col="red",pch=5); 
legend("topright", c("Human", "Alligator", "Turtle", "Chicken", "Zebrafinch", "Budgerigar", "mean.48.Birds"), pch=c(1,3,4,2,6,7,5), col=c("black","black","black","red","red","red","red"))
dev.off()

humU6<-read.table("data/R/U6-human-bitscores.dat",header = F, sep = "\t")
chickU6<-read.table("data/R/U6-chicken-bitscores.dat",header = F, sep = "\t")
humSRP<-read.table("data/R/SRP-human-bitscores.dat",header = F, sep = "\t")
chickSRP<-read.table("data/R/SRP-chicken-bitscores.dat",header = F, sep = "\t")
humY<-read.table("data/R/Y_RNA-human-bitscores.dat",header = F, sep = "\t")
chickY<-read.table("data/R/Y_RNA-chicken-bitscores.dat",header = F, sep = "\t")

allU6<-c(humU6$V1,chickU6$V1)
breaks<-seq(floor(min(allU6)),ceiling(max(allU6)), len=40)
hU6<-hist(humU6$V1,plot=F,breaks=breaks)
cU6<-hist(chickU6$V1,plot=F,breaks=breaks)

allSRP<-c(humSRP$V1,chickSRP$V1)
breaks<-seq(floor(min(allSRP)),ceiling(max(allSRP)), len=40)
hSRP<-hist(humSRP$V1,plot=F,breaks=breaks)
cSRP<-hist(chickSRP$V1,plot=F,breaks=breaks)

allY<-c(humY$V1,chickY$V1)
breaks<-seq(floor(min(allY)),ceiling(max(allY)), len=40)
hY<-hist(humY$V1,plot=F,breaks=breaks)
cY<-hist(chickY$V1,plot=F,breaks=breaks)

pdf(file="paper/figures/high-copy-numbers-U6-SRP.pdf",width=8, height=10)
op<-par(las=1,cex=1.5,mfcol=c(3,1))
plot(hU6$mids[hU6$counts>0],hU6$counts[hU6$counts>0],ylim=c(0,max(hU6$counts[hU6$counts>0])),col="black",pch=1,main="U6 bitscore distributions",xlab="bit-score",ylab="Freq.")
lines(hU6$mids[hU6$counts>0],hU6$counts[hU6$counts>0],col="black")
points(cU6$mids[cU6$counts>0],cU6$counts[cU6$counts>0],col="red",pch=2)
lines(cU6$mids[cU6$counts>0],cU6$counts[cU6$counts>0],col="red")

plot(hSRP$mids[hSRP$counts>0],hSRP$counts[hSRP$counts>0],col="black",pch=1,main="SRP bitscore distributions",xlab="bit-score",ylab="Freq.")
lines(hSRP$mids[hSRP$counts>0],hSRP$counts[hSRP$counts>0],col="black")
points(cSRP$mids[cSRP$counts>0],cSRP$counts[cSRP$counts>0],col="red",pch=2)
lines(cSRP$mids[cSRP$counts>0],cSRP$counts[cSRP$counts>0],col="red")
legend("topright", c("Human", "Chicken"), pch=c(1,2), col=c("black","red"))

plot(hY$mids[hY$counts>0],hY$counts[hY$counts>0],col="black",pch=1,main="Y RNA bitscore distributions",xlab="bit-score",ylab="Freq.")
lines(hY$mids[hY$counts>0],hY$counts[hY$counts>0],col="black")
points(cY$mids[cY$counts>0],cY$counts[cY$counts>0],col="red",pch=2)
lines(cY$mids[cY$counts>0],cY$counts[cY$counts>0],col="red")
legend("topright", c("Human", "Chicken"), pch=c(1,2), col=c("black","red"))

dev.off()

######################################################################
rnaExp   <-read.table("data/RNA-seq/summed-vals.dat",header = F, sep = "\t")
rnaExpIndMc<-read.csv2("data/RNA-seq/mccarthy_expression_tissue.dat")
rnaExpIndUl<-read.csv2("data/RNA-seq/ulitsky_expression_tissue.dat")

len<-0;
c<-0; 
tf<-(rnaExpIndMc[,2]>2)
thresh<-0.25
for (i in 2:16) {
    c<-c(c,as.numeric(rnaExpIndMc[,i]))
    len<-len+length(rnaExpIndMc[rnaExpIndMc[,i]>thresh,i])
    tf<-(tf | rnaExpIndMc[,i]>thresh)
}
for (i in 2:7 ) {
    c<-c(c,as.numeric(rnaExpIndUl[,i]))
    tf<-(tf | rnaExpIndMc[,i]>thresh)
}
length(rnaExpIndMc[tf,1])/length(rnaExpIndMc[,1])
length(rnaExp$V3[rnaExp$V3>10])/length(rnaExp$V3)

breaks<- seq(0, max(log10(rnaExp$V3+1)), length.out=200)


pdf(file="paper/figures/expression-distribution.pdf",width=5, height=5)
par(mfcol=c(2,1))
hist(log10(rnaExp$V3+1),breaks=breaks,xaxt = "n",xlab="log10(sum(RNA)+1)", main="Summed RNA over all tissues")
axis(1, 0:6, 10^(0:6))
hist(log10(c+1),        breaks=breaks,xaxt = "n",xlab="log10(RNA+1)",      main="RNA over all tissues")
axis(1, 0:6, 10^(0:6))
dev.off()

length(rnaExp$V3[rnaExp$V3>10])



