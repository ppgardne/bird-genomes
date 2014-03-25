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



