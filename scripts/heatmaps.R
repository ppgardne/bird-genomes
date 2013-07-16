#R CMD BATCH --no-save heatmaps.R

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



