
install.packages("ggplot2")

library(edgeR)
library(pheatmap)
library(ggplot2)


dir="C:\\Users\\scikuangren\\Desktop\\miRNA\\1_diff"
setwd(dir)
getwd()
tcga<-read.table("miRNA.txt",header = T,sep = "\t",check.names = F)
tcga=as.matrix(tcga)
rownames(tcga)=tcga[,1]
GeneExp=tcga[,2:ncol(tcga)]
TCGA=matrix(as.numeric(as.matrix(GeneExp)),nrow=nrow(GeneExp),dimnames=list(rownames(GeneExp),colnames(GeneExp)))
TCGA=avereps(TCGA)
TCGA=TCGA[rowMeans(TCGA)>1,]
design=c(rep("normal",3),rep("cancer",162))                       
mydesign <- model.matrix(~design)
mydgelist <- DGEList(counts=TCGA,group=design)
mydgelist<- calcNormFactors(mydgelist)
mydgelist<- estimateCommonDisp(mydgelist)
mydgelist <- estimateTagwiseDisp(mydgelist,trend = "movingave")
mytest <- exactTest(mydgelist,pair = c("normal","cancer"))
Allgene<-topTags(mytest,n=10000000)
Allgene=Allgene$table
iddata<-mydgelist$pseudo.counts

write.table(Allgene,"Allgene.xls",sep="\t",quote = F)
Diffgene = Allgene[(Allgene$FDR < 0.05 & (Allgene$logFC>2 | Allgene$logFC<(-2))),]
write.table(Diffgene, "Diffgene.xls",sep="\t",quote=F)
Upgene = Allgene[(Allgene$FDR < 0.05 & (Allgene$logFC>2)),]
write.table(Upgene, "Upgene.xls",sep="\t",quote=F)
Downgene = Allgene[(Allgene$FDR < 0.05 & (Allgene$logFC<(-2))),]
write.table(Downgene, "Downgene.xls",sep="\t",quote=F)
Normalizegeneexp=rbind(id=colnames(iddata),iddata)
write.table(Normalizegeneexp,"Normalizegeneexp.txt",sep="\t",quote=F,col.names=F)   
Diffgeneexp=rbind(id=colnames(iddata),iddata[rownames(Diffgene),])
write.table(Diffgeneexp,"Diffgeneexp.txt",sep="\t",quote=F,col.names=F)

heatmap=m6A_data
rownames(heatmap)=newGeneLists
heatmap=log2(heatmap+1)
Type=c(rep("N",conNum),rep("T",treatNum))  
names(Type)=colnames(heatmap)
Type=as.data.frame(Type)
png(file="heatmap.png",width = 2400,   height =1800,res = 300)
pheatmap(heatmap, 
         annotation=Type, 
         color = colorRampPalette(c("green", "white", "red"))(50),
         cluster_cols =F,
         show_colnames = F,
         scale="row",
         fontsize = 10,
         fontsize_row=8,
         fontsize_col=3)
dev.off()