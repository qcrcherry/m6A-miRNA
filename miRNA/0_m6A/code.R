

rt=read.table(ESCA-count.xls),sep='\t',header=T,check.names = F)
rt=na.omit(rt)
exp=as.matrix(rt[-1])
exp=avereps(exp)
exp=2^exp-1
exp=exp[,order(substr(colnames(exp),14,15),decreasing = T)]
m6A_gene=read.table('m6A_gene.txt',sep='\t',check.names = F)
m6A_gene=m6A_gene$V1
m6A_gene=m6A_gene[m6A_gene%in%rownames(exp)]
m6A_data=exp[m6A_gene,]
rownames(m6A_data)[4]='KIAA1429'
write.table(m6A_data,file="m6Aexp.txt",sep="\t",quote=F,col.names = NA)

table(substr(colnames(m6A_data),14,15)<10)
outTab=data.frame()
newGeneLists=c()
conNum=table(substr(colnames(m6A_data),14,15)<10)[1]
treatNum=table(substr(colnames(m6A_data),14,15)<10)[2]

#################################################

grade=c(rep(1,conNum),rep(2,treatNum)) 
for(i in row.names(m6A_data)){
  geneName=unlist(strsplit(i,"\\|",))[1]
  geneName=gsub("\\/", "_", geneName)
  rt=rbind(expression=m6A_data[i,],grade=grade)
  rt=as.matrix(t(rt))
  wilcoxTest<-wilcox.test(expression ~ grade, data=rt)
  conGeneMeans=mean(m6A_data[i,1:conNum])
  treatGeneMeans=mean(m6A_data[i,(conNum+1):ncol(m6A_data)])
  logFC=log2(treatGeneMeans)-log2(conGeneMeans)  
  pvalue=wilcoxTest$p.value
  conMed=median(m6A_data[i,1:conNum])
  treatMed=median(m6A_data[i,(conNum+1):ncol(m6A_data)])
  diffMed=treatMed-conMed
  outTab=rbind(outTab,cbind(gene=i,conMean=conGeneMeans,treatMean=treatGeneMeans,logFC=logFC,pValue=pvalue))
  if(pvalue<0.001){
    newGeneLists=c(newGeneLists,paste0(i,"***"))
  }else if(pvalue<0.01){
    newGeneLists=c(newGeneLists,paste0(i,"**"))
  }else if(pvalue<0.05){
    newGeneLists=c(newGeneLists,paste0(i,"*"))
  }else{
    newGeneLists=c(newGeneLists,i)
  }
}

write.table(outTab,file="diff.xls",sep="\t",row.names=F,quote=F)

heatmap=m6A_data
rownames(heatmap)=newGeneLists
heatmap=log2(heatmap+1)
Type=c(rep("N",conNum),rep("T",treatNum))   
names(Type)=colnames(heatmap)
Type=as.data.frame(Type)
pdf("heatmap.pdf",height=5,width=10)
pheatmap(heatmap, 
         annotation=Type, 
         color = colorRampPalette(c("green", "white", "red"))(50),
         cluster_cols =F,
         show_colnames = F,
         scale="row",
         fontsize = 10,
         fontsize_row=10,
         fontsize_col=3)
dev.off()


#install.packages("vioplot")

viomap=t(m6A_data)
pdf("vioplot.pdf",height=6,width=9)           
par(las=1,mar=c(4,5,3,3))
x=c(1:ncol(viomap))
y=c(1:ncol(viomap))
plot(x,y,
     xlim=c(0,37),ylim=c(min(viomap),max(viomap)*1.1),
     main="",xlab="", ylab="Gene expression",
     pch=21,
     cex.lab=1.5,
     col="white",
     xaxt="n")

for(i in 1:ncol(viomap)){
  normalData=viomap[1:conNum,i]
  tumorData=viomap[(conNum+1):(conNum+treatNum),i]
  vioplot(normalData,at=3*(i-1),lty=1,add = T,col = 'blue')
  vioplot(tumorData,at=3*(i-1)+1,lty=1,add = T,col = 'red')
  wilcoxTest=wilcox.test(normalData,tumorData)
  p=round(wilcoxTest$p.value,3)
  mx=max(c(normalData,tumorData))
  lines(c(x=3*(i-1)+0.2,x=3*(i-1)+0.8),c(mx,mx))
  text(x=3*(i-1)+0.5, y=mx+1.5, labels=ifelse(p<0.001, paste0("p<0.001"), paste0("p=",p)), cex = 0.8)
}
text(seq(1,37,3),-5,xpd = NA,labels=colnames(viomap),cex = 1,srt = 45,pos=2)
dev.off()
