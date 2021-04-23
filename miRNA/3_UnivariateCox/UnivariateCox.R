

library(survival)

dir="C:\\Users\\scikuangren\\Desktop\\miRNA\\3_UnivariateCox"

setwd(dir)

inputfile="merger_data.txt"

miRNA<-read.table(inputfile,header=T,sep="\t",row.names = 1,check.names = F)

miRNAEXP=log2(miRNA[,3:ncol(miRNA)]+1)
miRNA=cbind(miRNA[,1:2],miRNAEXP)

coxR=data.frame()
coxf<-function(x){
fmla1 <- as.formula(Surv(survival_time,status)~miRNA[,x])
mycox <- coxph(fmla1,data=miRNA)
}
for(a in colnames(miRNA[,3:ncol(miRNA)])){
mycox=coxf(a)
coxResult = summary(mycox)
coxR=rbind(coxR,cbind(miRNAname=a,HR=coxResult$coefficients[,"exp(coef)"],
P=coxResult$coefficients[,"Pr(>|z|)"]))
}
write.table(coxR,"coxResult.txt",sep="\t",row.names=F,quote=F)




