library(survival)

dir="C:\\Users\\scikuangren\\Desktop\\miRNA\\12_Survival"
setwd(dir)
inputfile2="multiCoxinput.txt"

inputmiRNA<- read.table(inputfile2,header=T,sep="\t",check.names = F,row.names = 1)

inputmiRNA$survival_time=inputmiRNA$survival_time/365



for(a in colnames(inputmiRNA[,3:ncol(inputmiRNA)])){
  exp<-factor(ifelse(inputmiRNA[,a]>median(inputmiRNA[,a]),"high","low"))
  kms<-survfit(Surv(survival_time,status)~exp,data=inputmiRNA)
  kmdffexp=survdiff(Surv(survival_time,status)~exp,data=inputmiRNA)
  pValue=round(1-pchisq(kmdffexp$chisq,df=1),8)
    myname=unlist(strsplit(a,"\\|"))[1]
    pdffile=paste(myname,".pdf")
  pdf(file = pdffile,10,8)
  plot(kms,lty=1.8,col=c("red","green"),
       xlab="Survival time in years",ylab="Survival probabilities",
       main=paste(myname,"(P=", pValue ,")",sep=""))
  legend("bottomleft",c("High expression","Low expression"),lty=1.8,col=c("red","green"))
  dev.off()
  }

inputfileRisk="risk_score.txt"

inputRisk<- read.table(inputfileRisk,header=T,sep="\t",check.names = F,row.names = 1)
kms2<-survfit(Surv(survival_time,status)~risk_level,data=inputRisk)
kmdffexp1=survdiff(Surv(survival_time,status)~risk_level,data=inputRisk)
pValue1=round(1-pchisq(kmdffexp$chisq,df=1),8)

pdf("risk.pdf")
plot(kms2,lty=1.8,col=c("red","green"),
     xlab="Survival time in years",ylab="Survival probabilities",
     main=paste("Risk score","(P=", pValue1 ,")",sep=""))
legend("bottomleft",c("High risk","Low risk"),lty=1.8,col=c("red","green"))
dev.off()
