library(survival)

dir="C:\\Users\\scikuangren\\Desktop\\miRNA\\7_MultivariateCox"
setwd(dir)
inputfile="multiCoxinput.txt"

miRNA<-read.table(inputfile,header=T,sep="\t",row.names = 1,check.names = F,stringsAsFactors = F) 

miRNAEXP=log2(miRNA[,3:ncol(miRNA)]+1)
miRNA=cbind(miRNA[,1:2],miRNAEXP)
miRNA[,"survival_time"]=miRNA[,"survival_time"]/365
fmla1 <- as.formula(Surv(survival_time,status)~.)
mycox <- coxph(fmla1,data=miRNA)
risk_score<-predict(mycox,type="risk",newdata=miRNA)
risk_level<-as.factor(ifelse(risk_score>median(risk_score),"High","Low"))
write.table(cbind(id=rownames(cbind(miRNA[,1:2],risk_score,risk_level)),cbind(miRNA[,1:2],risk_score,risk_level)),"risk_score.txt",sep="\t",quote=F,row.names=F)

summary(mycox)


install.packages("survminer")

library(survminer)

pdf("forest1.pdf",12,8)
ggforest(mycox,fontsize = 1)
dev.off()

