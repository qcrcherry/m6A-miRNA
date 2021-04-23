library(ggplot2)

dir="C:\\Users\\scikuangren\\Desktop\\miRNA\\9_Riskscore"
setwd(dir)
data<-read.table("risk_score.txt",header = T,sep = "\t",check.names = F)

ggplot(data,aes(x=risk_level,y=risk_score,fill=risk_level))+geom_boxplot()+xlab("Risk")+ylab("Risk Score")+labs(fill="Risk")


ggplot(data,aes(x=risk_level,y=risk_score,fill=risk_level))+geom_dotplot(binaxis = "y",binwidth = 0.2,stackdir = "center")+xlab("Risk")+ylab("Risk Score")+labs(fill="Risk")


attach(data)
pdf("3.pdf")
plot(risk_score,survival_time,col=ifelse(status==1,"red","green"),xlab = "",ylab="Survival time",pch=16)
legend("topright", c("Death", "Alive"), pch=16, col=c("red","green"))
dev.off()


ggplot(data,aes(x=risk_score,y=survival_time))+geom_area(aes(fill=risk_level))+geom_line()+geom_hline(yintercept = 0)+xlab("Risk Score")+ylab("Survival time")+labs(fill="Risk")



