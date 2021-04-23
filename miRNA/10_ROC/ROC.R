install.packages("timeROC")


library(survival)
library(timeROC)

dir="C:\\Users\\scikuangren\\Desktop\\miRNA\\10_ROC"
setwd(dir)
miRNA<-read.table("risk_score.txt",header=T,sep="\t")
predict_3_year<- 3
predict_5_year<- 5

ROC<-timeROC(T=miRNA$survival_time,delta=miRNA$status,
                  marker=miRNA$risk_score,cause=1,
                  weighting="marginal",
                  times=c(predict_3_year,predict_5_year),ROC=TRUE)

pdf("ROC.pdf")
plot(ROC,time=predict_3_year,title=FALSE,lwd=3)
plot(ROC,time=predict_5_year,col="blue",add=TRUE,title=FALSE,lwd=3)
legend("bottomright",
       c(paste("AUC of 3 year survival: ",round(ROC$AUC[1],3)),
         paste("AUC of 5 year survival: ",round(ROC$AUC[2],3))),col=c("red","blue"),lwd=3)
dev.off()

