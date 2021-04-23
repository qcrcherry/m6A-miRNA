


library(glmnet)
library(survival)
dir="C:\\Users\\scikuangren\\Desktop\\miRNA\\5_lasso"

setwd(dir)

inputfile="lasso_input.txt"

miRNA<-read.table(inputfile,header=T,sep="\t",row.names = 1,check.names = F,stringsAsFactors = F) 

miRNAEXP=log2(miRNA[,3:ncol(miRNA)]+1)
miRNA=cbind(miRNA[,1:2],miRNAEXP)
miRNA[,"survival_time"]=miRNA[,"survival_time"]/365
v1<-as.matrix(miRNA[,c(3:ncol(miRNA))])
v2 <- as.matrix(Surv(miRNA$survival_time,miRNA$status))

myfit <- glmnet(v1, v2, family = "cox")
pdf("lambda.pdf")
plot(myfit, xvar = "lambda", label = TRUE)
dev.off()

myfit2 <- cv.glmnet(v1, v2, family="cox")
pdf("min.pdf")
plot(myfit2)
abline(v=log(c(myfit2$lambda.min,myfit2$lambda.1se)),lty="dashed")
dev.off()

myfit2$lambda.min
coe <- coef(myfit, s = myfit2$lambda.min)
act_index <- which(coe != 0)
act_coe <- coe[act_index]
row.names(coe)[act_index]



