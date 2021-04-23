
dir="C:\\Users\\scikuangren\\Desktop\\miRNA\\2_UnivariateCoxData"
setwd(dir)
inputfile1="clinical_data.txt" 
inputfile2="miRNA_data.txt" 

time_data<-read.table(inputfile1,header = T,sep = "\t",check.names = F)
geneEXP<-read.table(inputfile2,header = T,sep = "\t",check.names = F)
merger_data<-merge(time_data,geneEXP,by="id")
write.table(merger_data,"merger_data.txt",sep = "\t",row.names = F,quote = F)