
dir="C:\\Users\\scikuangren\\Desktop\\miRNA\\4_lassoData"
setwd(dir)
inputfile1="id.txt" 
inputfile2="input.txt" 

time_data<-read.table(inputfile1,header = T,sep = "\t",check.names = F)
geneEXP<-read.table(inputfile2,header = T,sep = "\t",check.names = F)
merger_data<-merge(time_data,geneEXP,by="id")
write.table(merger_data,"output.txt",sep = "\t",row.names = F,quote = F)