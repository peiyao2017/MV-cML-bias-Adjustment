setwd("C:/Users/wangp/Downloads/gwas_sim/BMI")
data=read.table(file="BMI_all_cov.txt",sep="\t",header=TRUE)
MET_20PC=cbind(data[,1:2],data[,15:ncol(data)])
noMET_20PC=cbind(data[,1:2],data[,3:14])
fwrite(MET_20PC,file='MET_20PC.txt',sep="\t")
fwrite(noMET_20PC,file='noMET_20PC.txt',sep="\t")


setwd("C:/Users/wangp/Downloads/gwas_sim/HC")
data=read.table(file="HC_all_cov.txt",sep="\t",header=TRUE)
MET_20PC=cbind(data[,1:2],data[,15:ncol(data)])
noMET_20PC=cbind(data[,1:2],data[,3:14])
fwrite(MET_20PC,file='MET_20PC.txt',sep="\t")
fwrite(noMET_20PC,file='noMET_20PC.txt',sep="\t")


setwd("C:/Users/wangp/Downloads/gwas_sim/WHR")
data=read.table(file="WHR_all_cov.txt",sep="\t",header=TRUE)
MET_20PC=cbind(data[,1:2],data[,15:ncol(data)])
noMET_20PC=cbind(data[,1:2],data[,3:14])
fwrite(MET_20PC,file='MET_20PC.txt',sep="\t")
fwrite(noMET_20PC,file='noMET_20PC.txt',sep="\t")


setwd("C:/Users/wangp/Downloads/gwas_sim/WC")
data=read.table(file="WC_all_cov.txt",sep="\t",header=TRUE)
MET_20PC=cbind(data[,1:2],data[,15:ncol(data)])
noMET_20PC=cbind(data[,1:2],data[,3:14])
fwrite(MET_20PC,file='MET_20PC.txt',sep="\t")
fwrite(noMET_20PC,file='noMET_20PC.txt',sep="\t")

setwd("C:/Users/wangp/Downloads/gwas_sim/Height")
data=read.table(file="Height_all_cov.txt",sep="\t",header=TRUE)
MET_20PC=cbind(data[,1:2],data[,15:ncol(data)])
noMET_20PC=cbind(data[,1:2],data[,3:14])
fwrite(MET_20PC,file='MET_20PC.txt',sep="\t")
fwrite(noMET_20PC,file='noMET_20PC.txt',sep="\t")



setwd("C:/Users/wangp/Downloads/gwas_sim/DBP")
data=read.table(file="DBP_all_cov.txt",sep="\t",header=TRUE)
MET_20PC=cbind(data[,1:2],data[,15:ncol(data)])
noMET_20PC=cbind(data[,1:2],data[,3:14])
fwrite(MET_20PC,file='MET_20PC.txt',sep="\t")
fwrite(noMET_20PC,file='noMET_20PC.txt',sep="\t")


setwd("C:/Users/wangp/Downloads/gwas_sim/SBP")
data=read.table(file="SBP_all_cov.txt",sep="\t",header=TRUE)
MET_20PC=cbind(data[,1:2],data[,15:ncol(data)])
noMET_20PC=cbind(data[,1:2],data[,3:14])
fwrite(MET_20PC,file='MET_20PC.txt',sep="\t")
fwrite(noMET_20PC,file='noMET_20PC.txt',sep="\t")






