setwd("/panfs/jay/groups/20/panwei/wan01299/collider_bias/M2/BMI/coefficient/")
data=read.table(file="BMI_all_cov.txt",sep="\t",header=TRUE)
names=numeric()
for(i in 1:20){
  names[i]=paste("MET",i,"_score.profile",sep="")
}
a=list()
for(i in 1:20){
  a[[i]]=read.table(file=names[i],header = TRUE,sep="")
}
for(i in 1:20){
  data[,14+i]=scale(data[,14+i],center = TRUE,scale = TRUE)-a[[i]][,ncol(a[[i]])]
}
write.table(data,file="BMI_all_cov_NoGene.txt",sep="\t",col.names=TRUE,quote=FALSE,row.names = FALSE)



setwd("/panfs/jay/groups/20/panwei/wan01299/collider_bias/M2/WHR/coefficient/")
data=read.table(file="WHR_all_cov.txt",sep="\t",header=TRUE)
names=numeric()
for(i in 1:20){
  names[i]=paste("MET",i,"_score.profile",sep="")
}
a=list()
for(i in 1:20){
  a[[i]]=read.table(file=names[i],header = TRUE,sep="")
}
for(i in 1:20){
  data[,14+i]=scale(data[,14+i],center = TRUE,scale = TRUE)-a[[i]][,ncol(a[[i]])]
}
write.table(data,file="WHR_all_cov_NoGene.txt",sep="\t",col.names=TRUE,quote=FALSE,row.names = FALSE)



setwd("/panfs/jay/groups/20/panwei/wan01299/collider_bias/M2/WC/coefficient/")
data=read.table(file="WC_all_cov.txt",sep="\t",header=TRUE)
names=numeric()
for(i in 1:20){
  names[i]=paste("MET",i,"_score.profile",sep="")
}
a=list()
for(i in 1:20){
  a[[i]]=read.table(file=names[i],header = TRUE,sep="")
}
for(i in 1:20){
  data[,14+i]=scale(data[,14+i],center = TRUE,scale = TRUE)-a[[i]][,ncol(a[[i]])]
}
write.table(data,file="WC_all_cov_NoGene.txt",sep="\t",col.names=TRUE,quote=FALSE,row.names = FALSE)



setwd("/panfs/jay/groups/20/panwei/wan01299/collider_bias/M2/HC/coefficient/")
data=read.table(file="HC_all_cov.txt",sep="\t",header=TRUE)
names=numeric()
for(i in 1:20){
  names[i]=paste("MET",i,"_score.profile",sep="")
}
a=list()
for(i in 1:20){
  a[[i]]=read.table(file=names[i],header = TRUE,sep="")
}
for(i in 1:20){
  data[,14+i]=scale(data[,14+i],center = TRUE,scale = TRUE)-a[[i]][,ncol(a[[i]])]
}
write.table(data,file="HC_all_cov_NoGene.txt",sep="\t",col.names=TRUE,quote=FALSE,row.names = FALSE)




setwd("/panfs/jay/groups/20/panwei/wan01299/collider_bias/M2/Height/coefficient/")
data=read.table(file="Height_all_cov.txt",sep="\t",header=TRUE)
names=numeric()
for(i in 1:20){
  names[i]=paste("MET",i,"_score.profile",sep="")
}
a=list()
for(i in 1:20){
  a[[i]]=read.table(file=names[i],header = TRUE,sep="")
}
for(i in 1:20){
  data[,14+i]=scale(data[,14+i],center = TRUE,scale = TRUE)-a[[i]][,ncol(a[[i]])]
}
write.table(data,file="Height_all_cov_NoGene.txt",sep="\t",col.names=TRUE,quote=FALSE,row.names = FALSE)



setwd("/panfs/jay/groups/20/panwei/wan01299/collider_bias/M2/SBP/coefficient/")
data=read.table(file="SBP_all_cov.txt",sep="\t",header=TRUE)
names=numeric()
for(i in 1:20){
  names[i]=paste("MET",i,"_score.profile",sep="")
}
a=list()
for(i in 1:20){
  a[[i]]=read.table(file=names[i],header = TRUE,sep="")
}
for(i in 1:20){
  data[,14+i]=scale(data[,14+i],center = TRUE,scale = TRUE)-a[[i]][,ncol(a[[i]])]
}
write.table(data,file="SBP_all_cov_NoGene.txt",sep="\t",col.names=TRUE,quote=FALSE,row.names = FALSE)

setwd("/panfs/jay/groups/20/panwei/wan01299/collider_bias/M2/DBP/coefficient/")
data=read.table(file="DBP_all_cov.txt",sep="\t",header=TRUE)
names=numeric()
for(i in 1:20){
  names[i]=paste("MET",i,"_score.profile",sep="")
}
a=list()
for(i in 1:20){
  a[[i]]=read.table(file=names[i],header = TRUE,sep="")
}
for(i in 1:20){
  data[,14+i]=scale(data[,14+i],center = TRUE,scale = TRUE)-a[[i]][,ncol(a[[i]])]
}
write.table(data,file="DBP_all_cov_NoGene.txt",sep="\t",col.names=TRUE,quote=FALSE,row.names = FALSE)
