library(data.table)

setwd("C:/Users/wangp/Downloads/HC")
data=read.table(file="GIANT_2015_HIP_COMBINED_EUR.txt",sep="\t",header = TRUE)
data_adj=read.table(file="GIANT_2015_HIPadjBMI_COMBINED_EUR.txt",header = TRUE,sep="\t")
data=na.omit(data)
data_adj=na.omit(data_adj)
data=data[data$p<=5e-8,]
rsid=data$MarkerName
rsid_adj=data_adj$MarkerName
for(i in 1:length(rsid)){
  tar=unlist(gregexpr(":",rsid[i]))[1]
  if(tar>0){
    rsid[i]=paste("rs",substr(rsid[i],start=tar+1,stop=length(rsid[i])),sep="")
  }
  print(i)
}

rsid %in% rsid_adj
BP=numeric()
chr=numeric()
for(i in 1:length(rsid)){
  BP[i]=data_adj[data_adj$MarkerName==rsid[i],]$Pos
  chr[i]=data_adj[data_adj$MarkerName==rsid[i],]$Chr
  print(i)
}
HC_validation_large=data.frame(BP=BP,chr=chr)
HC_validation_large=na.omit(HC_validation_large)
write.table(HC_validation_large,file="HC_validation_large.txt",sep="\t",row.names = FALSE,col.names = TRUE,quote=FALSE)