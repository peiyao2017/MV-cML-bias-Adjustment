setwd("/panfs/jay/groups/20/panwei/wan01299/collider_bias/WHR")

library(qqman)
library(data.table)
library(tidyr)
library(dplyr)
library(R.utils)
Ncov=10
filename=numeric()
for(i in 1:Ncov){
  filename[i]=paste("resultMAInfo",i,".txt",sep="")
}
outputname=numeric()
for(i in 1:Ncov){
  outputname[i]=paste('Res_ldblockM2_',i,'.RData',sep="")
}
 

result1=read.table("/panfs/jay/groups/20/panwei/wan01299/collider_bias/M2/WHR/coefficient/responseM2.WHR.glm.linear",header = FALSE)
result_met=data.frame(chr=result1$V1,BP=result1$V2,rs=result1$V3,beta=result1$V9,se=result1$V10,p=result1$V12)
colnames(result_met)=c("chr","BP","rs","beta","se","p")
bed=read.table(file="fourier_ls-all.bed.txt",header=TRUE,sep="\t")
bed$chr=as.numeric(substr(bed$chr,4,5))
sig_before=result_met[result_met$p<=5e-08,]
bed$id=rep(0,nrow(bed))
for(i in 1:nrow(bed)){
  bed$id[i]=i
}


for(u in 1:Ncov){
  print(u)
  result2=read.table(file=filename[u],header = TRUE)
  result2$BP=result_met$BP
  result_adj=data.frame(chr=result2$chr,BP=result2$BP,rs=result2$rs,beta=result2$beta_adj,se=result2$sd_adj,p=result2$p_adj)
  colnames(result_adj)=c("chr","BP","rs","beta","se","p")
  sig_after=result_adj[result_adj$p<=5e-08,]
  validation=read.table(file="WHR_validation.txt",sep="\t",header=TRUE)
  validation=na.omit(validation)
  if(nrow(sig_after)>0){
    before=matrix(0,nrow=nrow(sig_before),ncol=nrow(bed))
    after=matrix(0,nrow=nrow(sig_after),ncol=nrow(bed))
    valid=matrix(0,nrow=nrow(validation),ncol=nrow(bed))
    BP_before=sig_before$BP
    BP_after=sig_after$BP
    chr_before=sig_before$chr
    chr_after=sig_after$chr
    start=bed$start
    stop=bed$stop
    id=bed$id
    chr=bed$chr
    BP_valid=validation$BP
    chr_valid=validation$chr
    
    
    for(i in 1:nrow(sig_before)){
      for(j in 1:nrow(bed)){
        if(chr_before[i]==chr[j]&(BP_before[i]>=start[j]&BP_before[i]<=stop[j])){
          before[i,j]=id[j]
        }
      }
    }
    
    
    for(i in 1:nrow(sig_after)){
      for(j in 1:nrow(bed)){
        if(chr_after[i]==chr[j]&(BP_after[i]>=start[j]&BP_after[i]<=stop[j])){
          after[i,j]=id[j]
        }
      }
    }
    
    
    for(i in 1:nrow(valid)){
      for(j in 1:nrow(bed)){
        if(chr_valid[i]==chr[j]&(BP_valid[i]>=start[j]&BP_valid[i]<=stop[j])){
          valid[i,j]=id[j]
        }
      }
    }
    
    
    a1=numeric()
    a2=numeric()
    a3=numeric()
    
    for(i in 1:nrow(before)){
      a1=c(a1,before[i,][before[i,]!=0])
    }
    
    for(i in 1:nrow(after)){
      a2=c(a2,after[i,][after[i,]!=0])
    }
    
    
    for(i in 1:nrow(valid)){
      a3=c(a3,valid[i,][valid[i,]!=0])
    }
    
    
    a1=unique(na.omit(a1))
    a2=unique(na.omit(a2))
    a3=unique(na.omit(a3))
    
    
    Nsnps_before=nrow(sig_before)
    before=length(a1)
    validation_small_before=NA
    validation_large_before=sum(a1 %in% a3)
    
    Nsnps_after=nrow(sig_after)
    after=length(a2)
    validation_small_after=NA
    validation_large_after=sum(a2 %in% a3)
    
    
    Res = list(Nsnps_before=Nsnps_before,
               before=before,
               validation_small_before=validation_small_before,
               validation_large_before=validation_large_before,
               Nsnps_after=Nsnps_after,
               after=after,
               validation_small_after=validation_small_after,
               validation_large_after=validation_large_after)
    save(Res,file=outputname[u])
    print(u)
  }
  if(nrow(sig_after)==0){
    before=matrix(0,nrow=nrow(sig_before),ncol=nrow(bed))
    valid=matrix(0,nrow=nrow(validation),ncol=nrow(bed))
    BP_before=sig_before$BP
    BP_after=sig_after$BP
    chr_before=sig_before$chr
    chr_after=sig_after$chr
    start=bed$start
    stop=bed$stop
    id=bed$id
    chr=bed$chr
    BP_valid=validation$BP
    chr_valid=validation$chr
    
    
    for(i in 1:nrow(sig_before)){
      for(j in 1:nrow(bed)){
        if(chr_before[i]==chr[j]&(BP_before[i]>=start[j]&BP_before[i]<=stop[j])){
          before[i,j]=id[j]
        }
      }
    }
    
    
    
    
    for(i in 1:nrow(valid)){
      for(j in 1:nrow(bed)){
        if(chr_valid[i]==chr[j]&(BP_valid[i]>=start[j]&BP_valid[i]<=stop[j])){
          valid[i,j]=id[j]
        }
      }
    }
    
    
    a1=numeric()
    a3=numeric()
    
    for(i in 1:nrow(before)){
      a1=c(a1,before[i,][before[i,]!=0])
    }
    
    
    
    for(i in 1:nrow(valid)){
      a3=c(a3,valid[i,][valid[i,]!=0])
    }
    
    
    a1=unique(na.omit(a1))
    a3=unique(na.omit(a3))
    
    
    Nsnps_before=nrow(sig_before)
    before=length(a1)
    validation_small_before=NA
    validation_large_before=sum(a1 %in% a3)
    
    Nsnps_after=0
    after=length=0
    validation_small_after=NA
    validation_large_after=0
    
    
    Res = list(Nsnps_before=Nsnps_before,
               before=before,
               validation_small_before=validation_small_before,
               validation_large_before=validation_large_before,
               Nsnps_after=Nsnps_after,
               after=after,
               validation_small_after=validation_small_after,
               validation_large_after=validation_large_after)
    save(Res,file=outputname[u])
    print(u)
  }
}


