Ncov=10
library(data.table)
library(tidyr)
library(dplyr)
library(R.utils)
library(hudson)
setwd("D:/BMI")
if(Ncov>0){
  N_manplots=Ncov
  Man_plot_name=numeric()
  for(i in 1:Ncov){
    Man_plot_name[i]=paste("D:/BMI/manhattanBMIfirstturn",i,sep="")
  }
  main_before=numeric()
  main_after=numeric()
  for(i in 1:Ncov){
    if(i>1){
      main_before[i]=paste("BMI, ",20, "covariates are included")
      main_after[i]=paste(i, "covariates are adjusted")
    }
    if(i==1){
      main_before[i]=paste("BMI, ",20, "covariate are included")
      main_after[i]=paste(i, "covariate is adjusted")
    }
  }
  filename_after=numeric()
  for(i in 1:Ncov){
    filename_after[i]=paste("resultMAInfo",i,".txt",sep="")
  }
  result1=read.table(file="response.BMI.glm.linear",header = FALSE)
  result_met=data.frame(CHR=result1$V1,POS=result1$V2,SNP=result1$V3,beta=result1$V9,se=result1$V10,p=result1$V12)
  head(result_met)
  result_met$p1=result_met$p
  for(i in 1:nrow(result_met)){
    if(result_met$p[i]==0){
      result_met$p1[i]=min(result_met$p[result_met$p!=0])+abs(rnorm(1,sd=5e-16))
    }
  }
  result_met$pvalue=result_met$p1  
  A=list()
  for(j in 1:Ncov){
    result2=read.table(file = filename_after[j],header = TRUE)
    result2$BP=result_met$POS
    result_adj=data.frame(CHR=result2$chr,POS=result2$BP,SNP=result2$rs,beta=result2$beta_adj,se=result2$sd_adj,p=result2$p_adj)
    head(result_adj)
    result_adj$p1=result_adj$p
    for(i in 1:nrow(result_adj)){
      if(result_adj$p[i]==0){
        result_adj$p1[i]=min(result_adj$p[result_adj$p!=0])+abs(rnorm(1,sd=5e-16))
      }
    }
    result_adj$pvalue=result_adj$p1
    A[[j]]=result_adj
    print(j)
  }
  
  for(i in 1:Ncov){
    gmirror(top=A[[i]],bottom = result_met,tline=5e-8,bline=5e-8,toptitle =main_after[i],bottomtitle = main_before[i],file = Man_plot_name[i],res=100)
  }
}
  
  
setwd("D:/HC")
if(Ncov>0){
  N_manplots=Ncov
  Man_plot_name=numeric()
  for(i in 1:Ncov){
    Man_plot_name[i]=paste("D:/HC/manhattanHCfirstturn",i,sep="")
  }
  main_before=numeric()
  main_after=numeric()
  for(i in 1:Ncov){
    if(i>1){
      main_before[i]=paste("HC, ",20, "covariates are included")
      main_after[i]=paste(i, "covariates are adjusted")
    }
    if(i==1){
      main_before[i]=paste("HC, ",20, "covariate are included")
      main_after[i]=paste(i, "covariate is adjusted")
    }
  }
  filename_after=numeric()
  for(i in 1:Ncov){
    filename_after[i]=paste("resultMAInfo",i,".txt",sep="")
  }
  result1=read.table(file="response.HC.glm.linear",header = FALSE)
  result_met=data.frame(CHR=result1$V1,POS=result1$V2,SNP=result1$V3,beta=result1$V9,se=result1$V10,p=result1$V12)
  head(result_met)
  result_met$p1=result_met$p
  for(i in 1:nrow(result_met)){
    if(result_met$p[i]==0){
      result_met$p1[i]=min(result_met$p[result_met$p!=0])+abs(rnorm(1,sd=5e-16))
    }
  }
  result_met$pvalue=result_met$p1  
  A=list()
  for(j in 1:Ncov){
    result2=read.table(file = filename_after[j],header = TRUE)
    result2$BP=result_met$POS
    result_adj=data.frame(CHR=result2$chr,POS=result2$BP,SNP=result2$rs,beta=result2$beta_adj,se=result2$sd_adj,p=result2$p_adj)
    head(result_adj)
    result_adj$p1=result_adj$p
    for(i in 1:nrow(result_adj)){
      if(result_adj$p[i]==0){
        result_adj$p1[i]=min(result_adj$p[result_adj$p!=0])+abs(rnorm(1,sd=5e-16))
      }
    }
    result_adj$pvalue=result_adj$p1
    A[[j]]=result_adj
    print(j)
  }
  
  for(i in 1:Ncov){
    gmirror(top=A[[i]],bottom = result_met,tline=5e-8,bline=5e-8,toptitle =main_after[i],bottomtitle = main_before[i],file = Man_plot_name[i],res=100)
  }
}


setwd("D:/WC")
if(Ncov>0){
  N_manplots=Ncov
  Man_plot_name=numeric()
  for(i in 1:Ncov){
    Man_plot_name[i]=paste("D:/WC/manhattanWCfirstturn",i,sep="")
  }
  main_before=numeric()
  main_after=numeric()
  for(i in 1:Ncov){
    if(i>1){
      main_before[i]=paste("WC, ",20, "covariates are included")
      main_after[i]=paste(i, "covariates are adjusted")
    }
    if(i==1){
      main_before[i]=paste("WC, ",20, "covariate are included")
      main_after[i]=paste(i, "covariate is adjusted")
    }
  }
  filename_after=numeric()
  for(i in 1:Ncov){
    filename_after[i]=paste("resultMAInfo",i,".txt",sep="")
  }
  result1=read.table(file="response.WC.glm.linear",header = FALSE)
  result_met=data.frame(CHR=result1$V1,POS=result1$V2,SNP=result1$V3,beta=result1$V9,se=result1$V10,p=result1$V12)
  head(result_met)
  result_met$p1=result_met$p
  for(i in 1:nrow(result_met)){
    if(result_met$p[i]==0){
      result_met$p1[i]=min(result_met$p[result_met$p!=0])+abs(rnorm(1,sd=5e-16))
    }
  }
  result_met$pvalue=result_met$p1  
  A=list()
  for(j in 1:Ncov){
    result2=read.table(file = filename_after[j],header = TRUE)
    result2$BP=result_met$POS
    result_adj=data.frame(CHR=result2$chr,POS=result2$BP,SNP=result2$rs,beta=result2$beta_adj,se=result2$sd_adj,p=result2$p_adj)
    head(result_adj)
    result_adj$p1=result_adj$p
    for(i in 1:nrow(result_adj)){
      if(result_adj$p[i]==0){
        result_adj$p1[i]=min(result_adj$p[result_adj$p!=0])+abs(rnorm(1,sd=5e-16))
      }
    }
    result_adj$pvalue=result_adj$p1
    A[[j]]=result_adj
    print(j)
  }
  
  for(i in 1:Ncov){
    gmirror(top=A[[i]],bottom = result_met,tline=5e-8,bline=5e-8,toptitle =main_after[i],bottomtitle = main_before[i],file = Man_plot_name[i],res=100)
  }
}


setwd("D:/WHR")
if(Ncov>0){
  N_manplots=Ncov
  Man_plot_name=numeric()
  for(i in 1:Ncov){
    Man_plot_name[i]=paste("D:/WHR/manhattanWHRfirstturn",i,sep="")
  }
  main_before=numeric()
  main_after=numeric()
  for(i in 1:Ncov){
    if(i>1){
      main_before[i]=paste("WHR, ",20, "covariates are included")
      main_after[i]=paste(i, "covariates are adjusted")
    }
    if(i==1){
      main_before[i]=paste("WHR, ",20, "covariate are included")
      main_after[i]=paste(i, "covariate is adjusted")
    }
  }
  filename_after=numeric()
  for(i in 1:Ncov){
    filename_after[i]=paste("resultMAInfo",i,".txt",sep="")
  }
  result1=read.table(file="response.WHR.glm.linear",header = FALSE)
  result_met=data.frame(CHR=result1$V1,POS=result1$V2,SNP=result1$V3,beta=result1$V9,se=result1$V10,p=result1$V12)
  head(result_met)
  result_met$p1=result_met$p
  for(i in 1:nrow(result_met)){
    if(result_met$p[i]==0){
      result_met$p1[i]=min(result_met$p[result_met$p!=0])+abs(rnorm(1,sd=5e-16))
    }
  }
  result_met$pvalue=result_met$p1  
  A=list()
  for(j in 1:Ncov){
    result2=read.table(file = filename_after[j],header = TRUE)
    result2$BP=result_met$POS
    result_adj=data.frame(CHR=result2$chr,POS=result2$BP,SNP=result2$rs,beta=result2$beta_adj,se=result2$sd_adj,p=result2$p_adj)
    head(result_adj)
    result_adj$p1=result_adj$p
    for(i in 1:nrow(result_adj)){
      if(result_adj$p[i]==0){
        result_adj$p1[i]=min(result_adj$p[result_adj$p!=0])+abs(rnorm(1,sd=5e-16))
      }
    }
    result_adj$pvalue=result_adj$p1
    A[[j]]=result_adj
    print(j)
  }
  
  for(i in 1:Ncov){
    gmirror(top=A[[i]],bottom = result_met,tline=5e-8,bline=5e-8,toptitle =main_after[i],bottomtitle = main_before[i],file = Man_plot_name[i],res=100)
  }
}


setwd("D:/Height")
if(Ncov>0){
  N_manplots=Ncov
  Man_plot_name=numeric()
  for(i in 1:Ncov){
    Man_plot_name[i]=paste("D:/Height/manhattanHeightfirstturn",i,sep="")
  }
  main_before=numeric()
  main_after=numeric()
  for(i in 1:Ncov){
    if(i>1){
      main_before[i]=paste("Height, ",20, "covariates are included")
      main_after[i]=paste(i, "covariates are adjusted")
    }
    if(i==1){
      main_before[i]=paste("Height, ",20, "covariate are included")
      main_after[i]=paste(i, "covariate is adjusted")
    }
  }
  filename_after=numeric()
  for(i in 1:Ncov){
    filename_after[i]=paste("resultMAInfo",i,".txt",sep="")
  }
  result1=read.table(file="response.Height.glm.linear",header = FALSE)
  result_met=data.frame(CHR=result1$V1,POS=result1$V2,SNP=result1$V3,beta=result1$V9,se=result1$V10,p=result1$V12)
  head(result_met)
  result_met$p1=result_met$p
  for(i in 1:nrow(result_met)){
    if(result_met$p[i]==0){
      result_met$p1[i]=min(result_met$p[result_met$p!=0])+abs(rnorm(1,sd=5e-16))
    }
  }
  result_met$pvalue=result_met$p1  
  A=list()
  for(j in 1:Ncov){
    result2=read.table(file = filename_after[j],header = TRUE)
    result2$BP=result_met$POS
    result_adj=data.frame(CHR=result2$chr,POS=result2$BP,SNP=result2$rs,beta=result2$beta_adj,se=result2$sd_adj,p=result2$p_adj)
    head(result_adj)
    result_adj$p1=result_adj$p
    for(i in 1:nrow(result_adj)){
      if(result_adj$p[i]==0){
        result_adj$p1[i]=min(result_adj$p[result_adj$p!=0])+abs(rnorm(1,sd=5e-16))
      }
    }
    result_adj$pvalue=result_adj$p1
    A[[j]]=result_adj
    print(j)
  }
  
  for(i in 1:Ncov){
    gmirror(top=A[[i]],bottom = result_met,tline=5e-8,bline=5e-8,toptitle =main_after[i],bottomtitle = main_before[i],file = Man_plot_name[i],res=100)
  }
}

setwd("D:/SBP")
if(Ncov>0){
  N_manplots=Ncov
  Man_plot_name=numeric()
  for(i in 1:Ncov){
    Man_plot_name[i]=paste("D:/SBP/manhattanSBPfirstturn",i,sep="")
  }
  main_before=numeric()
  main_after=numeric()
  for(i in 1:Ncov){
    if(i>1){
      main_before[i]=paste("SBP, ",20, "covariates are included")
      main_after[i]=paste(i, "covariates are adjusted")
    }
    if(i==1){
      main_before[i]=paste("SBP, ",20, "covariate are included")
      main_after[i]=paste(i, "covariate is adjusted")
    }
  }
  filename_after=numeric()
  for(i in 1:Ncov){
    filename_after[i]=paste("resultMAInfo",i,".txt",sep="")
  }
  result1=read.table(file="response.SBP.glm.linear",header = FALSE)
  result_met=data.frame(CHR=result1$V1,POS=result1$V2,SNP=result1$V3,beta=result1$V9,se=result1$V10,p=result1$V12)
  head(result_met)
  result_met$p1=result_met$p
  for(i in 1:nrow(result_met)){
    if(result_met$p[i]==0){
      result_met$p1[i]=min(result_met$p[result_met$p!=0])+abs(rnorm(1,sd=5e-16))
    }
  }
  result_met$pvalue=result_met$p1  
  A=list()
  for(j in 1:Ncov){
    result2=read.table(file = filename_after[j],header = TRUE)
    result2$BP=result_met$POS
    result_adj=data.frame(CHR=result2$chr,POS=result2$BP,SNP=result2$rs,beta=result2$beta_adj,se=result2$sd_adj,p=result2$p_adj)
    head(result_adj)
    result_adj$p1=result_adj$p
    for(i in 1:nrow(result_adj)){
      if(result_adj$p[i]==0){
        result_adj$p1[i]=min(result_adj$p[result_adj$p!=0])+abs(rnorm(1,sd=5e-16))
      }
    }
    result_adj$pvalue=result_adj$p1
    A[[j]]=result_adj
    print(j)
  }
  
  for(i in 1:Ncov){
    gmirror(top=A[[i]],bottom = result_met,tline=5e-8,bline=5e-8,toptitle =main_after[i],bottomtitle = main_before[i],file = Man_plot_name[i],res=100)
  }
}

setwd("D:/DBP")
if(Ncov>0){
  N_manplots=Ncov
  Man_plot_name=numeric()
  for(i in 1:Ncov){
    Man_plot_name[i]=paste("D:/DBP/manhattanDBPfirstturn",i,sep="")
  }
  main_before=numeric()
  main_after=numeric()
  for(i in 1:Ncov){
    if(i>1){
      main_before[i]=paste("DBP, ",20, "covariates are included")
      main_after[i]=paste(i, "covariates are adjusted")
    }
    if(i==1){
      main_before[i]=paste("DBP, ",20, "covariate are included")
      main_after[i]=paste(i, "covariate is adjusted")
    }
  }
  filename_after=numeric()
  for(i in 1:Ncov){
    filename_after[i]=paste("resultMAInfo",i,".txt",sep="")
  }
  result1=read.table(file="response.DBP.glm.linear",header = FALSE)
  result_met=data.frame(CHR=result1$V1,POS=result1$V2,SNP=result1$V3,beta=result1$V9,se=result1$V10,p=result1$V12)
  head(result_met)
  result_met$p1=result_met$p
  for(i in 1:nrow(result_met)){
    if(result_met$p[i]==0){
      result_met$p1[i]=min(result_met$p[result_met$p!=0])+abs(rnorm(1,sd=5e-16))
    }
  }
  result_met$pvalue=result_met$p1  
  A=list()
  for(j in 1:Ncov){
    result2=read.table(file = filename_after[j],header = TRUE)
    result2$BP=result_met$POS
    result_adj=data.frame(CHR=result2$chr,POS=result2$BP,SNP=result2$rs,beta=result2$beta_adj,se=result2$sd_adj,p=result2$p_adj)
    head(result_adj)
    result_adj$p1=result_adj$p
    for(i in 1:nrow(result_adj)){
      if(result_adj$p[i]==0){
        result_adj$p1[i]=min(result_adj$p[result_adj$p!=0])+abs(rnorm(1,sd=5e-16))
      }
    }
    result_adj$pvalue=result_adj$p1
    A[[j]]=result_adj
    print(j)
  }
  
  for(i in 1:Ncov){
    gmirror(top=A[[i]],bottom = result_met,tline=5e-8,bline=5e-8,toptitle =main_after[i],bottomtitle = main_before[i],file = Man_plot_name[i],res=100)
  }
}


