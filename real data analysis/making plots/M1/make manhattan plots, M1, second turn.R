Ncov=10
library(data.table)
library(tidyr)
library(dplyr)
library(R.utils)
library(hudson)





setwd("D:/BMI")
if(Ncov>0){
  subtitle=numeric()
  N_manplots=Ncov
  Man_plot_name=numeric()
  for(i in 1:Ncov){
    Man_plot_name[i]=paste("D:/BMI/manhattanBMIsecondturn",i,sep="")
  }
  main_before=numeric()
  main_after=numeric()
  for(i in 1:Ncov){
    if(i>1){
      main_before[i]=paste("BMI, ",i, "covariates are included")
      main_after[i]=paste(i, "covariates are adjusted")
    }
    if(i==1){
      main_before[i]=paste("BMI, ",i, "covariate is included")
      main_after[i]=paste(i, "covariate is adjusted")
    }
  }
  filename_before=numeric()
  filename_after=numeric()
  for(i in 1:Ncov){
    filename_before[i]=paste("response_",i,".BMI.glm.linear",sep="")
    filename_after[i]=paste("resultMAInfo_second_turn",i,".txt",sep="")
  }
  A=list()
  B=list()
  for(j in 1:Ncov){
    result1=read.table(file=filename_before[j],header = FALSE)
    result_met=data.frame(CHR=result1$V1,POS=result1$V2,SNP=result1$V3,beta=result1$V9,se=result1$V10,p=result1$V12)
    head(result_met)
    result_met$p1=result_met$p
    for(i in 1:nrow(result_met)){
      if(result_met$p[i]==0){
        result_met$p1[i]=min(result_met$p[result_met$p!=0])+abs(rnorm(1,sd=5e-16))
      }
    }
    result_met$pvalue=result_met$p1  
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
    B[[j]]=result_met
    print(j)
  }
  
  for(i in 1:Ncov){
    gmirror(top=A[[i]],bottom = B[[i]],tline=5e-8,bline=5e-8,toptitle =main_after[i],bottomtitle = main_before[i],file = Man_plot_name[i],res=100)
  }
}

setwd("D:/WC")
if(Ncov>0){
  subtitle=numeric()
  N_manplots=Ncov
  Man_plot_name=numeric()
  for(i in 1:Ncov){
    Man_plot_name[i]=paste("D:/WC/manhattanWCsecondturn",i,sep="")
  }
  main_before=numeric()
  main_after=numeric()
  for(i in 1:Ncov){
    if(i>1){
      main_before[i]=paste("WC, ",i, "covariates are included")
      main_after[i]=paste(i, "covariates are adjusted")
    }
    if(i==1){
      main_before[i]=paste("WC, ",i, "covariate is included")
      main_after[i]=paste(i, "covariate is adjusted")
    }
  }
  filename_before=numeric()
  filename_after=numeric()
  for(i in 1:Ncov){
    filename_before[i]=paste("response_",i,".WC.glm.linear",sep="")
    filename_after[i]=paste("resultMAInfo_second_turn",i,".txt",sep="")
  }
  A=list()
  B=list()
  for(j in 1:Ncov){
    result1=read.table(file=filename_before[j],header = FALSE)
    result_met=data.frame(CHR=result1$V1,POS=result1$V2,SNP=result1$V3,beta=result1$V9,se=result1$V10,p=result1$V12)
    head(result_met)
    result_met$p1=result_met$p
    for(i in 1:nrow(result_met)){
      if(result_met$p[i]==0){
        result_met$p1[i]=min(result_met$p[result_met$p!=0])+abs(rnorm(1,sd=5e-16))
      }
    }
    result_met$pvalue=result_met$p1  
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
    B[[j]]=result_met
    print(j)
  }
  
  for(i in 1:Ncov){
    gmirror(top=A[[i]],bottom = B[[i]],tline=5e-8,bline=5e-8,toptitle =main_after[i],bottomtitle = main_before[i],file = Man_plot_name[i],res=100)
  }
}

setwd("D:/Height")
if(Ncov>0){
  subtitle=numeric()
  N_manplots=Ncov
  Man_plot_name=numeric()
  for(i in 1:Ncov){
    Man_plot_name[i]=paste("D:/Height/manhattanHeightsecondturn",i,sep="")
  }
  main_before=numeric()
  main_after=numeric()
  for(i in 1:Ncov){
    if(i>1){
      main_before[i]=paste("Height, ",i, "covariates are included")
      main_after[i]=paste(i, "covariates are adjusted")
    }
    if(i==1){
      main_before[i]=paste("Height, ",i, "covariate is included")
      main_after[i]=paste(i, "covariate is adjusted")
    }
  }
  filename_before=numeric()
  filename_after=numeric()
  for(i in 1:Ncov){
    filename_before[i]=paste("response_",i,".Height.glm.linear",sep="")
    filename_after[i]=paste("resultMAInfo_second_turn",i,".txt",sep="")
  }
  A=list()
  B=list()
  for(j in 1:Ncov){
    result1=read.table(file=filename_before[j],header = FALSE)
    result_met=data.frame(CHR=result1$V1,POS=result1$V2,SNP=result1$V3,beta=result1$V9,se=result1$V10,p=result1$V12)
    head(result_met)
    result_met$p1=result_met$p
    for(i in 1:nrow(result_met)){
      if(result_met$p[i]==0){
        result_met$p1[i]=min(result_met$p[result_met$p!=0])+abs(rnorm(1,sd=5e-16))
      }
    }
    result_met$pvalue=result_met$p1  
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
    B[[j]]=result_met
    print(j)
  }
  
  for(i in 1:Ncov){
    gmirror(top=A[[i]],bottom = B[[i]],tline=5e-8,bline=5e-8,toptitle =main_after[i],bottomtitle = main_before[i],file = Man_plot_name[i],res=100)
  }
}

setwd("D:/WHR")
if(Ncov>0){
  subtitle=numeric()
  N_manplots=Ncov
  Man_plot_name=numeric()
  for(i in 1:Ncov){
    Man_plot_name[i]=paste("D:/WHR/manhattanWHRsecondturn",i,sep="")
  }
  main_before=numeric()
  main_after=numeric()
  for(i in 1:Ncov){
    if(i>1){
      main_before[i]=paste("WHR, ",i, "covariates are included")
      main_after[i]=paste(i, "covariates are adjusted")
    }
    if(i==1){
      main_before[i]=paste("WHR, ",i, "covariate is included")
      main_after[i]=paste(i, "covariate is adjusted")
    }
  }
  filename_before=numeric()
  filename_after=numeric()
  for(i in 1:Ncov){
    filename_before[i]=paste("response_",i,".WHR.glm.linear",sep="")
    filename_after[i]=paste("resultMAInfo_second_turn",i,".txt",sep="")
  }
  A=list()
  B=list()
  for(j in 1:Ncov){
    result1=read.table(file=filename_before[j],header = FALSE)
    result_met=data.frame(CHR=result1$V1,POS=result1$V2,SNP=result1$V3,beta=result1$V9,se=result1$V10,p=result1$V12)
    head(result_met)
    result_met$p1=result_met$p
    for(i in 1:nrow(result_met)){
      if(result_met$p[i]==0){
        result_met$p1[i]=min(result_met$p[result_met$p!=0])+abs(rnorm(1,sd=5e-16))
      }
    }
    result_met$pvalue=result_met$p1  
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
    B[[j]]=result_met
    print(j)
  }
  
  for(i in 1:Ncov){
    gmirror(top=A[[i]],bottom = B[[i]],tline=5e-8,bline=5e-8,toptitle =main_after[i],bottomtitle = main_before[i],file = Man_plot_name[i],res=100)
  }
}


setwd("D:/HC")
if(Ncov>0){
  subtitle=numeric()
  N_manplots=Ncov
  Man_plot_name=numeric()
  for(i in 1:Ncov){
    Man_plot_name[i]=paste("D:/HC/manhattanHCsecondturn",i,sep="")
  }
  main_before=numeric()
  main_after=numeric()
  for(i in 1:Ncov){
    if(i>1){
      main_before[i]=paste("HC, ",i, "covariates are included")
      main_after[i]=paste(i, "covariates are adjusted")
    }
    if(i==1){
      main_before[i]=paste("HC, ",i, "covariate is included")
      main_after[i]=paste(i, "covariate is adjusted")
    }
  }
  filename_before=numeric()
  filename_after=numeric()
  for(i in 1:Ncov){
    filename_before[i]=paste("response_",i,".HC.glm.linear",sep="")
    filename_after[i]=paste("resultMAInfo_second_turn",i,".txt",sep="")
  }
  A=list()
  B=list()
  for(j in 1:Ncov){
    result1=read.table(file=filename_before[j],header = FALSE)
    result_met=data.frame(CHR=result1$V1,POS=result1$V2,SNP=result1$V3,beta=result1$V9,se=result1$V10,p=result1$V12)
    head(result_met)
    result_met$p1=result_met$p
    for(i in 1:nrow(result_met)){
      if(result_met$p[i]==0){
        result_met$p1[i]=min(result_met$p[result_met$p!=0])+abs(rnorm(1,sd=5e-16))
      }
    }
    result_met$pvalue=result_met$p1  
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
    B[[j]]=result_met
    print(j)
  }
  
  for(i in 1:Ncov){
    gmirror(top=A[[i]],bottom = B[[i]],tline=5e-8,bline=5e-8,toptitle =main_after[i],bottomtitle = main_before[i],file = Man_plot_name[i],res=100)
  }
}


setwd("D:/SBP")
if(Ncov>0){
  subtitle=numeric()
  N_manplots=Ncov
  Man_plot_name=numeric()
  for(i in 1:Ncov){
    Man_plot_name[i]=paste("D:/SBP/manhattanSBPsecondturn",i,sep="")
  }
  main_before=numeric()
  main_after=numeric()
  for(i in 1:Ncov){
    if(i>1){
      main_before[i]=paste("SBP, ",i, "covariates are included")
      main_after[i]=paste(i, "covariates are adjusted")
    }
    if(i==1){
      main_before[i]=paste("SBP, ",i, "covariate is included")
      main_after[i]=paste(i, "covariate is adjusted")
    }
  }
  filename_before=numeric()
  filename_after=numeric()
  for(i in 1:Ncov){
    filename_before[i]=paste("response_",i,".SBP.glm.linear",sep="")
    filename_after[i]=paste("resultMAInfo_second_turn",i,".txt",sep="")
  }
  A=list()
  B=list()
  for(j in 1:Ncov){
    result1=read.table(file=filename_before[j],header = FALSE)
    result_met=data.frame(CHR=result1$V1,POS=result1$V2,SNP=result1$V3,beta=result1$V9,se=result1$V10,p=result1$V12)
    head(result_met)
    result_met$p1=result_met$p
    for(i in 1:nrow(result_met)){
      if(result_met$p[i]==0){
        result_met$p1[i]=min(result_met$p[result_met$p!=0])+abs(rnorm(1,sd=5e-16))
      }
    }
    result_met$pvalue=result_met$p1  
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
    B[[j]]=result_met
    print(j)
  }
  
  for(i in 1:Ncov){
    gmirror(top=A[[i]],bottom = B[[i]],tline=5e-8,bline=5e-8,toptitle =main_after[i],bottomtitle = main_before[i],file = Man_plot_name[i],res=100)
  }
}


setwd("D:/DBP")
if(Ncov>0){
  subtitle=numeric()
  N_manplots=Ncov
  Man_plot_name=numeric()
  for(i in 1:Ncov){
    Man_plot_name[i]=paste("D:/DBP/manhattanDBPsecondturn",i,sep="")
  }
  main_before=numeric()
  main_after=numeric()
  for(i in 1:Ncov){
    if(i>1){
      main_before[i]=paste("DBP, ",i, "covariates are included")
      main_after[i]=paste(i, "covariates are adjusted")
    }
    if(i==1){
      main_before[i]=paste("DBP, ",i, "covariate is included")
      main_after[i]=paste(i, "covariate is adjusted")
    }
  }
  filename_before=numeric()
  filename_after=numeric()
  for(i in 1:Ncov){
    filename_before[i]=paste("response_",i,".DBP.glm.linear",sep="")
    filename_after[i]=paste("resultMAInfo_second_turn",i,".txt",sep="")
  }
  A=list()
  B=list()
  for(j in 1:Ncov){
    result1=read.table(file=filename_before[j],header = FALSE)
    result_met=data.frame(CHR=result1$V1,POS=result1$V2,SNP=result1$V3,beta=result1$V9,se=result1$V10,p=result1$V12)
    head(result_met)
    result_met$p1=result_met$p
    for(i in 1:nrow(result_met)){
      if(result_met$p[i]==0){
        result_met$p1[i]=min(result_met$p[result_met$p!=0])+abs(rnorm(1,sd=5e-16))
      }
    }
    result_met$pvalue=result_met$p1  
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
    B[[j]]=result_met
    print(j)
  }
  
  for(i in 1:Ncov){
    gmirror(top=A[[i]],bottom = B[[i]],tline=5e-8,bline=5e-8,toptitle =main_after[i],bottomtitle = main_before[i],file = Man_plot_name[i],res=100)
  }
}