
for(d in c( 4 )){
  for(rho in c(0  )){
    if(10>0){
      N=20
      setwd("D:/art sim/New/no_pleiotropy/uvmr_mvmr")
      
      name1_1=numeric()
 
      filename1=paste("D:/art sim/New/no_pleiotropy/uvmr_mvmr/result/","d=",d,"betaXYnot0rho",rho,"Slope.txt",sep="")
      filename2=paste("D:/art sim/New/no_pleiotropy/uvmr_mvmr/result/","d=",d,"betaXYnot0rho",rho,"point_effect_estimates.RData",sep="") 
      for(i in 1:N){
        name1_1[i]=paste("sim",rho,"_",i,"betaXYnot",0,".txt",sep="")
      }
      result=list()
      for(i in 1:N){
        a1=read.table(file=name1_1[i],header = TRUE,sep="\t")
        result[[i]]=a1
      }
      for(i in 1:N){
        if(i==1){
          a=result[[i]]
        }
        if(i>1){
          a=rbind(a,result[[i]])
        }
      }
      all_snps=list()
      for(i in 1:1000){
        all_snps[[i]]=a[1:1000,]
      }
      for(i in 1:1000){
        for(j in  1:1000){
          all_snps[[i]][j,]=a[i+(j-1)*1000,]
        }
        print(i)
      }
      
      null_snps=list()
      for(i in 1:900){
        null_snps[[i]]=all_snps[[i+100]]
      }
      
      all_snps_no_y=list()
      for(i in c(1:950)){
        if(i<=50){
          all_snps_no_y[[i]]=all_snps[[i]]
        }
        if(i>50){
          all_snps_no_y[[i]]=all_snps[[i+50]]
        }
      }
      
      
      all_snps_H_no_y=list()
      for(i in 1:50){
        all_snps_H_no_y[[i]]=all_snps[[i]]
      }
      
      all_snps_y_only=list()
      for(i in 1:50){
        all_snps_y_only[[i]]=all_snps[[50+i]]
      }
       
      
      all_snps_y=list()
      for(i in 1:50){
        all_snps_y[[i]]=all_snps[[50+i]]
      }
      
    
      if(10>0){
        H_only=data.frame(beta_true=rep(0,50),beta_before=rep(0,50),sd_before=rep(0,50),se_before=rep(0,50),error_before=rep(0,50),
                          beta_cml=rep(0,50),sd_cml=rep(0,50),se_cml=rep(0,50),error_cml=rep(0,50),
                          beta_uv_cml=rep(0,50),sd_uv_cml=rep(0,50),se_uv_cml=rep(0,50),error_uv_cml=rep(0,50),
                          beta_egger=rep(0,50),sd_egger=rep(0,50),se_egger=rep(0,50),error_egger=rep(0,50),
                          beta_uv_egger=rep(0,50),sd_uv_egger=rep(0,50),se_uv_egger=rep(0,50),error_uv_egger=rep(0,50),
                          beta_uv_DD=rep(0,50),sd_uv_DD=rep(0,50),se_uv_DD=rep(0,50),error_uv_DD=rep(0,50),
                          beta_uv_SH=rep(0,50),sd_uv_SH=rep(0,50),se_uv_SH=rep(0,50),error_uv_SH=rep(0,50)
                          
                          
        )
        
        
        
        for(i1 in 1:50){
          H_only$beta_true[i1]=mean(all_snps[[i1]]$beta_true)
          H_only$beta_before[i1]=mean(all_snps[[i1]]$beta_before)
          H_only$sd_before[i1]=sd(all_snps[[i1]]$beta_before)
          H_only$se_before[i1]=mean(all_snps[[i1]]$se_before)
          H_only$error_before[i1]=mean(all_snps[[i1]]$p_before<0.05)
          H_only$se_error_before[i1]=sd(all_snps[[i1]]$p_before<0.05) 
          
          H_only$beta_cml[i1]=mean(all_snps[[i1]]$beta_after)
          H_only$sd_cml[i1]=sd(all_snps[[i1]]$beta_after)
          H_only$se_cml[i1]=mean(all_snps[[i1]]$se_after3)
          H_only$error_cml[i1]=mean(all_snps[[i1]]$p_after3<0.05)
          H_only$se_error_cml[i1]=sd(all_snps[[i1]]$p_after3<0.05)
          
          H_only$beta_uv_cml[i1]=mean(all_snps[[i1]]$beta_uv_cml)
          H_only$sd_uv_cml[i1]=sd(all_snps[[i1]]$beta_uv_cml)
          H_only$se_uv_cml[i1]=mean(all_snps[[i1]]$sd_uv_cml)
          H_only$error_uv_cml[i1]=mean(all_snps[[i1]]$p_uv_cml<0.05)
          H_only$se_error_uv_cml[i1]=sd(all_snps[[i1]]$p_uv_cml<0.05)
          
          H_only$beta_egger[i1]=mean(all_snps[[i1]]$beta_egger)
          H_only$sd_egger[i1]=sd(all_snps[[i1]]$beta_egger)
          H_only$se_egger[i1]=mean(all_snps[[i1]]$se_egger)
          H_only$error_egger[i1]=mean(all_snps[[i1]]$p_egger<0.05)
          H_only$se_error_egger[i1]=sd(all_snps[[i1]]$p_egger<0.05)
          
          
          H_only$beta_uv_egger[i1]=mean(all_snps[[i1]]$beta_uv_egger)
          H_only$sd_uv_egger[i1]=sd(all_snps[[i1]]$beta_uv_egger)
          H_only$se_uv_egger[i1]=mean(all_snps[[i1]]$sd_uv_egger)
          H_only$error_uv_egger[i1]=mean(all_snps[[i1]]$p_uv_egger<0.05)
          H_only$se_error_uv_egger[i1]=sd(all_snps[[i1]]$p_uv_egger<0.05)
          
          
          H_only$beta_uv_DD[i1]=mean(all_snps[[i1]]$beta_uv_DD)
          H_only$sd_uv_DD[i1]=sd(all_snps[[i1]]$beta_uv_DD)
          H_only$se_uv_DD[i1]=mean(all_snps[[i1]]$sd_uv_DD)
          H_only$error_uv_DD[i1]=mean(all_snps[[i1]]$p_uv_DD<0.05)
          H_only$se_error_uv_DD[i1]=sd(all_snps[[i1]]$p_uv_DD<0.05)
          
          
          H_only$beta_uv_SH[i1]=mean(all_snps[[i1]]$beta_uv_SH)
          H_only$sd_uv_SH[i1]=sd(all_snps[[i1]]$beta_uv_SH)
          H_only$se_uv_SH[i1]=mean(all_snps[[i1]]$sd_uv_SH)
          H_only$error_uv_SH[i1]=mean(all_snps[[i1]]$p_uv_SH<0.05)
          H_only$se_error_uv_SH[i1]=sd(all_snps[[i1]]$p_uv_SH<0.05)
          
        }
        
       
        round(H_only[sample(c(1:50),4,replace = FALSE),],2)
        
        b_true_name=numeric()
        b_egger_name=numeric()
        b_uv_egger_name=numeric()
        b_cml_name=numeric()
        b_uv_cml_name=numeric()
        b_DD_name=numeric()
        b_SH_name=numeric()
        se_b_true_name=numeric()
        se_b_egger_name=numeric()
        se_b_uv_egger_name=numeric()
        se_b_cml_name=numeric()
        se_b_uv_cml_name=numeric()
        se_b_DD_name=numeric()
        se_b_SH_name=numeric()
        
     
        for(i1 in 1:d){
          b_true_name[i1]=paste("true_slop",i1,sep="")
          b_egger_name[i1]=paste("egger_slop",i1,sep="")
          b_uv_egger_name[i1]=paste("egger_uv_slop",i1,sep="")
          b_cml_name[i1]=paste("cml_slop",i1,sep="")
          b_uv_cml_name[i1]=paste("cml_uv_slop",i1,sep="")
          b_DD_name[i1]=paste("DD_uv_slop",i1,sep="")
          b_SH_name[i1]=paste("SH_uv_slop",i1,sep="")
          se_b_true_name[i1]=paste("true_slop_se",i1,sep="")
          se_b_egger_name[i1]=paste("egger_slop_se",i1,sep="")
          se_b_uv_egger_name[i1]=paste("egger_uv_slop_se",i1,sep="")
          se_b_cml_name[i1]=paste("cml_slop_se",i1,sep="")
          se_b_uv_cml_name[i1]=paste("cml_uv_slop_se",i1,sep="")
          se_b_DD_name[i1]=paste("DD_uv_slop_se",i1,sep="")
          se_b_SH_name[i1]=paste("SH_uv_slop_se",i1,sep="")
        }
        
        slope_estimated=data.frame(b_cml=sapply(all_snps[[1]][b_cml_name],mean),b_cml_sd=sapply(all_snps[[1]][b_cml_name],sd),b_cml_se= sapply(all_snps[[1]][se_b_cml_name],mean) ,
                                   b_uv_cml=sapply(all_snps[[1]][b_uv_cml_name],mean),b_uv_cml_sd=sapply(all_snps[[1]][b_uv_cml_name],sd),b_uv_cml_se= sapply(all_snps[[1]][se_b_uv_cml_name],mean) ,
                                   b_egger=sapply(all_snps[[1]][b_egger_name],mean),b_egger_sd=sapply(all_snps[[1]][b_egger_name],sd),b_egger_se= sapply(all_snps[[1]][se_b_egger_name],mean) ,
                                   b_uv_egger=sapply(all_snps[[1]][b_uv_egger_name],mean),b_uv_egger_sd=sapply(all_snps[[1]][b_uv_egger_name],sd),b_uv_egger_se= sapply(all_snps[[1]][se_b_uv_egger_name],mean) ,
                                   b_DD=sapply(all_snps[[1]][b_DD_name],mean),b_DD_sd=sapply(all_snps[[1]][b_DD_name],sd),b_DD_se= sapply(all_snps[[1]][se_b_DD_name],mean) ,
                                   b_SH=sapply(all_snps[[1]][b_SH_name],mean),b_SH_sd=sapply(all_snps[[1]][b_SH_name],sd),b_SH_se= sapply(all_snps[[1]][se_b_SH_name],mean) ,
                                   b_true=sapply(all_snps[[1]][b_true_name],mean),b_true_sd=sapply(all_snps[[1]][b_true_name],sd)*0,b_true_se= sapply(all_snps[[1]][se_b_true_name],mean)*0 )
        slope_estimated=round(slope_estimated,digits=3)
        
      }
      
      
      
      
     
      point_est=list(H_only )
      
      write.table(slope_estimated,file = filename1,sep="\t",col.names = TRUE,quote = FALSE,row.names = FALSE)
      save(point_est,file=filename2)
      
    }
  }
  print(d)
}

a=round(H_only[sample(c(1:50),4,replace = FALSE),],2)
 
b1=paste("&$",a$beta_before[1],"$&$",a$beta_cml[1],"$&$",a$beta_uv_cml[1],"$&$",a$beta_egger[1],"$&$",a$beta_uv_egger[1],"$&$",a$beta_uv_DD[1],"$&$",a$beta_uv_SH[1],"$&$",
         a$error_before[1],"$&$",a$error_cml[1],"$&$",a$error_uv_cml[1],"$&$",a$error_egger[1],"$&$",a$error_uv_egger[1],"$&$",a$error_uv_DD[1],"$&$",a$error_uv_SH[1],"$\\")
       


b2=paste("&$(",a$sd_before[1],")$&$(",a$sd_cml[1],")$&$(",a$sd_uv_cml[1],")$&$(",a$sd_egger[1],
         ")$&$(",a$sd_uv_egger[1],")$&$(",a$sd_uv_DD[1],")$&$(",a$sd_uv_SH[1],")$&$(",
         a$se_error_before[1],")$&$(",a$se_error_cml[1],")$&$(",a$se_error_uv_cml[1],")$&$(",a$se_error_egger[1],
         ")$&$(",a$se_error_uv_egger[1],")$&$(",a$se_error_uv_DD[1],")$&$(",a$se_error_uv_SH[1],")$\\")

b3=paste("&$(",a$se_before[1],")$&$(",a$se_cml[1],")$&$(",a$se_uv_cml[1],")$&$(",a$se_egger[1],
         ")$&$(",a$se_uv_egger[1],")$&$(",a$se_uv_DD[1],")$&$(",a$se_uv_SH[1],")$&$",
         "NA","$&$", "NA","$&$", "NA","$&$", "NA",
         "$&$", "NA", "$&$", "NA" , "$&$" ,  "NA", "$\\")
 


library(vioplot)
load("D:/art sim/New/no_pleiotropy/uvmr_mvmr/result/d=2betaXYnot0rho0point_effect_estimates.RData")
a=point_est[[1]]
no= a$beta_before
cml=a$beta_cml
uv_cml=a$beta_uv_cml
egger=a$beta_egger
uv_egger=a$beta_uv_egger
DD=a$beta_uv_DD
SH=a$beta_uv_SH
 


vioplot(no,cml,egger,uv_cml,uv_egger,DD,SH ,names=c("No correction","MVMR-cML","MVMR-Egger","UVMR-cML","UVMR-Egger","DHO","SH" ),xlab="bias-correction methods",
        ylab="point estimates",horizontal = FALSE,
        main="Mean effect estimates of SNPs affecting covariates only, before and after bias correction",col=5)
abline(a=0,b=0 ,col="#3333FF",lwd=1,lty=2)


 