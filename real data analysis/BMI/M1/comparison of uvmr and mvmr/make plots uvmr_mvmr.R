
 

Ncov=2
if(Ncov>0){
  setwd("D:/art sim/New/newGWAS/BMIMETPC/M1/second_turn/uvmr_mvmr/")
  library(hudson)
  
  method=c("UVMR-Median","UVMR-cML","UVMR-Lasso","UVMR-IVW","UVMR-Egger","DHO","SH")
  subtitle=method
  Nmethod=length(method)
  filename=c("resultUVMRMedian2.txt","resultUVMRcML2.txt","resultUVMRLasso2.txt","resultUVMRIVW2.txt","resultUVMREgger2.txt","resultDD2.txt","resultSH2.txt")
  effect_plot_title=numeric()
  variance_plot_title=numeric()
  manhattan_plot_title_before=numeric()
  manhattan_plot_title_after=numeric()
  man_plot_name=numeric()
  effect_plot_name=numeric()
  variance_plot_name=numeric()
  outplace="/art sim/New/newGWAS/BMIMETPC/M1/second_turn/uvmr_mvmr/"
  for(i in 1:Nmethod){
    man_plot_name[i]=paste(outplace,"manhattan/",method[i],"mahattan",Ncov,sep="")
    effect_plot_name[i]=paste(outplace,"effect/",method[i],"effect",Ncov,".jpeg",sep="")
    variance_plot_name[i]=paste(outplace,"effect/",method[i],"variance",Ncov,".jpeg",sep="")
  }
  for(i in 1:Nmethod){
    effect_plot_title[i]=paste("SNP effect estimates of BMI, ", method[i], sep="")
    variance_plot_title[i]=paste("Standard error of SNP effect estimates, ", method[i], sep="")
    manhattan_plot_title_before[i]=paste("Manhattan plot of BMI before correction",sep="")
    manhattan_plot_title_after[i]=paste("Manhattan plot of BMI after ", method[i]," correction",sep="")
  }
  
  result1=read.table(file="response_2.BMI.glm.linear",header = FALSE)
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
  for(j in 1:Nmethod){
    result2=read.table(file = filename[j],header = TRUE)
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
  
  
  for(i in 1:Nmethod){
    gmirror(top=A[[i]],bottom = result_met,tline=5e-8,bline=5e-8,toptitle =manhattan_plot_title_after[i],bottomtitle = manhattan_plot_title_before[i],file = man_plot_name[i],res=100)
  }
  
  for(i in 1:Nmethod){
    jpeg(file = variance_plot_name[i])
    result_adj=A[[i]]
    a=data.frame(se1=result_met$se,se2=result_adj$se)
    atotal=a[(result_met$p<=5e-08)|(result_adj$p<=5e-08),]
    a1=a[(result_met$p<=5e-08)&(result_adj$p>5e-08),]
    a2=a[(result_met$p>5e-08)&(result_adj$p<=5e-08),]
    a3=a[(result_met$p<=5e-08)&(result_adj$p<=5e-08),]
    
    plot(atotal[,1],atotal[,2],main=variance_plot_title[i], xlab="before correction ",ylab="after correction ",type ="n")
    abline(a=0,b=1,lwd=1)
    points(a1[,1],a1[,2],pch=16,col="#666666",cex=1.5)
    points(a3[,1],a3[,2],pch=15,col="#CC9900",cex=1.5)
    points(a2[,1],a2[,2],pch=17,col="#9933FF",cex=1.5)
    
    legend("bottomright",pch=c(16,17,15),col=c(col="#666666","#9933FF","#CC9900"),cex=1.5,legend = c("before","after","both"))
    dev.off()
  }
  for(i in 1:Nmethod){
    jpeg(file = effect_plot_name[i])
    result_adj=A[[i]]
    b=data.frame(beta1=result_met$beta,beta2=result_adj$beta)
    btotal=b[(result_met$p<=5e-08)|(result_adj$p<=5e-08),]
    b1=b[(result_met$p<=5e-08)&(result_adj$p>5e-08),]
    b2=b[(result_met$p>5e-08)&(result_adj$p<=5e-08),]
    b3=b[(result_met$p<=5e-08)&(result_adj$p<=5e-08),]
    a=data.frame(se1=result_met$se,se2=result_adj$se)
    atotal=a[(result_met$p<=5e-08)|(result_adj$p<=5e-08),]
    a1=a[(result_met$p<=5e-08)&(result_adj$p>5e-08),]
    a2=a[(result_met$p>5e-08)&(result_adj$p<=5e-08),]
    a3=a[(result_met$p<=5e-08)&(result_adj$p<=5e-08),]
    y11=numeric()
    y12=numeric()
    y13=numeric()
    y21=numeric()
    y22=numeric()
    y23=numeric()
    x11=numeric()
    x12=numeric()
    x13=numeric()
    x21=numeric()
    x22=numeric()
    x23=numeric()
    for(j in 1:nrow(a1)){
      y11[j]=(b1[j,][2]-a1[j,][2])[1,1]
      y21[j]=(b1[j,][2]+a1[j,][2])[1,1]
      x11[j]=(b1[j,][1]-a1[j,][1])[1,1]
      x21[j]=(b1[j,][1]+a1[j,][1])[1,1]
    }
    for(j in 1:nrow(a2)){
      y12[j]=(b2[j,][2]-a2[j,][2])[1,1]
      y22[j]=(b2[j,][2]+a2[j,][2])[1,1]
      x12[j]=(b2[j,][1]-a2[j,][1])[1,1]
      x22[j]=(b2[j,][1]+a2[j,][1])[1,1]
    }
    for(j in 1:nrow(a3)){
      y13[j]=(b3[j,][2]-a3[j,][2])[1,1]
      y23[j]=(b3[j,][2]+a3[j,][2])[1,1]
      x13[j]=(b3[j,][1]-a3[j,][1])[1,1]
      x23[j]=(b3[j,][1]+a3[j,][1])[1,1]
    }
    y1=c(y11,y12,y13)
    y2=c(y21,y22,y23)
    x1=c(x11,x12,x13)
    x2=c(x21,x22,x23)
    y1=na.omit(y1)
    y2=na.omit(y2)
    x1=na.omit(x1)
    x2=na.omit(x2)
    plot(btotal[,1],btotal[,2],main=effect_plot_title[i],xlab="before correction",ylab="after correction",type ="n",xlim=c(min(x1),max(x2)),ylim=c(min(y1),max(y2)))
    abline(a=0,b=1,lwd=1)
    
    for(j in 1:nrow(a1)){
      segments(x0=(b1[j,][1]-a1[j,][1])[1,1],y0=(b1[j,][2])[1,1],x1=(b1[j,][1]+a1[j,][1])[1,1],y1=(b1[j,][2])[1,1],lwd=0.2,col="#999999")
      segments(x0=(b1[j,][1])[1,1],y0=(b1[j,][2]-a1[j,][2])[1,1],x1=(b1[j,][1])[1,1],y1=(b1[j,][2]+a1[j,][2])[1,1],lwd=0.2,col="#999999")
    }
    for(j in 1:nrow(a3)){
      segments(x0=(b3[j,][1]-a2[j,][1])[1,1],y0=(b3[j,][2])[1,1],x1=(b3[j,][1]+a2[j,][1])[1,1],y1=(b3[j,][2])[1,1],lwd=0.2,col="#CC9900")
      segments(x0=(b3[j,][1])[1,1],y0=(b3[j,][2]-a2[j,][2])[1,1],x1=(b3[j,][1])[1,1],y1=(b3[j,][2]+a2[j,][2])[1,1],lwd=0.2,col="#CC9900")
    }
    for(j in 1:nrow(a2)){
      segments(x0=(b2[j,][1]-a2[j,][1])[1,1],y0=(b2[j,][2])[1,1],x1=(b2[j,][1]+a2[j,][1])[1,1],y1=(b2[j,][2])[1,1],lwd=0.2,col="#9933FF")
      segments(x0=(b2[j,][1])[1,1],y0=(b2[j,][2]-a2[j,][2])[1,1],x1=(b2[j,][1])[1,1],y1=(b2[j,][2]+a2[j,][2])[1,1],lwd=0.2,col="#9933FF")
    }
    
    points(b1[,1],b1[,2],pch=16,col="#666666", cex=1.5)
    points(b3[,1],b3[,2],pch=15,col="#CC9900", cex=1.5)
    points(b2[,1],b2[,2],pch=17,col="#9933FF", cex=1.5)
    
    legend("bottomright",pch=c(16,17,15),col=c(col="#666666","#9933FF","#CC9900"),cex=1.5,legend = c("before","after","both"))
    dev.off()
  }
  
  
} 