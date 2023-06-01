

Ncov=10
library(data.table)
library(tidyr)
library(dplyr)
library(R.utils)
setwd("D:/BMI")
if(Ncov>0){
  qqplot_name_after=numeric()
  for(i in 1:Ncov){
    qqplot_name_after[i]=paste("D:/BMI/qqplotBMIM1firstturn",i,".jpeg",sep="")
  }
  filename_after=numeric()
  for(i in 1:Ncov){
    filename_after[i]=paste("resultMAInfo",i,".txt",sep="")
  }
  result1=read.table(file="response.BMI.glm.linear",header = FALSE)
  result_met=data.frame(CHR=result1$V1,POS=result1$V2,SNP=result1$V3,beta=result1$V9,se=result1$V10,p=result1$V12)
  A=list()
  
  for(j in 1:Ncov){
    result2=read.table(file = filename_after[j],header = TRUE)
    result2$BP=result_met$POS
    result_adj=data.frame(CHR=result2$chr,POS=result2$BP,SNP=result2$rs,beta=result2$beta_adj,se=result2$sd_adj,p=result2$p_adj)
    A[[j]]=result_adj
    print(j)
  }
  lambda=c(round(median(qchisq(p=result_met$p,df=1,lower.tail = FALSE))/qchisq(p=0.5,df=1,lower.tail = FALSE),2))
  lambda_after=numeric()
  for(i in 1:Ncov){
    lambda_after[i]=median(qchisq(p=A[[i]]$p,df=1,lower.tail = FALSE))/qchisq(p=0.5,df=1,lower.tail = FALSE)
  }
  lambda_after=round(lambda_after,2)
  
  main=numeric()
  for(i in 1:10){
    if(i>1){
      main[i]=paste("QQ plot of BMI. ",i," PCs are adjusted" ,sep="")
    }
    if(i==1){
      main[i]=paste("QQ plot of BMI. ",i," PC is adjusted" ,sep="")
    }
  }
  
  sub=numeric()
  for(i in 1:Ncov){
    sub[i]=paste("lambda before adjustment is ", lambda,", lambda after adjustment is ", lambda_after[i])
  }
  
  
  
  y1=matrix(0,nrow=Ncov,ncol=nrow(A[[1]]))
  y2=matrix(0,nrow=Ncov,ncol=nrow(A[[1]]))
  y1a=matrix(0,nrow=Ncov,ncol=nrow(A[[1]]))
  y2a=matrix(0,nrow=Ncov,ncol=nrow(A[[1]]))
  addpy1=list()
  addpy2=list()
  addpx1=list()
  addpx2=list()
  x=sort(-log10(seq(from=1/nrow(A[[1]]),to=1,by=1/nrow(A[[1]]))))
  for(i in 1:Ncov){
    y1[i,]=sort(-log10(result_met$p))
    y1a[i,]=y1[i,]
    y1a[i,][y1a[i,]==Inf]=max(y1a[i,][y1a[i,]!=Inf])
    y2[i,]=sort(-log10(A[[i]]$p))
    y2a[i,]=y2[i,]
    y2a[i,][y2a[i,]==Inf]=max(y2a[i,][y2a[i,]!=Inf])
    addpy1[[i]]=rep(max(y1[i,][y1[i,]!=Inf]),length(y1[i,][y1[i,]==Inf]))
    addpy2[[i]]=rep(max(y1[i,][y1[i,]!=Inf]),length(y2[i,][y2[i,]==Inf]))
    addpx1[[i]]=x[y1[i,]==Inf]
    addpx2[[i]]=x[y2[i,]==Inf]
    }
  
 
  for(i in 1:Ncov){
    ystart=min(y1a[i,],y2a[i,])
    ystop=max(y1a[i,],y2a[i,])
    xstart=min(x)
    xstop=max(x)
    jpeg(filename =qqplot_name_after[i])
    plot(x,y1[i,],type = "n",main=main[i],sub=sub[i],xlab="Expected -log10(p-value)",ylab="Observed -log10(p-value)",ylim=c(ystart,ystop),xlim=c(xstart,xstop))
    abline(a=0,b=1,col=2)
    points(x,y1[i,],pch = 16,col="#CCCCCC",cex=2)
    points(x,y2[i,],pch = 17,col="#9933FF",cex=1)
    points(addpx1[[i]],addpy1[[i]],pch=1,col="#CCCCCC",cex=2)
    points(addpx2[[i]],addpy2[[i]],pch=24,col="#9933FF",cex=1)
    legend("topleft",legend = c("before","after","0 p-value before","0 p-value after"),pch=c(16,17,1,24),cex=c(1,1,1,1),col=c("#CCCCCC","#9933FF","#CCCCCC","#9933FF"))
    dev.off()
    }
}

setwd("D:/HC")
if(Ncov>0){
  qqplot_name_after=numeric()
  for(i in 1:Ncov){
    qqplot_name_after[i]=paste("D:/HC/qqplotHCM1firstturn",i,".jpeg",sep="")
  }
  filename_after=numeric()
  for(i in 1:Ncov){
    filename_after[i]=paste("resultMAInfo",i,".txt",sep="")
  }
  result1=read.table(file="response.HC.glm.linear",header = FALSE)
  result_met=data.frame(CHR=result1$V1,POS=result1$V2,SNP=result1$V3,beta=result1$V9,se=result1$V10,p=result1$V12)
  A=list()
  
  for(j in 1:Ncov){
    result2=read.table(file = filename_after[j],header = TRUE)
    result2$BP=result_met$POS
    result_adj=data.frame(CHR=result2$chr,POS=result2$BP,SNP=result2$rs,beta=result2$beta_adj,se=result2$sd_adj,p=result2$p_adj)
    A[[j]]=result_adj
    print(j)
  }
  lambda=c(round(median(qchisq(p=result_met$p,df=1,lower.tail = FALSE))/qchisq(p=0.5,df=1,lower.tail = FALSE),2))
  lambda_after=numeric()
  for(i in 1:Ncov){
    lambda_after[i]=median(qchisq(p=A[[i]]$p,df=1,lower.tail = FALSE))/qchisq(p=0.5,df=1,lower.tail = FALSE)
  }
  lambda_after=round(lambda_after,2)
  
  main=numeric()
  for(i in 1:10){
    if(i>1){
      main[i]=paste("QQ plot of HC. ",i," PCs are adjusted" ,sep="")
    }
    if(i==1){
      main[i]=paste("QQ plot of HC. ",i," PC is adjusted" ,sep="")
    }
  }
  
  sub=numeric()
  for(i in 1:Ncov){
    sub[i]=paste("lambda before adjustment is ", lambda,", lambda after adjustment is ", lambda_after[i])
  }
  
  
  
  y1=matrix(0,nrow=Ncov,ncol=nrow(A[[1]]))
  y2=matrix(0,nrow=Ncov,ncol=nrow(A[[1]]))
  y1a=matrix(0,nrow=Ncov,ncol=nrow(A[[1]]))
  y2a=matrix(0,nrow=Ncov,ncol=nrow(A[[1]]))
  addpy1=list()
  addpy2=list()
  addpx1=list()
  addpx2=list()
  x=sort(-log10(seq(from=1/nrow(A[[1]]),to=1,by=1/nrow(A[[1]]))))
  for(i in 1:Ncov){
    y1[i,]=sort(-log10(result_met$p))
    y1a[i,]=y1[i,]
    y1a[i,][y1a[i,]==Inf]=max(y1a[i,][y1a[i,]!=Inf])
    y2[i,]=sort(-log10(A[[i]]$p))
    y2a[i,]=y2[i,]
    y2a[i,][y2a[i,]==Inf]=max(y2a[i,][y2a[i,]!=Inf])
    addpy1[[i]]=rep(max(y1[i,][y1[i,]!=Inf]),length(y1[i,][y1[i,]==Inf]))
    addpy2[[i]]=rep(max(y1[i,][y1[i,]!=Inf]),length(y2[i,][y2[i,]==Inf]))
    addpx1[[i]]=x[y1[i,]==Inf]
    addpx2[[i]]=x[y2[i,]==Inf]
  }
  
  
  for(i in 1:Ncov){
    ystart=min(y1a[i,],y2a[i,])
    ystop=max(y1a[i,],y2a[i,])
    xstart=min(x)
    xstop=max(x)
    jpeg(filename =qqplot_name_after[i])
    plot(x,y1[i,],type = "n",main=main[i],sub=sub[i],xlab="Expected -log10(p-value)",ylab="Observed -log10(p-value)",ylim=c(ystart,ystop),xlim=c(xstart,xstop))
    abline(a=0,b=1,col=2)
    points(x,y1[i,],pch = 16,col="#CCCCCC",cex=2)
    points(x,y2[i,],pch = 17,col="#9933FF",cex=1)
    points(addpx1[[i]],addpy1[[i]],pch=1,col="#CCCCCC",cex=2)
    points(addpx2[[i]],addpy2[[i]],pch=24,col="#9933FF",cex=1)
    legend("topleft",legend = c("before","after","0 p-value before","0 p-value after"),pch=c(16,17,1,24),cex=c(1,1,1,1),col=c("#CCCCCC","#9933FF","#CCCCCC","#9933FF"))
    dev.off()
  }
}

setwd("D:/WC")
if(Ncov>0){
  qqplot_name_after=numeric()
  for(i in 1:Ncov){
    qqplot_name_after[i]=paste("D:/WC/qqplotWCM1firstturn",i,".jpeg",sep="")
  }
  filename_after=numeric()
  for(i in 1:Ncov){
    filename_after[i]=paste("resultMAInfo",i,".txt",sep="")
  }
  result1=read.table(file="response.WC.glm.linear",header = FALSE)
  result_met=data.frame(CHR=result1$V1,POS=result1$V2,SNP=result1$V3,beta=result1$V9,se=result1$V10,p=result1$V12)
  A=list()
  
  for(j in 1:Ncov){
    result2=read.table(file = filename_after[j],header = TRUE)
    result2$BP=result_met$POS
    result_adj=data.frame(CHR=result2$chr,POS=result2$BP,SNP=result2$rs,beta=result2$beta_adj,se=result2$sd_adj,p=result2$p_adj)
    A[[j]]=result_adj
    print(j)
  }
  lambda=c(round(median(qchisq(p=result_met$p,df=1,lower.tail = FALSE))/qchisq(p=0.5,df=1,lower.tail = FALSE),2))
  lambda_after=numeric()
  for(i in 1:Ncov){
    lambda_after[i]=median(qchisq(p=A[[i]]$p,df=1,lower.tail = FALSE))/qchisq(p=0.5,df=1,lower.tail = FALSE)
  }
  lambda_after=round(lambda_after,2)
  
  main=numeric()
  for(i in 1:10){
    if(i>1){
      main[i]=paste("QQ plot of WC. ",i," PCs are adjusted" ,sep="")
    }
    if(i==1){
      main[i]=paste("QQ plot of WC. ",i," PC is adjusted" ,sep="")
    }
  }
  
  sub=numeric()
  for(i in 1:Ncov){
    sub[i]=paste("lambda before adjustment is ", lambda,", lambda after adjustment is ", lambda_after[i])
  }
  
  
  
  y1=matrix(0,nrow=Ncov,ncol=nrow(A[[1]]))
  y2=matrix(0,nrow=Ncov,ncol=nrow(A[[1]]))
  y1a=matrix(0,nrow=Ncov,ncol=nrow(A[[1]]))
  y2a=matrix(0,nrow=Ncov,ncol=nrow(A[[1]]))
  addpy1=list()
  addpy2=list()
  addpx1=list()
  addpx2=list()
  x=sort(-log10(seq(from=1/nrow(A[[1]]),to=1,by=1/nrow(A[[1]]))))
  for(i in 1:Ncov){
    y1[i,]=sort(-log10(result_met$p))
    y1a[i,]=y1[i,]
    y1a[i,][y1a[i,]==Inf]=max(y1a[i,][y1a[i,]!=Inf])
    y2[i,]=sort(-log10(A[[i]]$p))
    y2a[i,]=y2[i,]
    y2a[i,][y2a[i,]==Inf]=max(y2a[i,][y2a[i,]!=Inf])
    addpy1[[i]]=rep(max(y1[i,][y1[i,]!=Inf]),length(y1[i,][y1[i,]==Inf]))
    addpy2[[i]]=rep(max(y1[i,][y1[i,]!=Inf]),length(y2[i,][y2[i,]==Inf]))
    addpx1[[i]]=x[y1[i,]==Inf]
    addpx2[[i]]=x[y2[i,]==Inf]
  }
  
  
  for(i in 1:Ncov){
    ystart=min(y1a[i,],y2a[i,])
    ystop=max(y1a[i,],y2a[i,])
    xstart=min(x)
    xstop=max(x)
    jpeg(filename =qqplot_name_after[i])
    plot(x,y1[i,],type = "n",main=main[i],sub=sub[i],xlab="Expected -log10(p-value)",ylab="Observed -log10(p-value)",ylim=c(ystart,ystop),xlim=c(xstart,xstop))
    abline(a=0,b=1,col=2)
    points(x,y1[i,],pch = 16,col="#CCCCCC",cex=2)
    points(x,y2[i,],pch = 17,col="#9933FF",cex=1)
    points(addpx1[[i]],addpy1[[i]],pch=1,col="#CCCCCC",cex=2)
    points(addpx2[[i]],addpy2[[i]],pch=24,col="#9933FF",cex=1)
    legend("topleft",legend = c("before","after","0 p-value before","0 p-value after"),pch=c(16,17,1,24),cex=c(1,1,1,1),col=c("#CCCCCC","#9933FF","#CCCCCC","#9933FF"))
    dev.off()
  }
}


setwd("D:/Height")
if(Ncov>0){
  qqplot_name_after=numeric()
  for(i in 1:Ncov){
    qqplot_name_after[i]=paste("D:/Height/qqplotHeightM1firstturn",i,".jpeg",sep="")
  }
  filename_after=numeric()
  for(i in 1:Ncov){
    filename_after[i]=paste("resultMAInfo",i,".txt",sep="")
  }
  result1=read.table(file="response.Height.glm.linear",header = FALSE)
  result_met=data.frame(CHR=result1$V1,POS=result1$V2,SNP=result1$V3,beta=result1$V9,se=result1$V10,p=result1$V12)
  A=list()
  
  for(j in 1:Ncov){
    result2=read.table(file = filename_after[j],header = TRUE)
    result2$BP=result_met$POS
    result_adj=data.frame(CHR=result2$chr,POS=result2$BP,SNP=result2$rs,beta=result2$beta_adj,se=result2$sd_adj,p=result2$p_adj)
    A[[j]]=result_adj
    print(j)
  }
  lambda=c(round(median(qchisq(p=result_met$p,df=1,lower.tail = FALSE))/qchisq(p=0.5,df=1,lower.tail = FALSE),2))
  lambda_after=numeric()
  for(i in 1:Ncov){
    lambda_after[i]=median(qchisq(p=A[[i]]$p,df=1,lower.tail = FALSE))/qchisq(p=0.5,df=1,lower.tail = FALSE)
  }
  lambda_after=round(lambda_after,2)
  
  main=numeric()
  for(i in 1:10){
    if(i>1){
      main[i]=paste("QQ plot of Height. ",i," PCs are adjusted" ,sep="")
    }
    if(i==1){
      main[i]=paste("QQ plot of Height. ",i," PC is adjusted" ,sep="")
    }
  }
  
  sub=numeric()
  for(i in 1:Ncov){
    sub[i]=paste("lambda before adjustment is ", lambda,", lambda after adjustment is ", lambda_after[i])
  }
  
  
  
  y1=matrix(0,nrow=Ncov,ncol=nrow(A[[1]]))
  y2=matrix(0,nrow=Ncov,ncol=nrow(A[[1]]))
  y1a=matrix(0,nrow=Ncov,ncol=nrow(A[[1]]))
  y2a=matrix(0,nrow=Ncov,ncol=nrow(A[[1]]))
  addpy1=list()
  addpy2=list()
  addpx1=list()
  addpx2=list()
  x=sort(-log10(seq(from=1/nrow(A[[1]]),to=1,by=1/nrow(A[[1]]))))
  for(i in 1:Ncov){
    y1[i,]=sort(-log10(result_met$p))
    y1a[i,]=y1[i,]
    y1a[i,][y1a[i,]==Inf]=max(y1a[i,][y1a[i,]!=Inf])
    y2[i,]=sort(-log10(A[[i]]$p))
    y2a[i,]=y2[i,]
    y2a[i,][y2a[i,]==Inf]=max(y2a[i,][y2a[i,]!=Inf])
    addpy1[[i]]=rep(max(y1[i,][y1[i,]!=Inf]),length(y1[i,][y1[i,]==Inf]))
    addpy2[[i]]=rep(max(y1[i,][y1[i,]!=Inf]),length(y2[i,][y2[i,]==Inf]))
    addpx1[[i]]=x[y1[i,]==Inf]
    addpx2[[i]]=x[y2[i,]==Inf]
  }
  
  
  for(i in 1:Ncov){
    ystart=min(y1a[i,],y2a[i,])
    ystop=max(y1a[i,],y2a[i,])
    xstart=min(x)
    xstop=max(x)
    jpeg(filename =qqplot_name_after[i])
    plot(x,y1[i,],type = "n",main=main[i],sub=sub[i],xlab="Expected -log10(p-value)",ylab="Observed -log10(p-value)",ylim=c(ystart,ystop),xlim=c(xstart,xstop))
    abline(a=0,b=1,col=2)
    points(x,y1[i,],pch = 16,col="#CCCCCC",cex=2)
    points(x,y2[i,],pch = 17,col="#9933FF",cex=1)
    points(addpx1[[i]],addpy1[[i]],pch=1,col="#CCCCCC",cex=2)
    points(addpx2[[i]],addpy2[[i]],pch=24,col="#9933FF",cex=1)
    legend("topleft",legend = c("before","after","0 p-value before","0 p-value after"),pch=c(16,17,1,24),cex=c(1,1,1,1),col=c("#CCCCCC","#9933FF","#CCCCCC","#9933FF"))
    dev.off()
  }
}


setwd("D:/WHR")
if(Ncov>0){
  qqplot_name_after=numeric()
  for(i in 1:Ncov){
    qqplot_name_after[i]=paste("D:/WHR/qqplotWHRM1firstturn",i,".jpeg",sep="")
  }
  filename_after=numeric()
  for(i in 1:Ncov){
    filename_after[i]=paste("resultMAInfo",i,".txt",sep="")
  }
  result1=read.table(file="response.WHR.glm.linear",header = FALSE)
  result_met=data.frame(CHR=result1$V1,POS=result1$V2,SNP=result1$V3,beta=result1$V9,se=result1$V10,p=result1$V12)
  A=list()
  
  for(j in 1:Ncov){
    result2=read.table(file = filename_after[j],header = TRUE)
    result2$BP=result_met$POS
    result_adj=data.frame(CHR=result2$chr,POS=result2$BP,SNP=result2$rs,beta=result2$beta_adj,se=result2$sd_adj,p=result2$p_adj)
    A[[j]]=result_adj
    print(j)
  }
  lambda=c(round(median(qchisq(p=result_met$p,df=1,lower.tail = FALSE))/qchisq(p=0.5,df=1,lower.tail = FALSE),2))
  lambda_after=numeric()
  for(i in 1:Ncov){
    lambda_after[i]=median(qchisq(p=A[[i]]$p,df=1,lower.tail = FALSE))/qchisq(p=0.5,df=1,lower.tail = FALSE)
  }
  lambda_after=round(lambda_after,2)
  
  main=numeric()
  for(i in 1:10){
    if(i>1){
      main[i]=paste("QQ plot of WHR. ",i," PCs are adjusted" ,sep="")
    }
    if(i==1){
      main[i]=paste("QQ plot of WHR. ",i," PC is adjusted" ,sep="")
    }
  }
  
  sub=numeric()
  for(i in 1:Ncov){
    sub[i]=paste("lambda before adjustment is ", lambda,", lambda after adjustment is ", lambda_after[i])
  }
  
  
  
  y1=matrix(0,nrow=Ncov,ncol=nrow(A[[1]]))
  y2=matrix(0,nrow=Ncov,ncol=nrow(A[[1]]))
  y1a=matrix(0,nrow=Ncov,ncol=nrow(A[[1]]))
  y2a=matrix(0,nrow=Ncov,ncol=nrow(A[[1]]))
  addpy1=list()
  addpy2=list()
  addpx1=list()
  addpx2=list()
  x=sort(-log10(seq(from=1/nrow(A[[1]]),to=1,by=1/nrow(A[[1]]))))
  for(i in 1:Ncov){
    y1[i,]=sort(-log10(result_met$p))
    y1a[i,]=y1[i,]
    y1a[i,][y1a[i,]==Inf]=max(y1a[i,][y1a[i,]!=Inf])
    y2[i,]=sort(-log10(A[[i]]$p))
    y2a[i,]=y2[i,]
    y2a[i,][y2a[i,]==Inf]=max(y2a[i,][y2a[i,]!=Inf])
    addpy1[[i]]=rep(max(y1[i,][y1[i,]!=Inf]),length(y1[i,][y1[i,]==Inf]))
    addpy2[[i]]=rep(max(y1[i,][y1[i,]!=Inf]),length(y2[i,][y2[i,]==Inf]))
    addpx1[[i]]=x[y1[i,]==Inf]
    addpx2[[i]]=x[y2[i,]==Inf]
  }
  
  
  for(i in 1:Ncov){
    ystart=min(y1a[i,],y2a[i,])
    ystop=max(y1a[i,],y2a[i,])
    xstart=min(x)
    xstop=max(x)
    jpeg(filename =qqplot_name_after[i])
    plot(x,y1[i,],type = "n",main=main[i],sub=sub[i],xlab="Expected -log10(p-value)",ylab="Observed -log10(p-value)",ylim=c(ystart,ystop),xlim=c(xstart,xstop))
    abline(a=0,b=1,col=2)
    points(x,y1[i,],pch = 16,col="#CCCCCC",cex=2)
    points(x,y2[i,],pch = 17,col="#9933FF",cex=1)
    points(addpx1[[i]],addpy1[[i]],pch=1,col="#CCCCCC",cex=2)
    points(addpx2[[i]],addpy2[[i]],pch=24,col="#9933FF",cex=1)
    legend("topleft",legend = c("before","after","0 p-value before","0 p-value after"),pch=c(16,17,1,24),cex=c(1,1,1,1),col=c("#CCCCCC","#9933FF","#CCCCCC","#9933FF"))
    dev.off()
  }
}


setwd("D:/SBP")
if(Ncov>0){
  qqplot_name_after=numeric()
  for(i in 1:Ncov){
    qqplot_name_after[i]=paste("D:/SBP/qqplotSBPM1firstturn",i,".jpeg",sep="")
  }
  filename_after=numeric()
  for(i in 1:Ncov){
    filename_after[i]=paste("resultMAInfo",i,".txt",sep="")
  }
  result1=read.table(file="response.SBP.glm.linear",header = FALSE)
  result_met=data.frame(CHR=result1$V1,POS=result1$V2,SNP=result1$V3,beta=result1$V9,se=result1$V10,p=result1$V12)
  A=list()
  
  for(j in 1:Ncov){
    result2=read.table(file = filename_after[j],header = TRUE)
    result2$BP=result_met$POS
    result_adj=data.frame(CHR=result2$chr,POS=result2$BP,SNP=result2$rs,beta=result2$beta_adj,se=result2$sd_adj,p=result2$p_adj)
    A[[j]]=result_adj
    print(j)
  }
  lambda=c(round(median(qchisq(p=result_met$p,df=1,lower.tail = FALSE))/qchisq(p=0.5,df=1,lower.tail = FALSE),2))
  lambda_after=numeric()
  for(i in 1:Ncov){
    lambda_after[i]=median(qchisq(p=A[[i]]$p,df=1,lower.tail = FALSE))/qchisq(p=0.5,df=1,lower.tail = FALSE)
  }
  lambda_after=round(lambda_after,2)
  
  main=numeric()
  for(i in 1:10){
    if(i>1){
      main[i]=paste("QQ plot of SBP. ",i," PCs are adjusted" ,sep="")
    }
    if(i==1){
      main[i]=paste("QQ plot of SBP. ",i," PC is adjusted" ,sep="")
    }
  }
  
  sub=numeric()
  for(i in 1:Ncov){
    sub[i]=paste("lambda before adjustment is ", lambda,", lambda after adjustment is ", lambda_after[i])
  }
  
  
  
  y1=matrix(0,nrow=Ncov,ncol=nrow(A[[1]]))
  y2=matrix(0,nrow=Ncov,ncol=nrow(A[[1]]))
  y1a=matrix(0,nrow=Ncov,ncol=nrow(A[[1]]))
  y2a=matrix(0,nrow=Ncov,ncol=nrow(A[[1]]))
  addpy1=list()
  addpy2=list()
  addpx1=list()
  addpx2=list()
  x=sort(-log10(seq(from=1/nrow(A[[1]]),to=1,by=1/nrow(A[[1]]))))
  for(i in 1:Ncov){
    y1[i,]=sort(-log10(result_met$p))
    y1a[i,]=y1[i,]
    y1a[i,][y1a[i,]==Inf]=max(y1a[i,][y1a[i,]!=Inf])
    y2[i,]=sort(-log10(A[[i]]$p))
    y2a[i,]=y2[i,]
    y2a[i,][y2a[i,]==Inf]=max(y2a[i,][y2a[i,]!=Inf])
    addpy1[[i]]=rep(max(y1[i,][y1[i,]!=Inf]),length(y1[i,][y1[i,]==Inf]))
    addpy2[[i]]=rep(max(y1[i,][y1[i,]!=Inf]),length(y2[i,][y2[i,]==Inf]))
    addpx1[[i]]=x[y1[i,]==Inf]
    addpx2[[i]]=x[y2[i,]==Inf]
  }
  
  
  for(i in 1:Ncov){
    ystart=min(y1a[i,],y2a[i,])
    ystop=max(y1a[i,],y2a[i,])
    xstart=min(x)
    xstop=max(x)
    jpeg(filename =qqplot_name_after[i])
    plot(x,y1[i,],type = "n",main=main[i],sub=sub[i],xlab="Expected -log10(p-value)",ylab="Observed -log10(p-value)",ylim=c(ystart,ystop),xlim=c(xstart,xstop))
    abline(a=0,b=1,col=2)
    points(x,y1[i,],pch = 16,col="#CCCCCC",cex=2)
    points(x,y2[i,],pch = 17,col="#9933FF",cex=1)
    points(addpx1[[i]],addpy1[[i]],pch=1,col="#CCCCCC",cex=2)
    points(addpx2[[i]],addpy2[[i]],pch=24,col="#9933FF",cex=1)
    legend("topleft",legend = c("before","after","0 p-value before","0 p-value after"),pch=c(16,17,1,24),cex=c(1,1,1,1),col=c("#CCCCCC","#9933FF","#CCCCCC","#9933FF"))
    dev.off()
  }
}

setwd("D:/DBP")
if(Ncov>0){
  qqplot_name_after=numeric()
  for(i in 1:Ncov){
    qqplot_name_after[i]=paste("D:/DBP/qqplotDBPM1firstturn",i,".jpeg",sep="")
  }
  filename_after=numeric()
  for(i in 1:Ncov){
    filename_after[i]=paste("resultMAInfo",i,".txt",sep="")
  }
  result1=read.table(file="response.DBP.glm.linear",header = FALSE)
  result_met=data.frame(CHR=result1$V1,POS=result1$V2,SNP=result1$V3,beta=result1$V9,se=result1$V10,p=result1$V12)
  A=list()
  
  for(j in 1:Ncov){
    result2=read.table(file = filename_after[j],header = TRUE)
    result2$BP=result_met$POS
    result_adj=data.frame(CHR=result2$chr,POS=result2$BP,SNP=result2$rs,beta=result2$beta_adj,se=result2$sd_adj,p=result2$p_adj)
    A[[j]]=result_adj
    print(j)
  }
  lambda=c(round(median(qchisq(p=result_met$p,df=1,lower.tail = FALSE))/qchisq(p=0.5,df=1,lower.tail = FALSE),2))
  lambda_after=numeric()
  for(i in 1:Ncov){
    lambda_after[i]=median(qchisq(p=A[[i]]$p,df=1,lower.tail = FALSE))/qchisq(p=0.5,df=1,lower.tail = FALSE)
  }
  lambda_after=round(lambda_after,2)
  
  main=numeric()
  for(i in 1:10){
    if(i>1){
      main[i]=paste("QQ plot of DBP. ",i," PCs are adjusted" ,sep="")
    }
    if(i==1){
      main[i]=paste("QQ plot of DBP. ",i," PC is adjusted" ,sep="")
    }
  }
  
  sub=numeric()
  for(i in 1:Ncov){
    sub[i]=paste("lambda before adjustment is ", lambda,", lambda after adjustment is ", lambda_after[i])
  }
  
  
  
  y1=matrix(0,nrow=Ncov,ncol=nrow(A[[1]]))
  y2=matrix(0,nrow=Ncov,ncol=nrow(A[[1]]))
  y1a=matrix(0,nrow=Ncov,ncol=nrow(A[[1]]))
  y2a=matrix(0,nrow=Ncov,ncol=nrow(A[[1]]))
  addpy1=list()
  addpy2=list()
  addpx1=list()
  addpx2=list()
  x=sort(-log10(seq(from=1/nrow(A[[1]]),to=1,by=1/nrow(A[[1]]))))
  for(i in 1:Ncov){
    y1[i,]=sort(-log10(result_met$p))
    y1a[i,]=y1[i,]
    y1a[i,][y1a[i,]==Inf]=max(y1a[i,][y1a[i,]!=Inf])
    y2[i,]=sort(-log10(A[[i]]$p))
    y2a[i,]=y2[i,]
    y2a[i,][y2a[i,]==Inf]=max(y2a[i,][y2a[i,]!=Inf])
    addpy1[[i]]=rep(max(y1[i,][y1[i,]!=Inf]),length(y1[i,][y1[i,]==Inf]))
    addpy2[[i]]=rep(max(y1[i,][y1[i,]!=Inf]),length(y2[i,][y2[i,]==Inf]))
    addpx1[[i]]=x[y1[i,]==Inf]
    addpx2[[i]]=x[y2[i,]==Inf]
  }
  
  
  for(i in 1:Ncov){
    ystart=min(y1a[i,],y2a[i,])
    ystop=max(y1a[i,],y2a[i,])
    xstart=min(x)
    xstop=max(x)
    jpeg(filename =qqplot_name_after[i])
    plot(x,y1[i,],type = "n",main=main[i],sub=sub[i],xlab="Expected -log10(p-value)",ylab="Observed -log10(p-value)",ylim=c(ystart,ystop),xlim=c(xstart,xstop))
    abline(a=0,b=1,col=2)
    points(x,y1[i,],pch = 16,col="#CCCCCC",cex=2)
    points(x,y2[i,],pch = 17,col="#9933FF",cex=1)
    points(addpx1[[i]],addpy1[[i]],pch=1,col="#CCCCCC",cex=2)
    points(addpx2[[i]],addpy2[[i]],pch=24,col="#9933FF",cex=1)
    legend("topleft",legend = c("before","after","0 p-value before","0 p-value after"),pch=c(16,17,1,24),cex=c(1,1,1,1),col=c("#CCCCCC","#9933FF","#CCCCCC","#9933FF"))
    dev.off()
  }
}


