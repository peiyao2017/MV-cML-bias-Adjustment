

setwd('D:/art sim/New/newGWAS/BMIMETPC/')
library(doParallel)

Ncov=5
library(mvtnorm)
library(doParallel)
library(MASS)
library(Matrix)
 
library(BEDMatrix)
G=BEDMatrix(path="D:/art sim/New/newGWAS/110K_QCed1.bed",simple_names = TRUE)
gene=read.table(file="D:/art sim/New/newGWAS/110K_QCed0.001.bim",header = FALSE)
gene2=read.table(file="D:/art sim/New/newGWAS/110K_QCed1.fam",header = FALSE)

IVID=intersect(gene$V2,colnames(G))


atotal=list()
btotal=list()
response_gwas_name=numeric()
for(i in 1:Ncov){
  response_gwas_name[i]=paste("response_",i,".BMI.glm.linear",sep="")
}
for(i in 1:Ncov){
  btotal[[i]]=read.table(file=response_gwas_name[i])
}
BP=btotal[[1]]$V2
Nsample=nrow(gene2)
nameX=numeric()
for(i in 1:Ncov){
  nameX[i]=paste("MET.PC",i,".glm.linear",sep="")
}

for(i in 1:Ncov){
  atotal[[i]]=read.table(file=nameX[i],header = FALSE)
}
gene1=atotal[[1]]
chr=atotal[[1]]$V1
betaGX1total=matrix(0,nrow=nrow(atotal[[1]]),ncol=length(atotal))
pvaluebetaGX1total=matrix(0,nrow=nrow(atotal[[1]]),ncol=length(atotal))
sdbetaGX1total=matrix(0,nrow=nrow(atotal[[1]]),ncol=length(atotal))
ZX1total=matrix(0,nrow=nrow(atotal[[1]]),ncol=length(atotal))

betaGYCtotal=matrix(0,nrow=nrow(btotal[[1]]),ncol = length(btotal))
pvaluebetaGYCtotal=matrix(0,nrow=nrow(btotal[[1]]),ncol = length(btotal))
sdbetaGYCtotal=matrix(0,nrow=nrow(btotal[[1]]),ncol = length(btotal))
ZYCtotal=matrix(0,nrow=nrow(btotal[[1]]),ncol = length(btotal))

for(i in 1:Ncov){
  betaGYCtotal[,i]=btotal[[i]]$V9
  pvaluebetaGYCtotal[,i]=btotal[[i]]$V12
  sdbetaGYCtotal[,i]=btotal[[i]]$V10
  ZYCtotal[,i]=btotal[[i]]$V11
}


ID=atotal[[1]]$V3
nsnpstotal=nrow(atotal[[1]])

for(i in 1:ncol(betaGX1total)){
  betaGX1total[,i]=atotal[[i ]]$V9
  pvaluebetaGX1total[,i]=atotal[[i ]]$V12
  sdbetaGX1total[,i]=atotal[[i ]]$V10
  ZX1total[,i]=atotal[[i ]]$V11
}


corr_total=list()
for(i in 1:Ncov){
  corr_total[[i]]=diag(1, nrow=i+1,ncol=i+1)
}


for(i in 1:length(btotal)){
  for(j in 1:(i+1) ){
    for(k in 1:(i+1)){
      if(j==1&k>1){
        x1=btotal[[i]]
        x2=atotal[[k-1]]
        y1=x1$V11[x1$V12>0.1&x2$V3%in%IVID&x2$V12>0.1&x2$V3%in%IVID]
        y2=x2$V11[x1$V12>0.1&x2$V3%in%IVID&x2$V12>0.1&x2$V3%in%IVID]
        corr_total[[i]][j,k]=cor(y1,y2)
      }
      if(k==1&j>1){
        x1=btotal[[i]]
        x2=atotal[[j-1]]
        y1=x1$V11[x1$V12>0.1&x2$V3%in%IVID&x2$V12>0.1&x2$V3%in%IVID]
        y2=x2$V11[x1$V12>0.1&x2$V3%in%IVID&x2$V12>0.1&x2$V3%in%IVID]
        corr_total[[i]][j,k]=cor(y1,y2)
      }
      if(k>1&j>1){
        x1=atotal[[k-1]]
        x2=atotal[[j-1]]
        y1=x1$V11[x1$V12>0.1&x2$V3%in%IVID&x2$V12>0.1&x2$V3%in%IVID]
        y2=x2$V11[x1$V12>0.1&x2$V3%in%IVID&x2$V12>0.1&x2$V3%in%IVID]
        corr_total[[i]][j,k]=cor(y1,y2)
      }
    }
  }
}

numcov=1:Ncov
name=numeric()
name1=numeric()
for(i in 1:length(numcov)){
  name[i]=paste("resultMAInfo_second_turn",numcov[i],".txt",sep="")
  name1[i]=paste("slope_and_covariate_effect_second_turn",numcov[i],".RData",sep="")
}






l=function(betaxtemp=betaGX_IV,betayctemp=betaGYC_IV,SIGtemp=SIG,betaytemp=betaGYC_IV,btemp=rep(0,times=ncol(betaGX_IV))){
  a1=numeric()
  a2=numeric()
  for(i in 1:length(SIGtemp)){
    a1[i]=t(c(betayctemp[i]-t(btemp)%*%betaxtemp[i,]-betaytemp[i],betaGX_IV[i,]-betaxtemp[i,]))%*%solve(SIGtemp[[i]])%*%c(betayctemp[i]-t(btemp)%*%betaxtemp[i,]-betaytemp[i],betaGX_IV[i,]-betaxtemp[i,])
    a2[i]=t(c(betayctemp[i]-t(btemp)%*%betaxtemp[i,]-0,betaGX_IV[i,]-betaxtemp[i,]))%*%solve(SIGtemp[[i]])%*%c(betayctemp[i]-t(btemp)%*%betaxtemp[i,]-0,betaGX_IV[i,]-betaxtemp[i,])
  }
  return(list(d=a2-a1,like=sum(a1)))
}

cml_MA=function(N=Nsample,p=ncol(betaGX1),m=nrow(betaGX_IV),betax=betaGX_IV,betayc=betaGYC_IV,SIG1=SIG,threshold=0.01,maxit=500){
  valid=list()
  Bic=numeric()
  btotal=matrix(0,ncol=length(0:(m-(p+1))),nrow=p)
  betaxesttotal=list()
  betayesttotal=matrix(0,ncol=length(0:(m-(p+1))),nrow=length(betayc))
  inv1=list()
  for(i in 1:m){
    inv1[[i]]=solve(SIG1[[i]])
  }
  for(k in 0:(m-(p+1))){
    b=rep(0,times=p)
    betayest=rep(0,times=m)
    betaxest=matrix(0,nrow=m,ncol=p)
    c=1
    repeat{
      c=c+1
      print(c)
      old=b
      for(i in 1:m){
        inv=inv1[[i]]
        A11=inv[1,1]
        A12=inv[1,(2:ncol(inv))]
        A22=inv[(2:ncol(inv)),(2:ncol(inv))]
        betayest[i]=betayc[i]-t(b)%*%betaxest[i,]+1/A11*t(betax[i,]-betaxest[i,])%*%A12
      } 
      total=data.frame(num=c(1:m),go=l(betaxtemp=betaxest,betayctemp=betayc,SIGtemp=SIG1,betaytemp=betayest,btemp=b)$d,betayest,betayc)
      total=total[order(-total$go),]
      for(i in 1:m){
        if(i>k){
          total$betayest[i]=0
        }
      }
      total=total[order(total$num),]
      betayest=total$betayest
      for(i in 1:m){
        inv=inv1[[i]]
        A11=inv[1,1]
        A12=inv[1,(2:ncol(inv))]
        A22=inv[(2:ncol(inv)),(2:ncol(inv))]
        betaxest[i,]=solve(A11*b%*%t(b)+A12%*%t(b)+b%*%t(A12)+A22)%*%((b%*%t(A12)+A22)%*%betax[i,]+(betayc[i]-betayest[i])*(A12+A11*b))
      }
      a=matrix(0,p,p)
      b=matrix(0,nrow=p,ncol=1)
      for(i in 1:m){
        inv=inv1[[i]]
        A11=inv[1,1]
        A12=inv[1,(2:ncol(inv))]
        a=a+betaxest[i,]%*%t(betaxest[i,])*A11
        b=b+t((A11*(betayc[i]-betayest[i])+t(A12)%*%(betax[i,]-betaxest[i,]))%*%betaxest[i,])
      }
      b=solve(a)%*%b
      new=b
      if(sum(abs(old-new))<threshold|c>maxit){
        break
      }
    }
    btotal[,k+1]=b
    betayesttotal[,k+1]=betayest
    betaxesttotal[[k+1]]=betaxest
    print(k)
    Bic[k+1]=l(betaxtemp=betaxest,betayctemp=betayc,SIGtemp=SIG1,betaytemp=betayest,btemp=b)$like+log(N)*k
    valid[[k+1]]=which(betayest==0)
    print(paste("current is",k))
  }
  b=btotal[,which.min(Bic)]
  validuse=valid[[which.min(Bic)]]
  betayest=betayesttotal[,which.min(Bic)]
  betaxest=betaxesttotal[[which.min(Bic)]]
  siginv=list()
  betaxvalid=matrix(0,nrow=length(validuse),ncol = ncol(betaxest))
  betaxIVvalid=matrix(0,nrow=length(validuse),ncol = ncol(betaxest))
  betaGYCvalid=numeric()
  for(i in 1:length(validuse)){
    siginv[[i]]=inv1[[validuse[i]]]
    betaxvalid[i,]=betaxest[validuse[i],]
    betaGYCvalid[i]=betayc[validuse[i]]
    betaxIVvalid[i,]=betax[validuse[i],]
  }
  B=list()
  for(i in 1:length(validuse)){
    a=matrix(0,nrow=ncol(betaxvalid),ncol=ncol(betaxvalid))
    a11=siginv[[i]][1,1]
    A12=as.matrix(siginv[[i]][1,(2:ncol(siginv[[i]]))])
    A22=as.matrix(siginv[[i]][(2:ncol(siginv[[i]])),(2:ncol(siginv[[i]]))])
    for(j in 1:nrow(a)){
      for(k in 1:nrow(a)){
        if(j==k){
          a[j,k]=(b[j]*betaxvalid[i,j]+t(b)%*%betaxvalid[i,]-betaGYCvalid[i])*a11-t(A12)%*%(betaxIVvalid[i,]-betaxvalid[i,])+A12[j]*betaxvalid[i,j]
        }
        if(j!=k){
          a[j,k]=b[j]*betaxvalid[i,k]*a11+A12[j]*betaxvalid[i,k]
        }
      }
    }
    B[[i]]=a
  }
  B1=B[[1]]
  for(i in 2:length(B)){
    B1=rbind(B1,B[[i]])
  }
  C=list()
  for(i in 1:length(validuse)){
    a=matrix(0,nrow=ncol(betaxvalid),ncol=ncol(betaxvalid))
    a11=siginv[[i]][1,1]
    A12=as.matrix(siginv[[i]][1,(2:ncol(siginv[[i]]))])
    A22=as.matrix(siginv[[i]][(2:ncol(siginv[[i]])),(2:ncol(siginv[[i]]))])
    for(j in 1:nrow(a)){
      for(k in 1:nrow(a)){
        if(j==k){
          a[j,k]=b[j]^2*a11+2*b[j]*A12[j]+A22[j,j]
        }
        if(j!=k){
          a[j,k]=b[j]*b[k]*a11+b[j]*A12[k]+b[k]*A12[j]+A22[j,k]
        }
      }
    }
    C[[i]]=a
  }
  C1=bdiag(C)
  A=betaxvalid[1,]%*%t(betaxvalid[1,])*siginv[[1]][1,1]
  for(i in 2:length(validuse)){
    a11=siginv[[i]][1,1]
    A=A+betaxvalid[i,]%*%t(betaxvalid[i,])*a11
  }
  Info1=rbind(C1,t(B1))
  Info2=rbind(B1,A)
  Info3=cbind(Info1,Info2)
  Info=solve(Info3)
  if(min(eigen(Info)$values)<0){
    Info=Info+(-1.001)*min(eigen(Info)$values)*diag(1,nrow=nrow(Info),ncol=ncol(Info))
  }
  Covb=Info[(nrow(Info)-d+1):nrow(Info),(ncol(Info)-d+1):ncol(Info)]
  niv=length(validuse)
  W=cbind(b,diag(1,nrow=d,ncol = d))
  V=list()
  for(i1 in 1:niv){
    V[[i1]]=cbind(betaxIVvalid[i1,],diag(0,nrow=d,ncol=d))
  }
  Omega1=list()
  for(i1 in 1:niv){
    Omega1[[i1]]=W%*%siginv[[i1]]
    if(i1==1){
      Omega2=V[[i1]]%*%siginv[[i1]]
    }
    if(i1>1){
      Omega2=cbind(Omega2,V[[i1]]%*%siginv[[i1]])
    }
  }
  Omega3=rbind(bdiag(Omega1),Omega2)
  Omega=-Info[(nrow(Info)-d+1):nrow(Info),]%*%Omega3
  return(list(b=b,Covb=Covb,betaxvalid=betaxvalid,betaxIVvalid=betaxIVvalid,validuse=validuse,Omega=Omega))
}  


for(u in 1:Ncov){
  
  if(u==2){
    d=u
    usecov=1:numcov[u]
    d=numcov[u]
    betaGX1= as.matrix(betaGX1total[,usecov])
    pvaluebetaGX1=as.matrix(pvaluebetaGX1total[,usecov])
    sdbetaGX1=as.matrix(sdbetaGX1total[,usecov])
    ZX1=as.matrix(ZX1total[,usecov])
    betaGYC=betaGYCtotal[,u]
    Nsnps=length(betaGYC)
    pvaluebetaGYC=pvaluebetaGYCtotal[,u]
    sdbetaGYC=sdbetaGYCtotal[,u]
    ZYC=ZYCtotal[,u]
    correlation1=corr_total[[u]][c(1,1+usecov),c(1,1+usecov)]
    gwas_X=list()
    for(i1 in 1:numcov[u]){
      gwas_X[[i1]]=data.frame(rs=ID,beta=betaGX1[,i1],s.d.=sdbetaGX1[,i1],Tstat=ZX1[,i1],p=pvaluebetaGX1[,i1])
      row.names(gwas_X[[i1]])=ID
    }
    gwas_y=data.frame(rs=ID,beta=betaGYC,s.d.=sdbetaGYC,Tstat=ZYC,p=pvaluebetaGYC)
    row.names(gwas_y)=ID
    if(numcov[u]==1){
      pvalX=5e-10
      IV=(pvaluebetaGX1[,1]<pvalX&(ID %in% IVID))
    }
    
    if(numcov[u]==2){
      pvalX=5e-10
      IV=(pvaluebetaGX1[,1]<pvalX&(ID %in% IVID))|(pvaluebetaGX1[,2]<pvalX&(ID %in% IVID))
    }
    
    if(numcov[u]>=3&numcov[u]<11){
      pvalX=5e-10
      IV=(pvaluebetaGX1[,1]<pvalX&(ID %in% IVID))|(pvaluebetaGX1[,2]<pvalX&(ID %in% IVID))
      for(i in 3:ncol(pvaluebetaGX1)){
        IV=IV|(pvaluebetaGX1[,i]<pvalX&(ID %in% IVID))
      }
    }
    
    if(numcov[u]>=11){
      pvalX=5e-10
      IV=(pvaluebetaGX1[,1]<pvalX&(ID %in% IVID))|(pvaluebetaGX1[,2]<pvalX&(ID %in% IVID))
      for(i in 3:ncol(pvaluebetaGX1)){
        IV=IV|(pvaluebetaGX1[,i]<pvalX&(ID %in% IVID))
      }
    }
    
    
    sum(IV)
    used=gene1$V3[IV]
    betaGX_IV2=as.matrix(betaGX1[IV,])
    betaGYC_IV=betaGYC[IV]
    sdbetaGX_IV2=as.matrix(sdbetaGX1[IV,])
    sdbetaGYC_IV=sdbetaGYC[IV]
    pvaluebetaGX_IV2=as.matrix(pvaluebetaGX1[IV,])
    pvaluebetaGYC_IV=pvaluebetaGYC[IV ]
    rs_IV=gene1$V3[IV]
    
    SIG=list()
    for(i in 1:sum(IV)){
      SIG[[i]]=diag(c(sdbetaGYC_IV[i],sdbetaGX_IV2[i,]))%*%correlation1%*%diag(c(sdbetaGYC_IV[i],sdbetaGX_IV2[i,]))
    }
    SIGU1=list()
    SIGU2=list()
    d=1
    for(i in 1:length(SIG)){
      SIGU1[[i]]=SIG[[i]][(1:2),(1:2)]
      SIGU2[[i]]=SIG[[i]][(-2),(-2)]
    }
    
    if(10>0){
    betaGX_IV=as.matrix(betaGX_IV2[,1])
    sdbetaGX_IV=as.matrix(sdbetaGX_IV2[,1])
    pvaluebetaGX_IV=as.matrix(pvaluebetaGX_IV2[,1])
    
    Want=cml_MA(N=Nsample,p=ncol(betaGX_IV ),m=nrow(betaGX_IV ),betax=betaGX_IV ,betayc=betaGYC_IV,SIG1=SIGU1,threshold=0.01,maxit=500)
    b_cml_uv1=Want$b
    covb_cml_uv1=Want$Covb
    library(MendelianRandomization) 
    library(indexevent)
    library(SlopeHunter)
    DP=200
    bDP_egger=matrix(0,nrow=DP,ncol=1)
    bDP_ivw=matrix(0,nrow=DP,ncol=1)
    bDP_lasso=matrix(0,nrow=DP,ncol=1)
    bDP_median=matrix(0,nrow=DP,ncol=1)
    bDP_DD=matrix(0,nrow=DP,ncol=1)
    bDP_SH=matrix(0,nrow=DP,ncol=1)
    
    
    for(j in 1: DP){
      betaX=matrix(0,nrow=nrow(betaGX_IV  ),ncol = d)
      betaY=numeric()
      for(i in 1:nrow(betaGX_IV  )){
        
        beta=mvrnorm(n=1,mu=c((betaGYC_IV  )[i],(betaGX_IV  )[i,]),Sigma = (SIGU1 )[[i]])
        betaX[i,]=beta[2:length(beta)]
        betaY[i]=beta[1]
      }
      inp=mr_mvinput(bx=betaX,bxse =sdbetaGX_IV ,by=betaY,byse =sdbetaGYC_IV )
      inp1=mr_input(bx=as.vector(betaX),bxse =as.vector(sdbetaGX_IV ),by=betaY,byse =sdbetaGYC_IV )
      egger=mr_mvegger(object = inp) 
      ivw=mr_mvivw(object = inp)
      lasso=mr_mvlasso(object = inp,lambda=seq(from=0,to=4,by=0.01))
      median1=mr_median(object = inp1,iterations=100)
      
      hunter=data.frame(xbeta=betaX, xse=sdbetaGX_IV , ybeta=betaY,  yse=sdbetaGYC_IV ,rs=rs_IV  )
      Hunter=hunt(dat=hunter,snp_col = "rs",xbeta_col = "xbeta",xse_col = "xse", ybeta_col = "ybeta",yse_col = "yse",xp_thresh = 0.001,init_pi = 0.6,init_sigmaIP = 1e-05,Bootstrapping = FALSE,M = 200,seed = 777,Plot = FALSE,show_adjustments = FALSE)
      DD=indexevent(seed = 2018,lambda = seq(0.25, 5, 0.25),B = 10, method = c("Hedges-Olkin"),prune = NULL,xbeta=betaX, xse=sdbetaGX_IV ,ybeta=betaY,yse=sdbetaGYC_IV ,weighted = T)
      bDP_egger[j,]=egger$Estimate
      bDP_ivw[j,]=ivw$Estimate
      bDP_lasso[j,]=lasso$Estimate
      bDP_median[j,]=median1$Estimate
      bDP_DD[j,]=DD$b
      bDP_SH[j,]=Hunter$b
    }
    b_egger_uv1=colMeans(bDP_egger)
    covb_egger_uv1=cov(bDP_egger)
    b_ivw_uv1=colMeans(bDP_ivw)
    covb_ivw_uv1=cov(bDP_ivw)
    b_lasso_uv1=colMeans(bDP_lasso)
    covb_lasso_uv1=cov(bDP_lasso) 
    b_median_uv1=colMeans(bDP_median)
    covb_median_uv1=cov(bDP_median) 
    b_DD_uv1=colMeans(bDP_DD)
    covb_DD_uv1=cov(bDP_DD) 
    b_SH_uv1=colMeans(bDP_SH)
    covb_SH_uv1=cov(bDP_SH) 
    }
    if(10>0){
      betaGX_IV=as.matrix(betaGX_IV2[,2])
      sdbetaGX_IV=as.matrix(sdbetaGX_IV2[,2])
      pvaluebetaGX_IV=as.matrix(pvaluebetaGX_IV2[,2])
      
      Want=cml_MA(N=Nsample,p=ncol(betaGX_IV ),m=nrow(betaGX_IV ),betax=betaGX_IV ,betayc=betaGYC_IV,SIG1=SIGU1,threshold=0.01,maxit=500)
      b_cml_uv2=Want$b
      covb_cml_uv2=Want$Covb
      library(MendelianRandomization) 
      library(indexevent)
      library(SlopeHunter)
      DP=200
      bDP_egger=matrix(0,nrow=DP,ncol=1)
      bDP_ivw=matrix(0,nrow=DP,ncol=1)
      bDP_lasso=matrix(0,nrow=DP,ncol=1)
      bDP_median=matrix(0,nrow=DP,ncol=1)
      bDP_DD=matrix(0,nrow=DP,ncol=1)
      bDP_SH=matrix(0,nrow=DP,ncol=1)
      
      
      for(j in 1: DP){
        betaX=matrix(0,nrow=nrow(betaGX_IV  ),ncol = d)
        betaY=numeric()
        for(i in 1:nrow(betaGX_IV  )){
          
          beta=mvrnorm(n=1,mu=c((betaGYC_IV  )[i],(betaGX_IV  )[i,]),Sigma = (SIGU1 )[[i]])
          betaX[i,]=beta[2:length(beta)]
          betaY[i]=beta[1]
        }
        inp=mr_mvinput(bx=betaX,bxse =sdbetaGX_IV ,by=betaY,byse =sdbetaGYC_IV )
        inp1=mr_input(bx=as.vector(betaX),bxse =as.vector(sdbetaGX_IV ),by=betaY,byse =sdbetaGYC_IV )
        egger=mr_mvegger(object = inp) 
        ivw=mr_mvivw(object = inp)
        lasso=mr_mvlasso(object = inp,lambda=seq(from=0,to=4,by=0.01))
        median1=mr_median(object = inp1,iterations=100)
        
        hunter=data.frame(xbeta=betaX, xse=sdbetaGX_IV , ybeta=betaY,  yse=sdbetaGYC_IV ,rs=rs_IV  )
        Hunter=hunt(dat=hunter,snp_col = "rs",xbeta_col = "xbeta",xse_col = "xse", ybeta_col = "ybeta",yse_col = "yse",xp_thresh = 0.001,init_pi = 0.6,init_sigmaIP = 1e-05,Bootstrapping = FALSE,M = 200,seed = 777,Plot = FALSE,show_adjustments = FALSE)
        DD=indexevent(seed = 2018,lambda = seq(0.25, 5, 0.25),B = 10, method = c("Hedges-Olkin"),prune = NULL,xbeta=betaX, xse=sdbetaGX_IV ,ybeta=betaY,yse=sdbetaGYC_IV ,weighted = T)
        bDP_egger[j,]=egger$Estimate
        bDP_ivw[j,]=ivw$Estimate
        bDP_lasso[j,]=lasso$Estimate
        bDP_median[j,]=median1$Estimate
        bDP_DD[j,]=DD$b
        bDP_SH[j,]=Hunter$b
      }
      b_egger_uv2=colMeans(bDP_egger)
      covb_egger_uv2=cov(bDP_egger)
      b_ivw_uv2=colMeans(bDP_ivw)
      covb_ivw_uv2=cov(bDP_ivw)
      b_lasso_uv2=colMeans(bDP_lasso)
      covb_lasso_uv2=cov(bDP_lasso) 
      b_median_uv2=colMeans(bDP_median)
      covb_median_uv2=cov(bDP_median) 
      b_DD_uv2=colMeans(bDP_DD)
      covb_DD_uv2=cov(bDP_DD) 
      b_SH_uv2=colMeans(bDP_SH)
      covb_SH_uv2=cov(bDP_SH) 
    }
    
     
    b_cml=c(b_cml_uv1,b_cml_uv2)
    b_egger=c(b_egger_uv1,b_egger_uv2)
    b_ivw=c(b_ivw_uv1,b_ivw_uv2)
    b_lasso=c(b_lasso_uv1,b_lasso_uv2)
    b_median=c(b_median_uv1,b_median_uv2)
    b_DD=c(b_DD_uv1,b_DD_uv2)
    b_SH=c(b_SH_uv1,b_SH_uv2)
    covb_cml=diag(c(covb_cml_uv1,covb_cml_uv2))
    covb_egger=diag(c(covb_egger_uv1[1,1],covb_egger_uv2[1,1]))
    covb_ivw=diag(c(covb_ivw_uv1[1,1],covb_ivw_uv2[1,1])) 
    covb_lasso=diag(c(covb_lasso_uv1[1,1],covb_lasso_uv2[1,1])) 
    covb_median=diag(c(covb_median_uv1[1,1],covb_median_uv2[1,1])) 
    covb_DD=diag(c(covb_DD_uv1[1,1],covb_DD_uv2[1,1])) 
    covb_SH=diag(c(covb_SH_uv1[1,1],covb_SH_uv2[1,1]))
    betaGYadj_cml=numeric()
    betaGYadj_egger=numeric()
    betaGYadj_lasso=numeric()
    betaGYadj_median=numeric()
    betaGYadj_ivw=numeric()
    betaGYadj_DD=numeric()
    betaGYadj_SH=numeric()
    varbetaGYadj_cml1=numeric()
    varbetaGYadj_egger1=numeric()
    varbetaGYadj_lasso1=numeric()
    varbetaGYadj_median1=numeric()
    varbetaGYadj_ivw1=numeric()
    varbetaGYadj_DD1=numeric()
    varbetaGYadj_SH1=numeric() 
    
    
    for(j1 in 1:Nsnps){
      betaGYadj_cml[j1]=betaGYC[j1]-b_cml%*%betaGX1[j1,]
      betaGYadj_egger[j1]=betaGYC[j1]-b_egger%*%betaGX1[j1,]
      betaGYadj_median[j1]=betaGYC[j1]-b_median%*%betaGX1[j1,]
      betaGYadj_lasso[j1]=betaGYC[j1]-b_lasso%*%betaGX1[j1,]
      betaGYadj_ivw[j1]=betaGYC[j1]-b_ivw%*%betaGX1[j1,]
      betaGYadj_DD[j1]=betaGYC[j1]-b_DD%*%betaGX1[j1,]
      betaGYadj_SH[j1]=betaGYC[j1]-b_SH%*%betaGX1[j1,]
      covbetax=diag(sdbetaGX1[1,])%*%correlation1[(2:ncol(correlation1)),(2:ncol(correlation1))]%*%diag(sdbetaGX1[1,])
      
      varbetaGYadj_cml1[j1]=sdbetaGYC[j1]^2+sum(diag(covbetax*covb_cml))+sum(b_cml^2*diag(covbetax))+t(betaGX1[j1,])%*%covb_cml%*%betaGX1[j1,]
      varbetaGYadj_egger1[j1]=sdbetaGYC[j1]^2+sum(diag(covbetax*covb_egger))+sum(b_egger^2*diag(covbetax))+t(betaGX1[j1,])%*%covb_egger%*%betaGX1[j1,]
      varbetaGYadj_ivw1[j1]=sdbetaGYC[j1]^2+sum(diag(covbetax*covb_ivw))+sum(b_ivw^2*diag(covbetax))+t(betaGX1[j1,])%*%covb_ivw%*%betaGX1[j1,]
      varbetaGYadj_lasso1[j1]=sdbetaGYC[j1]^2+sum(diag(covbetax*covb_lasso))+sum(b_lasso^2*diag(covbetax))+t(betaGX1[j1,])%*%covb_lasso%*%betaGX1[j1,]
      varbetaGYadj_median1[j1]=sdbetaGYC[j1]^2+sum(diag(covbetax*covb_median))+sum(b_median^2*diag(covbetax))+t(betaGX1[j1,])%*%covb_median%*%betaGX1[j1,]
      varbetaGYadj_DD1[j1]=sdbetaGYC[j1]^2+sum(diag(covbetax*covb_DD))+sum(b_DD^2*diag(covbetax))+t(betaGX1[j1,])%*%covb_DD%*%betaGX1[j1,]
      varbetaGYadj_SH1[j1]=sdbetaGYC[j1]^2+sum(diag(covbetax*covb_SH))+sum(b_SH^2*diag(covbetax))+t(betaGX1[j1,])%*%covb_SH%*%betaGX1[j1,]
    }
    SDbetaGYadj_cml1=sqrt(varbetaGYadj_cml1) 
    SDbetaGYadj_egger1=sqrt(varbetaGYadj_egger1) 
    SDbetaGYadj_ivw1=sqrt(varbetaGYadj_ivw1) 
    SDbetaGYadj_lasso1=sqrt(varbetaGYadj_lasso1)
    SDbetaGYadj_median1=sqrt(varbetaGYadj_median1)
    SDbetaGYadj_DD1=sqrt(varbetaGYadj_DD1)
    SDbetaGYadj_SH1=sqrt(varbetaGYadj_SH1)
    
    
    
    z_cml1=abs(betaGYadj_cml/SDbetaGYadj_cml1)
    z_egger1=abs(betaGYadj_egger/SDbetaGYadj_egger1)
    z_lasso1=abs(betaGYadj_lasso/SDbetaGYadj_lasso1)
    z_median1=abs(betaGYadj_median/SDbetaGYadj_median1)
    z_ivw1=abs(betaGYadj_ivw/SDbetaGYadj_ivw1)
    z_DD1=abs(betaGYadj_DD/SDbetaGYadj_DD1)
    z_SH1=abs(betaGYadj_SH/SDbetaGYadj_SH1)
    
    pvalueadj_cml1=(1-pnorm(z_cml1))*2
    pvalueadj_egger1=(1-pnorm(z_egger1))*2
    pvalueadj_median1=(1-pnorm(z_median1))*2
    pvalueadj_ivw1=(1-pnorm(z_ivw1))*2
    pvalueadj_lasso1=(1-pnorm(z_lasso1))*2
    pvalueadj_DD1=(1-pnorm(z_DD1))*2
    pvalueadj_SH1=(1-pnorm(z_SH1))*2
    
    
    
    
    
    
    
    
    
    resultUVMRcML=data.frame(chr=gene1$V1,rs=gene1$V3,beta=betaGYC,sd=sdbetaGYC,p=pvaluebetaGYC,beta_adj=betaGYadj_cml,sd_adj=SDbetaGYadj_cml1,p_adj=pvalueadj_cml1)
    resultUVMREgger=data.frame(chr=gene1$V1,rs=gene1$V3,beta=betaGYC,sd=sdbetaGYC,p=pvaluebetaGYC,beta_adj=betaGYadj_egger,sd_adj=SDbetaGYadj_egger1,p_adj=pvalueadj_egger1) 
    resultUVMRLasso=data.frame(chr=gene1$V1,rs=gene1$V3,beta=betaGYC,sd=sdbetaGYC,p=pvaluebetaGYC,beta_adj=betaGYadj_lasso,sd_adj=SDbetaGYadj_lasso1,p_adj=pvalueadj_lasso1) 
    resultUVMRMedian=data.frame(chr=gene1$V1,rs=gene1$V3,beta=betaGYC,sd=sdbetaGYC,p=pvaluebetaGYC,beta_adj=betaGYadj_median,sd_adj=SDbetaGYadj_median1,p_adj=pvalueadj_median1) 
    resultUVMRIVW=data.frame(chr=gene1$V1,rs=gene1$V3,beta=betaGYC,sd=sdbetaGYC,p=pvaluebetaGYC,beta_adj=betaGYadj_ivw,sd_adj=SDbetaGYadj_ivw1,p_adj=pvalueadj_ivw1) 
    resultDD=data.frame(chr=gene1$V1,rs=gene1$V3,beta=betaGYC,sd=sdbetaGYC,p=pvaluebetaGYC,beta_adj=betaGYadj_DD,sd_adj=SDbetaGYadj_DD1,p_adj=pvalueadj_DD1) 
    resultSH=data.frame(chr=gene1$V1,rs=gene1$V3,beta=betaGYC,sd=sdbetaGYC,p=pvaluebetaGYC,beta_adj=betaGYadj_SH,sd_adj=SDbetaGYadj_SH1,p_adj=pvalueadj_SH1) 
    
    
    slope=matrix(0,ncol=7*2,nrow=1)
    slope_se=matrix(0,ncol=7*2,nrow=1)
    for(i in 1:nrow(slope)){
      slope[i,]=c( b_cml,b_egger,b_ivw,b_lasso,b_median,b_DD,b_SH)
      slope_se[i,]=c( sqrt(diag(covb_cml)),sqrt(diag(covb_egger)),sqrt(diag(covb_ivw)),sqrt(diag(covb_lasso)),sqrt(diag(covb_median)),sqrt(diag(covb_DD)),sqrt(diag(covb_SH)))
    }
    
   
    name1=numeric()
    name2=numeric()
    name3=numeric()
    name4=numeric()
    name5=numeric()
    name6=numeric()
    name7=numeric()
    
    
    for(i in 1:2){
      name1[i]=paste("cml_slop",i,sep="")
      name2[i]=paste("egger_slop",i,sep="")
      name3[i]=paste("ivw_slop",i,sep="")
      name4[i]=paste("lasso_slop",i,sep="")
      name5[i]=paste("median_slop",i,sep="")
      name6[i]=paste("DD_slop",i,sep="")
      name7[i]=paste("SH_slop",i,sep="")
    }
    name1se=numeric()
    name2se=numeric()
    name3se=numeric()
    name4se=numeric()
    name5se=numeric()
    name6se=numeric()
    name7se=numeric()
    
    
    for(i in 1:2){
      
      name1se[i]=paste("cml_slop_se",i,sep="")
      name2se[i]=paste("egger_slop_se",i,sep="")
      name3se[i]=paste("ivw_slop_se",i,sep="")
      name4se[i]=paste("lasso_slop_se",i,sep="")
      name5se[i]=paste("median_slop_se",i,sep="")
      name6se[i]=paste("DD_slop_se",i,sep="")
      name7se[i]=paste("SH_slop_se",i,sep="")
    }
    slope=as.data.frame(slope)
    colnames(slope)=c(name1, name2,name3,name4,name5,name6,name7)
    slope_se=as.data.frame(slope_se)
    colnames(slope_se)=c(name1se, name2se,name3se,name4se,name5se,name6se,name7se)
    
     
    
    
    write.table(resultUVMRcML,file=paste("resultUVMRcML",u,".txt",sep=""),sep="\t",col.names = TRUE,row.names = FALSE,quote = FALSE)
    write.table(resultUVMRLasso,file=paste("resultUVMRLasso",u,".txt",sep=""),sep="\t",col.names = TRUE,row.names = FALSE,quote = FALSE)
    write.table(resultUVMRMedian,file=paste("resultUVMRMedian",u,".txt",sep=""),sep="\t",col.names = TRUE,row.names = FALSE,quote = FALSE)
    write.table(resultUVMRIVW,file=paste("resultUVMRIVW",u,".txt",sep=""),sep="\t",col.names = TRUE,row.names = FALSE,quote = FALSE)
    write.table(resultUVMREgger,file=paste("resultUVMREgger",u,".txt",sep=""),sep="\t",col.names = TRUE,row.names = FALSE,quote = FALSE)
    write.table(resultDD,file=paste("resultDD",u,".txt",sep=""),sep="\t",col.names = TRUE,row.names = FALSE,quote = FALSE)
    write.table(resultSH,file=paste("resultSH",u,".txt",sep=""),sep="\t",col.names = TRUE,row.names = FALSE,quote = FALSE)
    slopetotal=list(slope=slope,slope_se=slope_se )
    save(slopetotal,file = paste("BMI_M1_UVMR_second_turn_slope",u,".RData",sep=""))
    
  }
   
}