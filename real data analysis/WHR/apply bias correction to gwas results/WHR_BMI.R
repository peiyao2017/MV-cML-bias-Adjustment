setwd("/panfs/jay/groups/20/panwei/wan01299/collider_bias/BMI/MV_Egger/new/new_gwas")

Ncov=1
d=Ncov
library(mvtnorm)
library(doParallel)
number_of_cores=detectCores()
myCluster <- makeCluster(number_of_cores-12, # number of cores to use
                         type = "PSOCK") # type of cluster
#notice that we need to leave one core for computer system
registerDoParallel(myCluster)

library(MASS)
library(Matrix)
library(Matrix)
library(BEDMatrix)
G=BEDMatrix(path="110K_QCed1.bed",simple_names = TRUE)
gene=read.table(file="110K_QCed0.001.bim",header = FALSE)
gene2=read.table(file='110K_QCed1.fam',header = FALSE)
gene1=read.table(file='110K_QCed1.bim',header = FALSE)
IVID=gene$V2
a=read.table("WHR_BMI.WHR.glm.linear")
b1=read.table("BMI_marginal.BMI.glm.linear")

atotal=list(b1 )
BP=a$V2
Nsample=109996


chr=a$V1
betaGX1total=matrix(0,nrow=nrow(a),ncol=Ncov)
pvaluebetaGX1total=matrix(0,nrow=nrow(a),ncol=Ncov)
sdbetaGX1total=matrix(0,nrow=nrow(a),ncol=Ncov)
ZX1total=matrix(0,nrow=nrow(a),ncol=Ncov)

betaGYCtotal=matrix(0,nrow=nrow(a),ncol = 1)
pvaluebetaGYCtotal=matrix(0,nrow=nrow(a),ncol = 1)
sdbetaGYCtotal=matrix(0,nrow=nrow(a),ncol = 1)
ZYCtotal=matrix(0,nrow=nrow(a),ncol = 1)

for(i in 1:1){
  betaGYCtotal[,i]=a$V9
  pvaluebetaGYCtotal[,i]=a$V12
  sdbetaGYCtotal[,i]=a$V10
  ZYCtotal[,i]=a$V11
}


ID=atotal[[1]]$V3
nsnpstotal=nrow(atotal[[1]])

for(i in 1:ncol(betaGX1total)){
  betaGX1total[,i]=atotal[[i ]]$V9
  pvaluebetaGX1total[,i]=atotal[[i ]]$V12
  sdbetaGX1total[,i]=atotal[[i ]]$V10
  ZX1total[,i]=atotal[[i ]]$V11
}


correlation1=diag(1,nrow=Ncov+1,ncol=Ncov+1)


 
BMI=read.table(file="CovWHR_BMI.txt",header = TRUE)$BMI
WHR=read.table(file="WHR.txt",header = TRUE)$WHR 

for(i in 1:ncol(correlation1)){
  for(j in 1:nrow(correlation1)){
    if(i==1&j>1){
      correlation1[i,j]=cor(WHR,BMI)
    }
    if(i>1&j==1){
      correlation1[i,j]=cor(WHR,BMI)
    }
  }
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
  d=p
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



betaGX1= as.matrix(betaGX1total )
pvaluebetaGX1=as.matrix(pvaluebetaGX1total )
sdbetaGX1=as.matrix(sdbetaGX1total )
ZX1=as.matrix(ZX1total )
betaGYC=betaGYCtotal 
pvaluebetaGYC=pvaluebetaGYCtotal 
sdbetaGYC=sdbetaGYCtotal 
ZYC=ZYCtotal 

gwas_X=list()
for(i1 in 1:Ncov){
  gwas_X[[i1]]=data.frame(rs=ID,beta=betaGX1[,i1],s.d.=sdbetaGX1[,i1],Tstat=ZX1[,i1],p=pvaluebetaGX1[,i1])
  row.names(gwas_X[[i1]])=ID
}
gwas_y=data.frame(rs=ID,beta=betaGYC,s.d.=sdbetaGYC,Tstat=ZYC,p=pvaluebetaGYC)
row.names(gwas_y)=ID
if(Ncov==1){
  pvalX=5e-7
  IV=(pvaluebetaGX1[,1]<pvalX&(ID %in% IVID))
}
if(Ncov==2){
  pvalX=5e-7
  IV=(pvaluebetaGX1[,1]<pvalX&(ID %in% IVID))|(pvaluebetaGX1[,2]<pvalX&(ID %in% IVID))
}


sum(IV)
used=gene1$V2[IV]
rs_IV=gene1$V2[IV]
betaGX_IV=as.matrix(betaGX1[IV,])
betaGYC_IV=betaGYC[IV]
sdbetaGX_IV=as.matrix(sdbetaGX1[IV,])
sdbetaGYC_IV=sdbetaGYC[IV]
pvaluebetaGX_IV=as.matrix(pvaluebetaGX1[IV,])
pvaluebetaGYC_IV=pvaluebetaGYC[IV]


SIG=list()
for(i in 1:sum(IV)){
  SIG[[i]]=diag(c(sdbetaGYC_IV[i],sdbetaGX_IV[i,]))%*%correlation1%*%diag(c(sdbetaGYC_IV[i],sdbetaGX_IV[i,]))
}




Want=cml_MA(N=Nsample,p=ncol(betaGX1),m=nrow(betaGX_IV),betax=betaGX_IV,betayc=betaGYC_IV,SIG1=SIG,threshold=0.001,maxit=500)

b=Want$b
covb=Want$Covb
if(min(eigen(covb)$values)<0){
  covb=covb+diag((-1.001)*min(eigen(covb)$values),nrow=d,ncol = d)
}
valid=Want$validuse
betaxvalid=Want$betaxvalid
betaxIVvalid=Want$betaxIVvalid
used=used[valid]
Omega=Want$Omega
print("function finish")

betaGYadj=numeric()
for(i in 1:nsnpstotal){
  betaGYadj[i]=betaGYC[i]-t(b)%*%betaGX1[i,]
  print(i)
}

sdbetaGYadj2=numeric()
ld=matrix(0,nrow=nsnpstotal,ncol=length(used))
row.names(ld)=ID
colnames(ld)=used


blocks=read.table("fourier_ls_all.bed.txt",header = TRUE)
blocks_used=list()
chr_used=chr[ID%in%used]
for(i in 1:length(used)){
  blocks_used[[i]]=numeric()
  blockstemp=blocks[blocks$chr==chr_used[i],]
  starttemp=blockstemp$start
  stoptemp=blockstemp$stop
  snpstemp=ID[chr==chr_used[i]]
  BP_temp=BP[ID%in%snpstemp]
  BP_tar=BP[ID==used[i]]
  for(i1 in 1:length(starttemp)){
    if(BP_tar>=starttemp[i1]&BP_tar<=stoptemp[i1]){
      startar=starttemp[i1]
      stoptar=stoptemp[i1]
    }
  }
  for(i1 in 1:length(snpstemp)){
    if(BP_temp[i1]>=startar&BP_temp[i1]<=stoptar){
      blocks_used[[i]]=c(blocks_used[[i]],snpstemp[i1])
    }
  }
  print(i)
}

ld=list()
for(i in 1:length(used)){
  ld[[i]]=numeric()
  for(i1 in 1:length(blocks_used[[i]])){
    ld[[i]][i1]=cor(G[,used[i]],G[,blocks_used[[i]][i1]],use = "pairwise.complete.obs")
    
  }
  names(ld[[i]])=blocks_used[[i]]
  print(i)
}




gwas_XIV=list()
for(i1 in 1:Ncov){
  gwas_XIV[[i1]]=data.frame(rs=used,beta=betaGX1[,i1][ID%in%used],s.d.=sdbetaGX1[,i1][ID%in%used],Tstat=ZX1[,i1][ID%in%used],p=pvaluebetaGX1[,i1][ID%in%used])
  row.names(gwas_XIV[[i1]])=used
}
gwas_yIV=data.frame(rs=used,beta=betaGYC[ID%in%used],s.d.=sdbetaGYC[ID%in%used],Tstat=ZYC[ID%in%used],p=pvaluebetaGYC[ID%in%used])
row.names(gwas_yIV)=used


niv=length(used)

varbeta3=foreach(i=1:nsnpstotal,.combine = "c")%dopar%{
  library(Matrix)
  ldi=rep(0,times=niv)
  for(i1 in 1:niv){
    if(ID[i]%in%names(ld[[i1]])){
      ldi[i1]=ld[[i1]][ID[i]]
    }
  }
  rs=0
  beta=0
  s.d.=0
  Tstat=0
  p=0
  GX=data.frame(rs,beta,s.d.,Tstat,p)
  GY=GX
  ZY=GX
  ZX=list()
  for(i1 in 1:length(used)){
    ZX[[i1]]=GX
    for(i2 in 1:d){
      ZX[[i1]][i2,]=gwas_XIV[[i2]][used[i1],]
    }
    ZY[i1,]=gwas_yIV[used[i1],]
  }
  
  for(i1 in 1:d){
    GX[i1,]=gwas_X[[i1]][i,]
  }
  GY[1,]=gwas_y[i,]
  sigmaG=diag(c(GY$s.d.,GX$s.d.))%*%correlation1%*%diag(c(GY$s.d.,GX$s.d.))
  sigmaGZY=list()
  for(i1 in 1:niv){
    sigmaGZY[[i1]]=matrix(ldi[i1]*c(GY$s.d.,GX$s.d.)*ZY$s.d.[i1]*correlation1[1,],nrow=d+1,ncol=1)
  }
  sigmaGZH=list()
  for(i1 in 1:niv){
    sigmaGZH[[i1]]=matrix(0,nrow=d+1,ncol=d)
    for(i2 in 1:d+1){
      if(i2==1){
        sigmaGZH[[i1]][i2,]=ldi[i1]*GY$s.d.*ZX[[i1]]$s.d.*correlation1[1,(2:(d+1))]
      }
      if(i2>1){
        sigmaGZH[[i1]][i2,]=ldi[i1]*GX$s.d.[i2-1]*ZX[[i1]]$s.d.*correlation1[i2,(2:(d+1))]
      }
    }
  }
  sigmaZ=bdiag(SIG[valid])
  for(i1 in 1:length(used)){
    if(i1==1){
      sigma12=cbind(sigmaGZY[[i1]],sigmaGZH[[i1]])
    }
    if(i1>1){
      sigma12=cbind(sigma12,sigmaGZY[[i1]],sigmaGZH[[i1]])
    }
  }
  SIGMA1=cbind(sigmaG,sigma12)
  SIGMA2=cbind(t(sigma12),sigmaZ)
  SIGMA=rbind(SIGMA1,SIGMA2)
  U1=bdiag(diag(1,nrow=d+1,ncol=d+1),Omega)%*%SIGMA%*%t(bdiag(diag(1,nrow=d+1,ncol=d+1),Omega))
  v=c(1,-b,-GX$beta)
  if(min(eigen(U1)$values)<0){
    U1=U1+(-min(eigen(U1)$values))*1.001*diag(1,nrow=nrow(U1),ncol=ncol(U1))
  }
  out=as.numeric(t(v)%*%U1%*%v)
}

sdbetaGYadj=sqrt(varbeta3)


for(i in 1:nsnpstotal){
  a=sdbetaGYC[i]^2
  b1=t(betaGX1[i,])%*%covb%*%betaGX1[i,]
  c=sum(b^2*as.vector(sdbetaGX1[i,]^2))
  d1=sum(diag(covb)*sdbetaGX1[i,]^2)
  sdbetaGYadj2[i]=sqrt(a+b1+c+d1)
}

z=abs(betaGYadj/sdbetaGYadj)
pvalueadj=(1-pnorm(z))*2
z2=abs(betaGYadj/sdbetaGYadj2)
pvalueadj2=(1-pnorm(z2))*2
resultMVMRcML=data.frame(chr=gene1$V1,rs=gene1$V2,beta=betaGYC,sd=sdbetaGYC,p=pvaluebetaGYC,beta_adj=betaGYadj,sd_adj=sdbetaGYadj,p_adj=pvalueadj,sd_adj2=sdbetaGYadj2,p_adj2=pvalueadj2)

library(MendelianRandomization)
library(indexevent)
library(SlopeHunter)
d=1 


inp=mr_mvinput(bx=as.matrix(betaGX_IV),bxse =as.matrix(sdbetaGX_IV),by=betaGYC_IV,byse =sdbetaGYC_IV)
inp1=mr_input(bx=as.vector(betaGX_IV),bxse =as.vector(sdbetaGX_IV),by=as.vector(betaGYC_IV),byse =as.vector(sdbetaGYC_IV))
egger=mr_mvegger(object = inp) 
ivw=mr_mvivw(object = inp)
lasso=mr_mvlasso(object = inp,lambda=seq(from=0,to=4,by=0.01))
median1=mr_median(object = inp1)
hunter=data.frame(xbeta=betaGX_IV, xse=sdbetaGX_IV, ybeta=betaGYC_IV,  yse=sdbetaGYC_IV,rs=rs_IV,px=pvaluebetaGX_IV,py=pvaluebetaGYC_IV)
Hunter=hunt(dat=hunter,snp_col = "rs",xbeta_col = "xbeta",xse_col = "xse", ybeta_col = "ybeta",yse_col = "yse",xp_thresh = 0.001,init_pi = 0.6,init_sigmaIP = 1e-05,Bootstrapping = TRUE,M = 200,seed = 777,Plot = FALSE,show_adjustments = FALSE)
DP=200
bDP_DD=matrix(0,nrow=DP,ncol=d)
for(j in 1: DP){
  betaX=matrix(0,nrow=nrow(as.matrix(betaGX_IV)),ncol = d)
  betaY=numeric()
  for(i in 1:nrow(as.matrix(betaGX_IV))){
    beta=mvrnorm(n=1,mu=c((betaGYC)[i],as.matrix(betaGX_IV)[i,]),Sigma = (SIG)[[i]])
    betaX[i,]=beta[2:length(beta)]
    betaY[i]=beta[1]
  }
  DD=indexevent(seed = 2018,lambda = seq(0.25, 10, 0.25),B = 10, method = c("Hedges-Olkin"),prune = NULL,xbeta=betaX, xse=sdbetaGX_IV,ybeta=betaY,yse=sdbetaGYC_IV,weighted = T)
  bDP_DD[j,]=DD$b
}
b_DD=colMeans(bDP_DD)
covb_DD=cov(bDP_DD)





b_egger=egger$Estimate
covb_egger=egger$StdError.Est^2
b_ivw=ivw$Estimate
covb_ivw=ivw$StdError^2
b_lasso=lasso$Estimate
covb_lasso=lasso$StdError^2
b_median=median1$Estimate
covb_median=median1$StdError^2

b_SH=Hunter$b
covb_SH=Hunter$bse^2 

inp=mr_mvinput(bx=as.matrix(betaGX_IV ),bxse =as.matrix(sdbetaGX_IV ),by=betaGYC_IV ,byse =sdbetaGYC_IV )
inp1=mr_input(bx=as.vector(betaGX_IV),bxse =as.vector(sdbetaGX_IV),by=as.vector(betaGYC_IV) ,byse =as.vector(sdbetaGYC_IV) )
egger=mr_mvegger(object = inp) 
ivw=mr_mvivw(object = inp)
lasso=mr_mvlasso(object = inp,lambda=seq(from=0,to=4,by=0.01))
median1=mr_median(object = inp1)
hunter=data.frame(xbeta=betaGX_IV, xse=sdbetaGX_IV, ybeta=betaGYC_IV,  yse=sdbetaGYC_IV,rs=rs_IV,px=pvaluebetaGX_IV,py=pvaluebetaGYC_IV)
Hunter=hunt(dat=hunter,snp_col = "rs",xbeta_col = "xbeta",xse_col = "xse", ybeta_col = "ybeta",yse_col = "yse",xp_thresh = 0.001,init_pi = 0.6,init_sigmaIP = 1e-05,Bootstrapping = TRUE,M = 200,seed = 777,Plot = FALSE,show_adjustments = FALSE)
DD=indexevent(seed = 2018,lambda = seq(0.25, 10, 0.25),B = 10, method = c("Hedges-Olkin"),prune = NULL,xbeta=betaGX_IV, xse=sdbetaGX_IV,ybeta=betaGYC_IV,yse=sdbetaGYC_IV,weighted = T)
b_temp_egger=egger$Estimate
covb_temp_egger=egger$StdError.Est^2
b_temp_ivw=ivw$Estimate
covb_temp_ivw=ivw$StdError^2
b_temp_lasso=lasso$Estimate
covb_temp_lasso=lasso$StdError^2
b_temp_median=median1$Estimate
covb_temp_median=median1$StdError^2
b_temp_DD=DD$b
covb_temp_DD=DD$b.se^2 
b_temp_SH=Hunter$b
covb_temp_SH=Hunter$bse^2 

betaGYadj_egger=numeric()
betaGYadj_SH=numeric()
betaGYadj_DD=numeric()
betaGYadj_lasso=numeric()
betaGYadj_median=numeric()
betaGYadj_ivw=numeric()
varbetaGYadj_egger1=numeric()
varbetaGYadj_SH1=numeric()
varbetaGYadj_DD1=numeric()
varbetaGYadj_lasso1=numeric()
varbetaGYadj_median1=numeric()
varbetaGYadj_ivw1=numeric()
varbetaGYadj_egger2=numeric()
varbetaGYadj_SH2=numeric()
varbetaGYadj_DD2=numeric()
varbetaGYadj_lasso2=numeric()
varbetaGYadj_median2=numeric()
varbetaGYadj_ivw2=numeric()


for(j1 in 1:nsnpstotal){
  betaGYadj_egger[j1]=betaGYC[j1]-b_egger%*%betaGX1[j1,]
  betaGYadj_SH[j1]=betaGYC[j1]-b_SH%*%betaGX1[j1,]
  betaGYadj_DD[j1]=betaGYC[j1]-b_DD%*%betaGX1[j1,]
  betaGYadj_median[j1]=betaGYC[j1]-b_median%*%betaGX1[j1,]
  betaGYadj_lasso[j1]=betaGYC[j1]-b_lasso%*%betaGX1[j1,]
  betaGYadj_ivw[j1]=betaGYC[j1]-b_ivw%*%betaGX1[j1,]
  covbetax=diag(x=c(sdbetaGX1[j1,]),nrow=d,ncol=d)%*%correlation1[2:(1+d),2:(1+d)]%*%diag(x=c(sdbetaGX1[j1,]),nrow=d,ncol=d)
  covbetayx=sdbetaGYC[j1]*sdbetaGX1[j1,]*correlation1[1,2:(d+1)]
  varbetaGYadj_egger1[j1]=sdbetaGYC[j1]^2+sum(covbetax*covb_egger)+t(b_egger)%*%covbetax%*%b_egger+t(betaGX1[j1,])%*%covb_egger%*%betaGX1[j1,]
  varbetaGYadj_SH1[j1]=sdbetaGYC[j1]^2+sum(covbetax*covb_SH)+t(b_SH)%*%covbetax%*%b_SH+t(betaGX1[j1,])%*%covb_SH%*%betaGX1[j1,]
  varbetaGYadj_DD1[j1]=sdbetaGYC[j1]^2+sum(covbetax*covb_DD)+t(b_DD)%*%covbetax%*%b_DD+t(betaGX1[j1,])%*%covb_DD%*%betaGX1[j1,]
  varbetaGYadj_ivw1[j1]=sdbetaGYC[j1]^2+sum(covbetax*covb_ivw)+t(b_ivw)%*%covbetax%*%b_ivw+t(betaGX1[j1,])%*%covb_ivw%*%betaGX1[j1,]
  varbetaGYadj_lasso1[j1]=sdbetaGYC[j1]^2+sum(covbetax*covb_lasso)+t(b_lasso)%*%covbetax%*%b_lasso+t(betaGX1[j1,])%*%covb_lasso%*%betaGX1[j1,]
  varbetaGYadj_median1[j1]=sdbetaGYC[j1]^2+sum(covbetax*covb_median)+t(b_median)%*%covbetax%*%b_median+t(betaGX1[j1,])%*%covb_median%*%betaGX1[j1,]
  varbetaGYadj_egger2[j1]=sdbetaGYC[j1]^2+sum(covbetax*covb_egger)+t(b_egger)%*%covbetax%*%b_egger+t(betaGX1[j1,])%*%covb_egger%*%betaGX1[j1,]-2*sum(b_egger*covbetayx)
  varbetaGYadj_SH2[j1]=sdbetaGYC[j1]^2+sum(covbetax*covb_SH)+t(b_SH)%*%covbetax%*%b_SH+t(betaGX1[j1,])%*%covb_SH%*%betaGX1[j1,]-2*sum(b_SH*covbetayx)
  varbetaGYadj_DD2[j1]=sdbetaGYC[j1]^2+sum(covbetax*covb_DD)+t(b_DD)%*%covbetax%*%b_DD+t(betaGX1[j1,])%*%covb_DD%*%betaGX1[j1,]-2*sum(b_DD*covbetayx)
  varbetaGYadj_ivw2[j1]=sdbetaGYC[j1]^2+sum(covbetax*covb_ivw)+t(b_ivw)%*%covbetax%*%b_ivw+t(betaGX1[j1,])%*%covb_ivw%*%betaGX1[j1,]-2*sum(b_ivw*covbetayx)
  varbetaGYadj_lasso2[j1]=sdbetaGYC[j1]^2+sum(covbetax*covb_lasso)+t(b_lasso)%*%covbetax%*%b_lasso+t(betaGX1[j1,])%*%covb_lasso%*%betaGX1[j1,]-2*sum(b_lasso*covbetayx)
  varbetaGYadj_median2[j1]=sdbetaGYC[j1]^2+sum(covbetax*covb_median)+t(b_median)%*%covbetax%*%b_median+t(betaGX1[j1,])%*%covb_median%*%betaGX1[j1,]-2*sum(b_median*covbetayx)
}

SDbetaGYadj_egger1=sqrt(varbetaGYadj_egger1) 
SDbetaGYadj_SH1=sqrt(varbetaGYadj_SH1)
SDbetaGYadj_DD1=sqrt(varbetaGYadj_DD1)
SDbetaGYadj_egger2=sqrt(varbetaGYadj_egger2) 
SDbetaGYadj_SH2=sqrt(varbetaGYadj_SH2)
SDbetaGYadj_DD2=sqrt(varbetaGYadj_DD2)



z_egger1=abs(betaGYadj_egger/SDbetaGYadj_egger1)
z_DD1=abs(betaGYadj_DD/SDbetaGYadj_DD1)
z_SH1=abs(betaGYadj_SH/SDbetaGYadj_SH1)
z_egger2=abs(betaGYadj_egger/SDbetaGYadj_egger2)
z_DD2=abs(betaGYadj_DD/SDbetaGYadj_DD2)
z_SH2=abs(betaGYadj_SH/SDbetaGYadj_SH2)



pvalueadj_egger1=(1-pnorm(z_egger1))*2
pvalueadj_DD1=(1-pnorm(z_DD1))*2
pvalueadj_SH1=(1-pnorm(z_SH1))*2
pvalueadj_egger2=(1-pnorm(z_egger2))*2
pvalueadj_DD2=(1-pnorm(z_DD2))*2
pvalueadj_SH2=(1-pnorm(z_SH2))*2

SDbetaGYadj_ivw1=sqrt(varbetaGYadj_ivw1) 
SDbetaGYadj_lasso1=sqrt(varbetaGYadj_lasso1)
SDbetaGYadj_median1=sqrt(varbetaGYadj_median1)
SDbetaGYadj_ivw2=sqrt(varbetaGYadj_ivw2) 
SDbetaGYadj_lasso2=sqrt(varbetaGYadj_lasso2)
SDbetaGYadj_median2=sqrt(varbetaGYadj_median2)



z_ivw1=abs(betaGYadj_ivw/SDbetaGYadj_ivw1)
z_median1=abs(betaGYadj_median/SDbetaGYadj_median1)
z_lasso1=abs(betaGYadj_lasso/SDbetaGYadj_lasso1)
z_ivw2=abs(betaGYadj_ivw/SDbetaGYadj_ivw2)
z_median2=abs(betaGYadj_median/SDbetaGYadj_median2)
z_lasso2=abs(betaGYadj_lasso/SDbetaGYadj_lasso2)



pvalueadj_ivw1=(1-pnorm(z_ivw1))*2
pvalueadj_median1=(1-pnorm(z_median1))*2
pvalueadj_lasso1=(1-pnorm(z_lasso1))*2
pvalueadj_ivw2=(1-pnorm(z_ivw2))*2
pvalueadj_median2=(1-pnorm(z_median2))*2
pvalueadj_lasso2=(1-pnorm(z_lasso2))*2

resultMVMREgger=data.frame(chr=gene1$V1,rs=gene1$V2,beta=betaGYC,sd=sdbetaGYC,p=pvaluebetaGYC,beta_adj=betaGYadj_egger,sd_adj=SDbetaGYadj_egger1,p_adj=pvalueadj_egger1,sd_adj2=SDbetaGYadj_egger2,p_adj2=pvalueadj_egger2)
resultMVMRLasso=data.frame(chr=gene1$V1,rs=gene1$V2,beta=betaGYC,sd=sdbetaGYC,p=pvaluebetaGYC,beta_adj=betaGYadj_lasso,sd_adj=SDbetaGYadj_lasso1,p_adj=pvalueadj_lasso1,sd_adj2=SDbetaGYadj_lasso2,p_adj2=pvalueadj_lasso2)
resultMVMRMedian=data.frame(chr=gene1$V1,rs=gene1$V2,beta=betaGYC,sd=sdbetaGYC,p=pvaluebetaGYC,beta_adj=betaGYadj_median,sd_adj=SDbetaGYadj_median1,p_adj=pvalueadj_median1,sd_adj2=SDbetaGYadj_median2,p_adj2=pvalueadj_median2)
resultMVMRIVW=data.frame(chr=gene1$V1,rs=gene1$V2,beta=betaGYC,sd=sdbetaGYC,p=pvaluebetaGYC,beta_adj=betaGYadj_ivw,sd_adj=SDbetaGYadj_ivw1,p_adj=pvalueadj_ivw1,sd_adj2=SDbetaGYadj_ivw2,p_adj2=pvalueadj_ivw2)
resultDD=data.frame(chr=gene1$V1,rs=gene1$V2,beta=betaGYC,sd=sdbetaGYC,p=pvaluebetaGYC,beta_adj=betaGYadj_DD,sd_adj=SDbetaGYadj_DD1,p_adj=pvalueadj_DD1,sd_adj2=SDbetaGYadj_DD2,p_adj2=pvalueadj_DD2)
resultSH=data.frame(chr=gene1$V1,rs=gene1$V2,beta=betaGYC,sd=sdbetaGYC,p=pvaluebetaGYC,beta_adj=betaGYadj_SH,sd_adj=SDbetaGYadj_SH1,p_adj=pvalueadj_SH1,sd_adj2=SDbetaGYadj_SH2,p_adj2=pvalueadj_SH2)






slope=matrix(0,ncol=7*d,nrow=1)
slope_se=matrix(0,ncol=7*d,nrow=1)
for(i in 1:nrow(slope)){
  slope[i,]=c(b,b_egger,b_DD,b_SH,b_ivw,b_lasso,b_median)
  slope_se[i,]=c(sqrt(covb),sqrt(covb_egger),sqrt(covb_DD),sqrt(covb_SH),sqrt(covb_ivw),sqrt(covb_lasso),sqrt(covb_median))
}

slope_temp=matrix(0,ncol=7*d,nrow=1)
slope_se_temp=matrix(0,ncol=7*d,nrow=1)
for(i in 1:nrow(slope)){
  slope_temp[i,]=c(b,b_temp_egger,b_temp_DD,b_temp_SH,b_temp_ivw,b_temp_lasso,b_temp_median)
  slope_se_temp[i,]=c(sqrt(covb),sqrt(covb_temp_egger),sqrt(covb_temp_DD),sqrt(covb_temp_SH),sqrt(covb_temp_ivw),sqrt(covb_temp_lasso),sqrt(covb_temp_median))
}

name2=numeric()
name3=numeric()
name4=numeric()
name5=numeric()
name6=numeric()
name7=numeric()
name8=numeric()
for(i in 1:d){
  
  name2[i]=paste("cml_slop",i,sep="")
  name3[i]=paste("egger_slop",i,sep="")
  name4[i]=paste("DD_slop",i,sep="")
  name5[i]=paste("SH_slop",i,sep="")
  name6[i]=paste("ivw_slop",i,sep="")
  name7[i]=paste("lasso_slop",i,sep="")
  name8[i]=paste("median_slop",i,sep="")
}

name2se=numeric()
name3se=numeric()
name4se=numeric()
name5se=numeric()
name6se=numeric()
name7se=numeric()
name8se=numeric()
for(i in 1:d){
  
  name2se[i]=paste("cml_slop_se",i,sep="")
  name3se[i]=paste("egger_slop_se",i,sep="")
  name4se[i]=paste("DD_slop_se",i,sep="")
  name5se[i]=paste("SH_slop_se",i,sep="")
  name6se[i]=paste("ivw_slop_se",i,sep="")
  name7se[i]=paste("lasso_slop_se",i,sep="")
  name8se[i]=paste("median_slop_se",i,sep="")
}
slope=as.data.frame(slope)
colnames(slope)=c( name2,name3,name4,name5,name6,name7,name8)
slope_se=as.data.frame(slope_se)
colnames(slope_se)=c( name2se,name3se,name4se,name5se,name6se,name7se,name8se)
slope_temp=as.data.frame(slope_temp)
colnames(slope_temp)=c( name2,name3,name4,name5,name6,name7,name8)
slope_se_temp=as.data.frame(slope_se_temp)
colnames(slope_se_temp)=c( name2se,name3se,name4se,name5se,name6se,name7se,name8se)

write.table(resultMVMRcML,file="resultMVMRcML1WHR.txt",sep="\t",col.names = TRUE,row.names = FALSE,quote = FALSE)
write.table(resultMVMRLasso,file="resultMVMRLasso1WHR.txt",sep="\t",col.names = TRUE,row.names = FALSE,quote = FALSE)
write.table(resultMVMRMedian,file="resultMedian1WHR.txt",sep="\t",col.names = TRUE,row.names = FALSE,quote = FALSE)
write.table(resultMVMRIVW,file="resultMVMRIVW1WHR.txt",sep="\t",col.names = TRUE,row.names = FALSE,quote = FALSE)
write.table(resultMVMREgger,file="resultMVMREgger1WHR.txt",sep="\t",col.names = TRUE,row.names = FALSE,quote = FALSE)
write.table(resultDD,file="resultDDWHR.txt",sep="\t",col.names = TRUE,row.names = FALSE,quote = FALSE)
write.table(resultSH,file="resultSHWHR.txt",sep="\t",col.names = TRUE,row.names = FALSE,quote = FALSE)



slopetotal=list(slope=slope,slope_se=slope_se,slope_temp=slope_temp,slope_se_temp=slope_se_temp)
save(slopetotal,file = "BMIslope1.RData")