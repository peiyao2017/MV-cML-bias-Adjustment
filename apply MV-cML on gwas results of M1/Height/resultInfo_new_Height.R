setwd('/panfs/jay/groups/20/panwei/wan01299/collider_bias/Height')


library(mvtnorm)
library(doParallel)
library(MASS)
library(Matrix)

gene=read.table(file="110K_QCed0.001.bim",header = FALSE)
gene2=read.table(file='/panfs/jay/groups/20/panwei/wan01299/collider_bias/M2/Height/110K_QCed1.fam',header = FALSE)
Ncov=10
IVID=gene$V2
a=list()
b=read.table("response.Height.glm.linear")
Nsample=nrow(gene2)
nameX=numeric()
for(i in 1:20){
  nameX[i]=paste("/panfs/jay/groups/20/panwei/wan01299/collider_bias/M2/Height/MET.PC",i,".glm.linear",sep="")
}

for(i in 1:20){
  a[[i]]=read.table(file=nameX[i],header = FALSE)
}
gene1=a[[1]]

betaGX1total=matrix(0,nrow=nrow(a[[1]]),ncol=length(a))
pvaluebetaGX1total=matrix(0,nrow=nrow(a[[1]]),ncol=length(a))
sdbetaGX1total=matrix(0,nrow=nrow(a[[1]]),ncol=length(a))

betaGYC=b$V9
pvaluebetaGYC=b$V12
sdbetaGYC=b$V10
ID=a[[1]]$V3
nsnpstotal=nrow(a[[1]])


covnames=numeric()
for(i in 1:length(a)){
  covnames[i]=paste("PC",i)
}

for(i in 1:ncol(betaGX1total)){
  betaGX1total[,i]=a[[i ]]$V9
  pvaluebetaGX1total[,i]=a[[i ]]$V12
  sdbetaGX1total[,i]=a[[i ]]$V10
  colnames(betaGX1total)=covnames
  colnames(pvaluebetaGX1total)=covnames
  colnames(sdbetaGX1total)=covnames
}

correlation=diag(1, nrow=length(a)+1,ncol=length(a)+1)

snpx=list()
for(i in 1:length(a)){
  snpx[[i]]=a[[i]]$V3[a[[i]]$V12>0.1586553]
}
snpy=b$V3[b$V12>0.1586553]


for(i in 1:(ncol(correlation)-1)){
  snps=intersect(snpx[[i]],snpy) 
  snps=intersect(snps,IVID)
  zx=a[[i]]$V11[a[[i]]$V3 %in% snps]
  zy=b$V11[b$V3 %in% snps]
  correlation[1,(i+1)]=cor(zx,zy)
  correlation[(i+1),1]=cor(zx,zy)
}



numcov=1:Ncov
name=numeric()
name1=numeric()
for(i in 1:length(numcov)){
  name[i]=paste("resultMAInfo",numcov[i],".txt",sep="")
  name1[i]=paste("slope_and_covariate_effect",numcov[i],".RData",sep="")
}

Nsig=numeric()
for(i in 1:ncol(pvaluebetaGX1total)){
  Nsig[i]=sum(pvaluebetaGX1total[,i]<=5e-8)
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

cml_MA=function(N=Nsample,p=ncol(betaGX1),m=nrow(betaGX_IV),betax=betaGX_IV,betayc=betaGYC_IV,SIG1=SIG,threshold=0.01,maxit=5000){
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
  Info1=rbind(A,B1)
  Info2=rbind(t(B1),C1)
  Info3=cbind(Info1,Info2)
  Info=solve(Info3)
  Covb=Info[(1:length(b)),(1:length(b))]
  return(list(b=b,Covb=Covb))
}

for(u in 10:length(numcov)){
  usecov=c(1:numcov[u])
  betaGX1= as.matrix(betaGX1total[,usecov])
  pvaluebetaGX1=as.matrix(pvaluebetaGX1total[,usecov])
  sdbetaGX1=as.matrix(sdbetaGX1total[,usecov])
  correlation1=correlation[c(1,1+usecov),c(1,1+usecov)]
  
  
  
  if(numcov[u]==1){
    pvalX=5e-8
    IV=(pvaluebetaGX1[,1]<pvalX&(ID %in% IVID))
  }
  if(numcov[u]==2){
    pvalX=5e-8
    IV=(pvaluebetaGX1[,1]<pvalX&(ID %in% IVID))|(pvaluebetaGX1[,2]<pvalX&(ID %in% IVID))
  }
  
  if(numcov[u]>=3&numcov[u]<11){
    pvalX=5e-7
    IV=(pvaluebetaGX1[,1]<pvalX&(ID %in% IVID))|(pvaluebetaGX1[,2]<pvalX&(ID %in% IVID))
    for(i in 3:ncol(pvaluebetaGX1)){
      IV=IV|(pvaluebetaGX1[,i]<pvalX&(ID %in% IVID))
    }
  }
  
  if(numcov[u]>=11){
    pvalX=1e-5    
    IV=(pvaluebetaGX1[,1]<pvalX&(ID %in% IVID))|(pvaluebetaGX1[,2]<pvalX&(ID %in% IVID))
    for(i in 3:ncol(pvaluebetaGX1)){
      IV=IV|(pvaluebetaGX1[,i]<pvalX&(ID %in% IVID))
    }
  }
  
  
  sum(IV)
  
  betaGX_IV=as.matrix(betaGX1[IV,])
  betaGYC_IV=betaGYC[IV]
  sdbetaGX_IV=as.matrix(sdbetaGX1[IV,])
  sdbetaGYC_IV=sdbetaGYC[IV]
  
  
  
  SIG=list()
  for(i in 1:sum(IV)){
    SIG[[i]]=diag(c(sdbetaGYC_IV[i],sdbetaGX_IV[i,]))%*%correlation1%*%diag(c(sdbetaGYC_IV[i],sdbetaGX_IV[i,]))
  }
  
  
  
  
  Want=cml_MA(N=Nsample,p=ncol(betaGX1),m=nrow(betaGX_IV),betax=betaGX_IV,betayc=betaGYC_IV,SIG1=SIG,threshold=0.01,maxit=500)
  
  b=Want$b
  covb=Want$Covb
  
  
  betaGYadj=numeric()
  for(i in 1:nsnpstotal){
    betaGYadj[i]=betaGYC[i]-t(b)%*%betaGX1[i,]
  }
  
  sdbetaGYadj=numeric()
  for(i in 1:nsnpstotal){
    a=sdbetaGYC[i]^2
    b1=t(betaGX1[i,])%*%covb%*%betaGX1[i,]
    c=sum(b^2*as.vector(sdbetaGX1[i,]^2))
    d=sum(diag(covb)*sdbetaGX1[i,]^2)
    sdbetaGYadj[i]=sqrt(a+b1+c+d)
  }
  
  z=abs(betaGYadj/sdbetaGYadj)
  pvalueadj=(1-pnorm(z))*2
  
  result=data.frame(chr=gene1$V1,rs=gene1$V3,beta=betaGYC,sd=sdbetaGYC,p=pvaluebetaGYC,beta_adj=betaGYadj,sd_adj=sdbetaGYadj,p_adj=pvalueadj)
  slope_and_covariate_effect_info=list(b=b, betaGX1=betaGX1,covb=covb)
  write.table(result,file = name[u],sep="\t")
  save(slope_and_covariate_effect_info,file=name1[u])
  print(u)
}