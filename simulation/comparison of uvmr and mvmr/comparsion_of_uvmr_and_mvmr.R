setwd("/panfs/jay/groups/20/panwei/wan01299/collider_bias/BMI/MV_Egger/new/mvmr_uvmr/")

library(doParallel)
number_of_cores=detectCores()
myCluster <- makeCluster(number_of_cores-10, # number of cores to use
                         type = "PSOCK") # type of cluster
#notice that we need to leave one core for computer system
registerDoParallel(myCluster)
library(mvtnorm)
library(MASS)
d=2
rhos=c(0)
name_effects=numeric()
Nsnps=1000

N=20000
Nsample=20000
for(i in 1:length(rhos)){
  name_effects[i]=paste("effects",rhos[i],"dim",d,".RData",sep="")
}
outnames=matrix(0,nrow=length(rhos),ncol = 20)
for(t in 1:length(rhos)){
  for(t1 in 1:20){
    outnames[t,t1]=paste("sim",rhos[t],"_",t1,"betaXYnot0",".txt",sep="")
  }
}





for(t in 1:length(rhos)){
  load(name_effects[t])
  betaGY=effect[[1]]*c(rep(0,times=50),rep(1,times=50),rep(0,times=900))
  betaGX=matrix(0,nrow=d,ncol = Nsnps)
  for(i in 1:d){
    betaGX[i,]=effect[[2]][i,]*c(rep(1,50),rep(0,50),rep(0,50),rep(0,850))
  }
  MAF=effect[[3]]
  betaXY=effect[[4]]
  SNP_effect=betaGY
  
  for(t1 in 1:20){
    output=foreach(j=1:50,.combine = "rbind")%dopar%{
      library(mvtnorm)
      library(MASS)
      Gmatrix=matrix(0,nrow=N,ncol=Nsnps)
      for(i in 1:ncol(Gmatrix)){
        for(j in 1:nrow(Gmatrix)){
          Gmatrix[j,i]=rbinom(n=1,size = 2,prob = MAF[i])
        }
      }
      LD=diag(1,nrow=nrow(Gmatrix),ncol=ncol(Gmatrix))
      for(i in 1:ncol(Gmatrix)){
        Gmatrix[,i]=scale(Gmatrix[,i],center = TRUE,scale=FALSE)
      }
      varG=numeric()
      for(i in 1:ncol(Gmatrix)){
        varG[i]=var(Gmatrix[,i])
      }
      X=matrix(0,nrow=N,ncol=d)
      betaUX=numeric()
      for(i in 1:d){
        betaUX[i]=sqrt(0.8*sum(betaGX[i,]^2*varG))
      }
      varEX=numeric()
      for(i in 1:d){
        varEX[i]=0.2*sum(betaGX[i,]^2*varG)
      }
      
      U=rnorm(n=N,mean=0,sd=1)
      EX=mvrnorm(n=N,mu=rep(0,times=d),Sigma = diag(varEX))
      for(i in 1:d){
        X[,i]=Gmatrix%*%betaGX[i,]+betaUX[i]*U+EX[,i]
      }
      X_sim=matrix(0,nrow=nrow(X),ncol = ncol(X))
      for(i in 1:d){
        X_sim[,i]=scale(X[,i],center = TRUE,scale=FALSE)
      }
      y_known=Gmatrix%*%betaGY+X_sim%*%betaXY
      betaUY=sqrt(0.8*var(y_known))
      varEY=0.2*var(y_known)
      EY=rnorm(n=N,mean=0,sd=sqrt(varEY))
      y=y_known+betaUY[1,1]*U+EY
      y_sim=scale(y,center = TRUE,scale=FALSE)
      G_sim=Gmatrix
      snp_names=numeric()
      for(i in 1:Nsnps){
        snp_names[i]=paste("rs",i)
      }
      colnames(G_sim)=snp_names
      colnames(Gmatrix)=snp_names
      library(Matrix)  
      gwas_X=list()
      for(i in 1:d){
        a=matrix(0,nrow=ncol(Gmatrix),ncol=5)
        for(j in 1:ncol(Gmatrix)){
          data_X=data.frame(res_X=X_sim[,i],G=G_sim[,j])
          m=lm(res_X~.-1,data=data_X)
          obj=summary(m)
          a[j,1]=colnames(G_sim)[j]
          a[j,2]=m$coefficients[length(m$coefficients)]
          a[j,3]=obj[[4]][,2][length(obj[[4]][,2])]
          a[j,4]=obj[[4]][,3][length(obj[[4]][,3])]
          a[j,5]=obj[[4]][,4][length(obj[[4]][,4])]
          print(j)
        }
        colnames(a)=c("rs","beta","s.d.","t-stat","p")
        gwas_X[[i]]=as.data.frame(a)
        print(i)
      }
      
      
      for(i in 1:d){
        for(j in 2:ncol(gwas_X[[i]])){
          gwas_X[[i]][,j]=as.numeric(gwas_X[[i]][,j])
        }
      }
      
      gwas_y=matrix(0,nrow=ncol(Gmatrix),ncol=5)
      for(j in 1:ncol(Gmatrix)){
        data_y=data.frame(res_y=y_sim,X_sim,G=G_sim[,j])
        m=lm(res_y~.-1,data=data_y)
        obj=summary(m)
        gwas_y[j,1]=colnames(G_sim)[j]
        gwas_y[j,2]=m$coefficients[length(m$coefficients)]
        gwas_y[j,3]=obj[[4]][,2][length(obj[[4]][,2])]
        gwas_y[j,4]=obj[[4]][,3][length(obj[[4]][,3])]
        gwas_y[j,5]=obj[[4]][,4][length(obj[[4]][,4])]
        print(j)
      }
      colnames(gwas_y)=c("rs","beta","s.d.","t-stat","p")
      gwas_y=as.data.frame(gwas_y)
      for(j in 2:ncol(gwas_y )){
        gwas_y[,j]=as.numeric(gwas_y[,j])
      }
      print("data ready")
      Nsample=20000
      nsnpstotal=ncol(G_sim)
      print("sample size")
      
      betaGX1= matrix(0,ncol=d,nrow=nrow(gwas_y))
      for(i in 1:d){
        betaGX1[,i]=gwas_X[[i]]$beta
      }
      pvaluebetaGX1=matrix(0,ncol=d,nrow=nrow(gwas_y))
      for(i in 1:d){
        pvaluebetaGX1[,i]=gwas_X[[i]]$p
      }
      sdbetaGX1=matrix(0,ncol=d,nrow=nrow(gwas_y))
      for(i in 1:d){
        sdbetaGX1[,i]=gwas_X[[i]]$s.d.
      }
      betaGYC=gwas_y$beta
      pvaluebetaGYC=gwas_y$p
      sdbetaGYC=gwas_y$s.d.
      correlation1=matrix(0,nrow=1+d,ncol=1+d)
      
      for(i in 1:ncol(correlation1)){
        for(j in 1:ncol(correlation1)){
          if(i==1&j==1){
            correlation1[i,j]=1
          }
          if(i==1&j>1){
            x1=(gwas_y$`t-stat`)[gwas_y$p>0.1&gwas_X[[j-1]]$p>0.1]
            x2=(gwas_X[[j-1]]$`t-stat`)[gwas_y$p>0.1&gwas_X[[j-1]]$p>0.1]
            correlation1[i,j]=cor(x1,x2)
          }
          if(i>1&j==1){
            x1=(gwas_y$`t-stat`)[gwas_y$p>0.1&gwas_X[[i-1]]$p>0.1]
            x2=(gwas_X[[i-1]]$`t-stat`)[gwas_y$p>0.1&gwas_X[[i-1]]$p>0.1]
            correlation1[i,j]=cor(x1,x2)
          }
          if(i>1&j>1){
            x1=(gwas_X[[i-1]]$`t-stat`)[gwas_X[[j-1]]$p>0.1&gwas_X[[i-1]]$p>0.1]
            x2=(gwas_X[[j-1]]$`t-stat`)[gwas_X[[j-1]]$p>0.1&gwas_X[[i-1]]$p>0.1]
            correlation1[i,j]=cor(x1,x2)
          }
        }
      }
      
      
      
      
      
      print("summary ready")
      
      if(d==1){
        pvalX=5e-8
        IV=(pvaluebetaGX1[,1]<pvalX)
      }
      
      
      if(d==2){
        pvalX=5e-8
        IV=(pvaluebetaGX1[,1]<pvalX)|(pvaluebetaGX1[,2]<pvalX)
      }
      
      if(d==4){
        pvalX=5e-8
        IV=(pvaluebetaGX1[,1]<pvalX)|(pvaluebetaGX1[,2]<pvalX)
        for(i in 3:ncol(pvaluebetaGX1)){
          IV=IV|(pvaluebetaGX1[,i]<pvalX)
        }
      }
      
      if(d==6|d==8){
        pvalX=5e-8
        IV=(pvaluebetaGX1[,1]<pvalX )|(pvaluebetaGX1[,2]<pvalX )
        for(i in 3:ncol(pvaluebetaGX1)){
          IV=IV|(pvaluebetaGX1[,i]<pvalX )
        }
      }
      
      d=2
      sum(IV)
      rs_IV=colnames(Gmatrix)[IV]
      used=which(IV==TRUE)
      betaGX_IV=as.matrix(betaGX1[IV,])
      betaGYC_IV=betaGYC[IV]
      sdbetaGX_IV=as.matrix(sdbetaGX1[IV,])
      sdbetaGYC_IV=sdbetaGYC[IV]
      for(i in 1:nrow(betaGX_IV)){
        if(betaGX_IV[i,1]<0){
          betaGX_IV[i,]=-betaGX_IV[i,]
          betaGYC_IV[i]=-betaGYC_IV[i]
        }
      }
      
      SIG=list()
      for(i in 1:sum(IV)){
        SIG[[i]]=diag(c(sdbetaGYC_IV[i],sdbetaGX_IV[i,]))%*%correlation1%*%diag(c(sdbetaGYC_IV[i],sdbetaGX_IV[i,]))
      }
      SIGU1=list()
      SIGU2=list()
       
      
      
      
      
      
      
      print("before function")
      
      
      l=function(betaxtemp=betaGX_IV,betayctemp=betaGYC_IV,SIGtemp=SIG,betaytemp=betaGYC_IV,btemp=rep(0,times=ncol(betaGX_IV))){
        a1=numeric()
        a2=numeric()
        for(i in 1:length(SIGtemp)){
          a1[i]=t(c(betayctemp[i]-t(btemp)%*%betaxtemp[i,]-betaytemp[i],betaGX_IV[i,]-betaxtemp[i,]))%*%solve(SIGtemp[[i]])%*%c(betayctemp[i]-t(btemp)%*%betaxtemp[i,]-betaytemp[i],betaGX_IV[i,]-betaxtemp[i,])
          a2[i]=t(c(betayctemp[i]-t(btemp)%*%betaxtemp[i,]-0,betaGX_IV[i,]-betaxtemp[i,]))%*%solve(SIGtemp[[i]])%*%c(betayctemp[i]-t(btemp)%*%betaxtemp[i,]-0,betaGX_IV[i,]-betaxtemp[i,])
        }
        return(list(d=a2-a1,like=sum(a1)))
      }
      
      cml_MA=function(N=Nsample,p=ncol(betaGX_IV),m=nrow(betaGX_IV),betax=betaGX_IV,betayc=betaGYC_IV,SIG1=SIGU1 ,threshold=0.01,maxit=500){
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
      
      
      
      
      Want=cml_MA(N=Nsample,p=ncol(betaGX1),m=nrow(betaGX_IV),betax=betaGX_IV,betayc=betaGYC_IV,SIG1=SIG,threshold=0.01,maxit=500)
      
      b=Want$b
      covb=Want$Covb
      if(min(eigen(covb)$values)<0){
        covb=covb+diag((-1.1)*min(eigen(covb)$values),nrow=d,ncol = d)
      }
      
      library("MendelianRandomization")
      DP=200
      bDP_egger=matrix(0,nrow=DP,ncol=d)
     
      for(j in 1: DP){
        betaX=matrix(0,nrow=nrow(betaGX_IV),ncol = d)
        betaY=numeric()
        for(i in 1:nrow(betaGX_IV)){
          beta=mvrnorm(n=1,mu=c(betaGYC_IV[i],betaGX_IV[i,]),Sigma = SIG[[i]])
          betaX[i,]=beta[2:length(beta)]
          betaY[i]=beta[1]
        }
        inp=mr_mvinput(bx=betaX,bxse =sdbetaGX_IV,by=betaY,byse =sdbetaGYC_IV)
        egger=mr_mvegger(object = inp) 
        bDP_egger[j,]=egger$Estimate
     
        
      }
      b_egger=colMeans(bDP_egger)
      covb_egger=cov(bDP_egger)
      
      
      betaGYadj_egger=numeric()
      varbetaGYadj_egger1=numeric()
      for(j1 in 1:Nsnps){
        betaGYadj_egger[j1]=betaGYC[j1]-b_egger%*%betaGX1[j1,]
         
        covbetax=diag(x=c(sdbetaGX1[j1,]),nrow=d,ncol=d)%*%correlation1[2:(1+d),2:(1+d)]%*%diag(x=c(sdbetaGX1[j1,]),nrow=d,ncol=d)
        
        varbetaGYadj_egger1[j1]=sdbetaGYC[j1]^2+sum(diag(covbetax*covb_egger))+sum(b_egger^2*diag(covbetax))+t(betaGX1[j1,])%*%covb_egger%*%betaGX1[j1,]
      }
      
      
      
      SDbetaGYadj_egger1=sqrt(varbetaGYadj_egger1) 
      
       
      
      
      z_egger1=abs(betaGYadj_egger/SDbetaGYadj_egger1)
      
      
      
      
      pvalueadj_egger1=(1-pnorm(z_egger1))*2
       
      
      
      valid=Want$validuse
      betaxvalid=Want$betaxvalid
      betaxIVvalid=Want$betaxIVvalid
      used=which(IV==TRUE)
      used=used[valid]
      Omega=Want$Omega
      print("function finish")
     
      varbeta=numeric()
      for(i in 1:Nsnps){
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
            ZX[[i1]][i2,]=gwas_X[[i2]][used[i1],]
          }
          ZY[i1,]=gwas_y[used[i1],]
        }
        
        for(i1 in 1:d){
          GX[i1,]=gwas_X[[i1]][i,]
        }
        GY[1,]=gwas_y[i,]
        sigmaG=diag(c(GY$s.d.,GX$s.d.))%*%correlation1%*%diag(c(GY$s.d.,GX$s.d.))
        ld=numeric()
        for(i1 in 1:length(used)){
          ld[i1]=LD[i,used[i1]]
        }
        sigmaGZY=list()
        for(i1 in 1:length(used)){
          sigmaGZY[[i1]]=matrix(ld[i1]*c(GY$s.d.,GX$s.d.)*ZY$s.d.[i1]*correlation1[1,],nrow=d+1,ncol=1)
        }
        sigmaGZH=list()
        for(i1 in 1:length(used)){
          sigmaGZH[[i1]]=matrix(0,nrow=d+1,ncol=d)
          for(i2 in 1:d+1){
            if(i2==1){
              sigmaGZH[[i1]][i2,]=ld[i1]*GY$s.d.*ZX[[i1]]$s.d.*correlation1[1,(2:(d+1))]
            }
            if(i2>1){
              sigmaGZH[[i1]][i2,]=ld[i1]*GX$s.d.[i2-1]*ZX[[i1]]$s.d.*correlation1[i2,(2:(d+1))]
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
          U1=U1+(-min(eigen(U1)$values))*1.1*diag(1,nrow=nrow(U1),ncol=ncol(U1))
        }
        varbeta[i]=t(v)%*%U1%*%v
        print(i)
      }
       
      
      betaGYadj=numeric()
      for(i in 1:nsnpstotal){
        betaGYadj[i]=betaGYC[i]-t(b)%*%betaGX1[i,]
      }
      
      sdbetaGYadj3=sqrt(varbeta)
      
       
      z3=abs(betaGYadj/sdbetaGYadj3)
      
      
      
      
      
      pvalueadj3=(1-pnorm(z3))*2
      library(indexevent)
      library(SlopeHunter)
      d=1
      used=which(IV==TRUE)
      betaGX_IV=as.matrix(betaGX1[IV,][,1])
      betaGYC_IV=betaGYC[IV]
      sdbetaGX_IV=as.matrix(sdbetaGX1[IV,][,1])
      sdbetaGYC_IV=sdbetaGYC[IV]
      for(i in 1:nrow(betaGX_IV)){
        if(betaGX_IV[i,1]<0){
          betaGX_IV[i,]=-betaGX_IV[i,]
          betaGYC_IV[i]=-betaGYC_IV[i]
        }
      }
       
      for(i in 1:sum(IV)){
        SIGU1[[i]]=SIG[[i]][(1:2),(1:2)]
        SIGU2[[i]]=SIG[[i]][(-2),(-2)] 
      }
      
      Want1=cml_MA(N=Nsample,p=ncol( betaGX_IV),m=nrow(betaGX_IV),betax=betaGX_IV,betayc=betaGYC_IV,SIG1=SIGU1,threshold=0.1,maxit=500)
      
      b1=Want1$b
      covb1=Want1$Covb
      
      bDP_uv_egger1=matrix(0,nrow=DP,ncol=ncol(betaGX_IV))
      bDP_uv_DD1=matrix(0,nrow=DP,ncol=ncol(betaGX_IV))
      bDP_uv_SH1=matrix(0,nrow=DP,ncol=ncol(betaGX_IV))
      for(j in 1: DP){
        betaX=matrix(0,nrow=nrow(betaGX_IV),ncol = ncol(betaGX_IV))
        betaY=numeric()
        for(i in 1:nrow(betaGX_IV)){
          beta=mvrnorm(n=1,mu=c(betaGYC_IV[i],betaGX_IV[i,]),Sigma = SIGU1[[i]])
          betaX[i,]=beta[2:length(beta)]
          betaY[i]=beta[1]
        }
        inp=mr_mvinput(bx=betaX,bxse =sdbetaGX_IV,by=betaY,byse =sdbetaGYC_IV)
        egger=mr_mvegger(object = inp) 
        mr_mvivw(object = inp)
        
        DD=indexevent(seed = 2018,lambda = seq(0.25, 5, 0.25),B = 10, method = c("Hedges-Olkin"),prune = NULL,xbeta=betaX, xse=sdbetaGX_IV,ybeta=betaY,yse=sdbetaGYC_IV,weighted = T)
        hunter=data.frame(xbeta=betaX, xse=sdbetaGX_IV, ybeta=betaY,  yse=sdbetaGYC_IV,rs=rs_IV  )
        Hunter=hunt(dat=hunter,snp_col = "rs",xbeta_col = "xbeta",xse_col = "xse", ybeta_col = "ybeta",yse_col = "yse",xp_thresh = 0.001,init_pi = 0.6,init_sigmaIP = 1e-05,Bootstrapping = FALSE,M = 200,seed = 777,Plot = FALSE,show_adjustments = FALSE)
        bDP_uv_egger1[j,]=egger$Estimate
        bDP_uv_DD1[j,]=DD$b
        bDP_uv_SH1[j,]=Hunter$b
        }
      b_uv_egger1=colMeans(bDP_uv_egger1)
      covb_uv_egger1=cov(bDP_uv_egger1)
      b_uv_DD1=colMeans(bDP_uv_DD1)
      covb_uv_DD1=cov(bDP_uv_DD1)
      b_uv_SH1=colMeans(bDP_uv_SH1)
      covb_uv_SH1=cov(bDP_uv_SH1)
      
      
      used=which(IV==TRUE)
      betaGX_IV=as.matrix(betaGX1[IV,][,2])
      betaGYC_IV=betaGYC[IV]
      sdbetaGX_IV=as.matrix(sdbetaGX1[IV,][,2])
      sdbetaGYC_IV=sdbetaGYC[IV]
      for(i in 1:nrow(betaGX_IV)){
        if(betaGX_IV[i,1]<0){
          betaGX_IV[i,]=-betaGX_IV[i,]
          betaGYC_IV[i]=-betaGYC_IV[i]
        }
      }
      
      
      
      
      Want2=cml_MA(N=Nsample,p=ncol( betaGX_IV),m=nrow(betaGX_IV),betax=betaGX_IV,betayc=betaGYC_IV,SIG1=SIGU1,threshold=0.001,maxit=500)
      
      b2=Want2$b
      covb2=Want2$Covb
      
    
      
      bDP_uv_egger2=matrix(0,nrow=DP,ncol=ncol(betaGX_IV))
      bDP_uv_DD2=matrix(0,nrow=DP,ncol=ncol(betaGX_IV))
      bDP_uv_SH2=matrix(0,nrow=DP,ncol=ncol(betaGX_IV))
      for(j in 1: DP){
        betaX=matrix(0,nrow=nrow(betaGX_IV),ncol = ncol(betaGX_IV))
        betaY=numeric()
        for(i in 1:nrow(betaGX_IV)){
          beta=mvrnorm(n=1,mu=c(betaGYC_IV[i],betaGX_IV[i,]),Sigma = SIGU2[[i]])
          betaX[i,]=beta[2:length(beta)]
          betaY[i]=beta[1]
        }
        
        inp=mr_mvinput(bx=betaX,bxse =sdbetaGX_IV,by=betaY,byse =sdbetaGYC_IV)
        egger=mr_mvegger(object = inp) 
        bDP_uv_egger2[j,]=egger$Estimate
        bDP_uv_DD2[j,]=DD$b
        bDP_uv_SH2[j,]=Hunter$b
        
      }
      b_uv_egger2=colMeans(bDP_uv_egger1)
      covb_uv_egger2=cov(bDP_uv_egger1)
      b_uv_DD2=colMeans(bDP_uv_DD1)
      covb_uv_DD2=cov(bDP_uv_DD1)
      b_uv_SH2=colMeans(bDP_uv_SH1)
      covb_uv_SH2=cov(bDP_uv_SH1)
      
      
      b_uv=c(b1,b2)
      covb_uv=diag(c(covb1,covb2))
      b_uv_egger=c(b_uv_egger1,b_uv_egger2)
      covb_uv_egger=diag(c(covb_uv_egger1,covb_uv_egger2))
      b_uv_DD=c(b_uv_DD1,b_uv_DD2)
      covb_uv_DD=diag(c(covb_uv_DD1,covb_uv_DD2))
      b_uv_SH=c(b_uv_SH1,b_uv_SH2)
      covb_uv_SH=diag(c(covb_uv_SH1,covb_uv_SH2))
      
      
      
      
      betaGYadj_uv=numeric()
      varbetaGYadj_uv=numeric()
      betaGYadj_uv_egger=numeric()
      varbetaGYadj_uv_egger=numeric()
      betaGYadj_uv_DD=numeric()
      varbetaGYadj_uv_DD=numeric()
      betaGYadj_uv_SH=numeric()
      varbetaGYadj_uv_SH=numeric()
      
      
      
      
      
      
      for(j1 in 1:Nsnps){
        betaGYadj_uv[j1]=betaGYC[j1]-b_uv%*%betaGX1[j1,]
        betaGYadj_uv_egger[j1]=betaGYC[j1]-b_uv_egger%*%betaGX1[j1,]
        betaGYadj_uv_DD[j1]=betaGYC[j1]-b_uv_DD%*%betaGX1[j1,]
        betaGYadj_uv_SH[j1]=betaGYC[j1]-b_uv_SH%*%betaGX1[j1,]
        
        covbetax=diag(x=c(sdbetaGX1[j1,]),nrow=ncol(sdbetaGX1),ncol=ncol(sdbetaGX1))%*%correlation1[2:(1+ncol(sdbetaGX1)),2:(1+ncol(sdbetaGX1))]%*%diag(x=c(sdbetaGX1[j1,]),nrow=ncol(sdbetaGX1),ncol=ncol(sdbetaGX1))
        varbetaGYadj_uv[j1]=sdbetaGYC[j1]^2+sum(diag(covbetax*covb_uv))+sum(b_uv^2*diag(covbetax))+t(betaGX1[j1,])%*%covb_uv%*%betaGX1[j1,]
        varbetaGYadj_uv_egger[j1]=sdbetaGYC[j1]^2+sum(diag(covbetax*covb_uv_egger))+sum(b_uv_egger^2*diag(covbetax))+t(betaGX1[j1,])%*%covb_uv_egger%*%betaGX1[j1,]
        varbetaGYadj_uv_DD[j1]=sdbetaGYC[j1]^2+sum(diag(covbetax*covb_uv_DD))+sum(b_uv_DD^2*diag(covbetax))+t(betaGX1[j1,])%*%covb_uv_DD%*%betaGX1[j1,]
        varbetaGYadj_uv_SH[j1]=sdbetaGYC[j1]^2+sum(diag(covbetax*covb_uv_SH))+sum(b_uv_SH^2*diag(covbetax))+t(betaGX1[j1,])%*%covb_uv_SH%*%betaGX1[j1,]
      }
      sdbetaGYadj_uv=sqrt(varbetaGYadj_uv)
      sdbetaGYadj_uv_egger=sqrt(varbetaGYadj_uv_egger)
      sdbetaGYadj_uv_DD=sqrt(varbetaGYadj_uv_DD)
      sdbetaGYadj_uv_SH=sqrt(varbetaGYadj_uv_SH)
      
      
      
      z_uv=abs(betaGYadj_uv/sdbetaGYadj_uv)
      pvalueadj_uv=(1-pnorm(z_uv))*2
      z_uv_egger=abs(betaGYadj_uv_egger/sdbetaGYadj_uv_egger)
      pvalueadj_uv_egger=(1-pnorm(z_uv_egger))*2
      z_uv_DD=abs(betaGYadj_uv_DD/sdbetaGYadj_uv_DD)
      pvalueadj_uv_DD=(1-pnorm(z_uv_DD))*2
      z_uv_SH=abs(betaGYadj_uv_SH/sdbetaGYadj_uv_SH)
      pvalueadj_uv_SH=(1-pnorm(z_uv_SH))*2
      
      
      
      varG_true=2*MAF*(1-MAF)
      M22total=list()
      covG=cov(Gmatrix)
      for(i in 1:Nsnps){
        covGtemp=covG[-i,-i]
        betaGXtemp=betaGX[,-i]
        M22total[[i]]=betaUX%*%t(betaUX )*var(U )+diag(varEX)+betaGXtemp%*%covGtemp%*%t(betaGXtemp)
        print(i)
      }
      b_true=matrix(0,nrow = Nsnps,ncol = length(betaUX))
      for(i in 1:Nsnps){
        b_true[i,]=-betaUY[1,1]*solve(M22total[[i]])%*%betaUX*var(U )
      }
      
      
      
      
      
      
      
      slope=matrix(0,ncol=7*ncol(b_true),nrow=Nsnps)
      slope_se=matrix(0,ncol=7*ncol(b_true),nrow=Nsnps)
      for(i in 1:nrow(slope)){
        slope[i,]=c(b_true[i,],b,b_egger,b_uv,b_uv_egger,b_uv_DD,b_uv_SH)
        slope_se[i,]=c(rep(0,times=ncol(b_true)),sqrt(diag(covb)),sqrt(diag(covb_egger)),sqrt(diag(covb_uv)),sqrt(diag(covb_uv_egger)),sqrt(diag(covb_uv_DD)),sqrt(diag(covb_uv_SH)))
      }
      
      
      name1=numeric()
      name2=numeric()
      name3=numeric()
      name4=numeric()
      name5=numeric()
      name6=numeric()
      name7=numeric()
      
      for(i in 1:ncol(b_true)){
        name1[i]=paste("true_slop",i,sep="")
        name2[i]=paste("cml_slop",i,sep="")
        name3[i]=paste("egger_slop",i,sep="")
        name4[i]=paste("cml_uv_slop",i,sep="")
        name5[i]=paste("egger_uv_slop",i,sep="")
        name6[i]=paste("DD_uv_slop",i,sep="")
        name7[i]=paste("SH_uv_slop",i,sep="")
      }
      name1se=numeric()
      name2se=numeric()
      name3se=numeric()
      name4se=numeric()
      name5se=numeric()
      name6se=numeric()
      name7se=numeric()
      for(i in 1:ncol(b_true)){
        name1se[i]=paste("true_slop_se",i,sep="")
        name2se[i]=paste("cml_slop_se",i,sep="")
        name3se[i]=paste("egger_slop_se",i,sep="")
        name4se[i]=paste("cml_uv_slop_se",i,sep="")
        name5se[i]=paste("egger_uv_slop_se",i,sep="")
        name6se[i]=paste("DD_uv_slop_se",i,sep="")
        name7se[i]=paste("SH_uv_slop_se",i,sep="")
      }
      slope=as.data.frame(slope)
      colnames(slope)=c(name1,name2,name3,name4,name5 ,name6,name7)
      slope_se=as.data.frame(slope_se)
      colnames(slope_se)=c(name1se,name2se,name3se ,name4se,name5se ,name6se,name7se)
      
      
      result1=data.frame(rs=gwas_y$rs,beta_before=gwas_y$beta,se_before=gwas_y$s.d.,p_before=gwas_y$p,
                         beta_after=betaGYadj, se_after3=sdbetaGYadj3, p_after3=pvalueadj3,
                         beta_egger=betaGYadj_egger, se_egger=SDbetaGYadj_egger1, p_egger=pvalueadj_egger1,
                         beta_uv_cml=betaGYadj_uv,sd_uv_cml=sdbetaGYadj_uv,p_uv_cml=pvalueadj_uv, 
                         beta_uv_egger=betaGYadj_uv_egger,sd_uv_egger=sdbetaGYadj_uv_egger,p_uv_egger=pvalueadj_uv_egger,
                         beta_uv_DD=betaGYadj_uv_DD,sd_uv_DD=sdbetaGYadj_uv_DD,p_uv_DD=pvalueadj_uv_DD, 
                         beta_uv_SH=betaGYadj_uv_SH,sd_uv_SH=sdbetaGYadj_uv_SH,p_uv_SH=pvalueadj_uv_SH, 
                         beta_true=SNP_effect  )
      result=cbind(result1,slope,slope_se) 
      
    }
    write.table(output,file=outnames[t,t1],sep="\t",row.names=FALSE,quote = FALSE)  
  }
}