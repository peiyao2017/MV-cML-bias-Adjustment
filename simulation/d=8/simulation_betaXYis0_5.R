setwd("/panfs/jay/groups/20/panwei/wan01299/collider_bias/BMI/MV_Egger/d=8/")

library(doParallel)
number_of_cores=detectCores()
myCluster <- makeCluster(number_of_cores-10, # number of cores to use
                         type = "PSOCK") # type of cluster
#notice that we need to leave one core for computer system
registerDoParallel(myCluster)
library(mvtnorm)
library(MASS)
d=8
rhos=c(-0.9)
name_effects=numeric()
Nsnps=1000
N=20000

for(i in 1:length(rhos)){
  name_effects[i]=paste("effects",rhos[i],"dim",d,".RData",sep="")
}
outnames=matrix(0,nrow=length(rhos),ncol = 20)
for(t in 1:length(rhos)){
  for(t1 in 1:20){
    outnames[t,t1]=paste("sim",rhos[t],"_",t1,"betaXYis0",".txt",sep="")
  }
}

for(t in 1:length(rhos)){
  load(name_effects[t])
  betaGY=effect[[1]]
  betaGX=effect[[2]]
  MAF=effect[[3]]
  betaXY=rep(0,times=d)
  SNP_effect=betaGY
  varG=MAF*(1-MAF)
  for(t1 in 1:20){
    output=foreach(j=1:50,.combine = "rbind")%dopar%{
      library(mvtnorm)
      library(MASS)
      
      repeat{
        X=matrix(0,nrow=N,ncol=d)
        betaUX=numeric()
        for(i in 1:d){
          betaUX[i]=0.8*sum(betaGX[i,]^2*varG)
        }
        varEX=numeric()
        for(i in 1:d){
          varEX[i]=0.2*sum(betaGX[i,]^2*varG)
        }
        Gmatrix=matrix(0,nrow=N,ncol=Nsnps)
        for(i in 1:ncol(Gmatrix)){
          for(j in 1:nrow(Gmatrix)){
            Gmatrix[j,i]=rbinom(n=1,size = 2,prob = MAF[i])
          }
        }
        for(i in 1:ncol(Gmatrix)){
          Gmatrix[,i]=scale(Gmatrix[,i],center = TRUE,scale=FALSE)
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
        betaUY=0.8*var(y_known)
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
              x1=gwas_y$`t-stat`[gwas_y$p>5e-1&gwas_X[[j-1]]$p>5e-1]
              x2=gwas_X[[j-1]]$`t-stat`[gwas_y$p>5e-1&gwas_X[[j-1]]$p>5e-1]
              correlation1[i,j]=cor(x1,x2)
            }
            if(i>1&j==1){
              x1=gwas_y$`t-stat`[gwas_y$p>5e-1&gwas_X[[i-1]]$p>5e-1]
              x2=gwas_X[[i-1]]$`t-stat`[gwas_y$p>5e-1&gwas_X[[i-1]]$p>5e-1]
              correlation1[i,j]=cor(x1,x2)
            }
            if(i>1&j>1){
              x1=gwas_X[[j-1]]$`t-stat`[gwas_X[[i-1]]$p>5e-1&gwas_X[[j-1]]$p>5e-1]
              x2=gwas_X[[i-1]]$`t-stat`[gwas_X[[i-1]]$p>5e-1&gwas_X[[j-1]]$p>5e-1]
              correlation1[i,j]=cor(x1,x2)
            }
          }
        }
        print("summary ready")
        if(d==2){
          pvalX=5e-2
          IV=(pvaluebetaGX1[,1]<pvalX)|(pvaluebetaGX1[,2]<pvalX)
        }
        
        if(d==4){
          pvalX=5e-2
          IV=(pvaluebetaGX1[,1]<pvalX)|(pvaluebetaGX1[,2]<pvalX)
          for(i in 3:ncol(pvaluebetaGX1)){
            IV=IV|(pvaluebetaGX1[,i]<pvalX)
          }
        }
        
        if(d==6|d==8){
          pvalX=5e-2
          IV=(pvaluebetaGX1[,1]<pvalX )|(pvaluebetaGX1[,2]<pvalX )
          for(i in 3:ncol(pvaluebetaGX1)){
            IV=IV|(pvaluebetaGX1[,i]<pvalX )
          }
        }
        
        
        sum(IV)
        
        betaGX_IV=as.matrix(betaGX1[IV,])
        betaGYC_IV=betaGYC[IV]
        sdbetaGX_IV=as.matrix(sdbetaGX1[IV,])
        sdbetaGYC_IV=sdbetaGYC[IV]
        if(betaGX_IV[1,1]>=0){
          betaGX_IV=betaGX_IV
          betaGYC_IV=betaGYC_IV
        }
        if(betaGX_IV[1,1]<0){
          betaGX_IV=-betaGX_IV
          betaGYC_IV=-betaGYC_IV
        }
        SIG=list()
        for(i in 1:sum(IV)){
          SIG[[i]]=diag(c(sdbetaGYC_IV[i],sdbetaGX_IV[i,]))%*%correlation1%*%diag(c(sdbetaGYC_IV[i],sdbetaGX_IV[i,]))
        }
        
        print("before function")
        library("MendelianRandomization")
        DP=200
        bDP=matrix(0,nrow=DP,ncol=d)
        for(j in 1: DP){
          betaX=matrix(0,nrow=nrow(betaGX_IV),ncol = d)
          betaY=numeric()
          for(i in 1:nrow(betaGX_IV)){
            beta=mvrnorm(n=1,mu=c(betaGYC[i],betaGX_IV[i,]),Sigma = SIG[[i]])
            betaX[i,]=beta[2:length(beta)]
            betaY[i]=beta[1]
          }
          inp=mr_mvinput(bx=betaX,bxse =sdbetaGX_IV,by=betaY,byse =sdbetaGYC_IV)
          bDP[j,]=mr_mvegger(object = inp)$Estimate
        }
        b_egger=colMeans(bDP)
        covb_egger=cov(bDP)/DP
        
        
        print("function finish")
        
        betaGYadj_egger=numeric()
        for(i in 1:nsnpstotal){
          betaGYadj_egger[i]=betaGYC[i]-t(b_egger)%*%betaGX1[i,]
        }
        
        sdbetaGYadj_egger=numeric()
        for(i in 1:nsnpstotal){
          a=sdbetaGYC[i]^2
          b1=t(betaGX1[i,])%*%covb_egger%*%betaGX1[i,]
          c=sum(b_egger^2*as.vector(sdbetaGX1[i,]^2))
          d1=sum(diag(covb_egger)*sdbetaGX1[i,]^2)
          sdbetaGYadj_egger[i]=sqrt(a+b1+c+d1)
        }
        
        z_egger=abs(betaGYadj_egger/sdbetaGYadj_egger)
        pvalueadj_egger=(1-pnorm(z_egger))*2
        
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
          Info1=rbind(A,B1)
          Info2=rbind(t(B1),C1)
          Info3=cbind(Info1,Info2)
          Info=solve(Info3)
          Covb=Info[(1:length(b)),(1:length(b))]
          return(list(b=b,Covb=Covb))
        }  
        
        
        
        
        Want=cml_MA(N=Nsample,p=ncol(betaGX1),m=nrow(betaGX_IV),betax=betaGX_IV,betayc=betaGYC_IV,SIG1=SIG,threshold=0.04,maxit=500)
        
        b=Want$b
        covb=Want$Covb
        print("function finish")
        
        betaGYadj=numeric()
        for(i in 1:nsnpstotal){
          betaGYadj[i]=betaGYC[i]-t(b)%*%betaGX1[i,]
        }
        
        sdbetaGYadj=numeric()
        for(i in 1:nsnpstotal){
          a=sdbetaGYC[i]^2
          b1=t(betaGX1[i,])%*%covb%*%betaGX1[i,]
          c=sum(b^2*as.vector(sdbetaGX1[i,]^2))
          d1=sum(diag(covb)*sdbetaGX1[i,]^2)
          sdbetaGYadj[i]=sqrt(a+b1+c+d1)
        }
        
        z=abs(betaGYadj/sdbetaGYadj)
        pvalueadj=(1-pnorm(z))*2
        
        
        
        if(min(eigen(covb)$values)>=0){
          break
        }
      }
      bias=numeric()
      bias_true=numeric()
      M22=solve(betaUX%*%t(betaUX)+diag(varEX))
      b_true=betaUY[1,1]*M22%*%betaUX
      for(i in 1:nrow(gwas_y)){
        bias[i]=t(b)%*%betaGX1[i,]
        bias_true[i]= t(b_true)%*%betaGX[,i]
      }
      result=data.frame(rs=gwas_y$rs,beta_before=gwas_y$beta,sd_before=gwas_y$s.d.,p_before=gwas_y$p,beta_after=betaGYadj,sd_after=sdbetaGYadj,p_after=pvalueadj,beta_egger=betaGYadj_egger,sd_egger=sdbetaGYadj_egger,p_egger=pvalueadj_egger,beta_true=SNP_effect,bias=bias,bias_true=bias_true)
    }
    write.table(output,file=outnames[t,t1],sep="\t",row.names=FALSE,quote = FALSE)  
  }
}
