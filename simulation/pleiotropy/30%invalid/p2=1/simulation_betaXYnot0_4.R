setwd("/panfs/jay/groups/20/panwei/wan01299/collider_bias/BMI/MV_Egger/new/real_snp_sim/pleiotropy/30invalid/d=1/")

library(doParallel)
number_of_cores=detectCores()
myCluster <- makeCluster(number_of_cores-12, # number of cores to use
                         type = "PSOCK") # type of cluster
#notice that we need to leave one core for computer system
registerDoParallel(myCluster)
library(mvtnorm)
library(MASS)
d=1
rhos=c(-0.5)
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
  betaGY=c(effect[[1]],rep(0,3000))
  betaGX=cbind(effect[[2]],matrix(0,nrow=d,ncol = 3000))
 
  betaXY=effect[[4]]
  SNP_effect=betaGY
  valid_percent=0.7
  nIV=30
  for(t1 in 1:20){
    
    output=foreach(jt=1:50,.combine = "rbind")%dopar%{
      library(mvtnorm)
      library(MASS)
      library(BEDMatrix)
     
       load("Geno30.RDATA")
      Gmatrix=Geno[[1]][sample(c(1:nrow(Geno[[1]])),size=Nsample,replace = FALSE),]
      simID=Geno[[2]]
      simIV=Geno[[3]]
       rm(Geno)

      
      LD=cor(Gmatrix[1:1000,1:1000])
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
      EX=mvrnorm(n=N,mu=rep(0,times=d),Sigma =  (varEX))
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
      snp_names=simID
      
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
      rs=colnames(G_sim)
      
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
            x1=(gwas_y$`t-stat`[1001:4000])[gwas_y$p[1001:4000]>0.1&gwas_X[[j-1]]$p[1001:4000]>0.1]
            x2=(gwas_X[[j-1]]$`t-stat`[1001:4000])[gwas_y$p[1001:4000]>0.1&gwas_X[[j-1]]$p[1001:4000]>0.1]
            correlation1[i,j]=cor(x1,x2)
          }
          if(i>1&j==1){
            x1=(gwas_y$`t-stat`[1001:4000])[gwas_y$p[1001:4000]>0.1&gwas_X[[i-1]]$p[1001:4000]>0.1]
            x2=(gwas_X[[i-1]]$`t-stat`[1001:4000])[gwas_y$p[1001:4000]>0.1&gwas_X[[i-1]]$p[1001:4000]>0.1]
            correlation1[i,j]=cor(x1,x2)
          }
          if(i>1&j>1){
            x1=(gwas_X[[i-1]]$`t-stat`[1001:4000])[gwas_X[[i-1]]$p[1001:4000]>0.1&gwas_X[[j-1]]$p[1001:4000]>0.1]
            x2=(gwas_X[[j-1]]$`t-stat`[1001:4000])[gwas_X[[i-1]]$p[1001:4000]>0.1&gwas_X[[j-1]]$p[1001:4000]>0.1]
            correlation1[i,j]=cor(x1,x2)
          }
        }
      }
      print("summary ready")
      
      IV=simID%in%simIV
      
      sum(IV)
      
      betaGX_IV=as.matrix(betaGX1[IV,])
      betaGYC_IV=betaGYC[IV]
      sdbetaGX_IV=as.matrix(sdbetaGX1[IV,])
      sdbetaGYC_IV=sdbetaGYC[IV]
      pvaluebetaGX_IV=pvaluebetaGX1[IV,]
      pvaluebetaGYC_IV=pvaluebetaGYC[IV]
      rs_IV=rs[IV]
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
      
      print("before function")
      library("MendelianRandomization")
      library(SlopeHunter)
      library(indexevent)
      DP=200
      bDP_egger=matrix(0,nrow=DP,ncol=d)
      bDP_ivw=matrix(0,nrow=DP,ncol=d)
      bDP_lasso=matrix(0,nrow=DP,ncol=d)
      bDP_median=matrix(0,nrow=DP,ncol=d)
      bDP_cml=matrix(0,nrow=DP,ncol=d)
      bDP_DD=matrix(0,nrow=DP,ncol=d)
      bDP_SH=matrix(0,nrow=DP,ncol=d)
      IV_DP=list()
      betaxvalid_DP=list()
      betaxIVvalid_DP=list()
      for(j in 1: DP){
        betaX=matrix(0,nrow=nrow(betaGX_IV),ncol = d)
        betaY=numeric()
        for(i in 1:nrow(betaGX_IV)){
          set.seed(NULL)
          beta=mvrnorm(n=1,mu=c(betaGYC_IV[i],betaGX_IV[i,]),Sigma = SIG[[i]])
          betaX[i,]=beta[2:length(beta)]
          betaY[i]=beta[1]
        }
        inp=mr_mvinput(bx=betaX,bxse =sdbetaGX_IV,by=betaY,byse =sdbetaGYC_IV)
        inp1=mr_input(bx=as.vector(betaX),bxse =as.vector(sdbetaGX_IV),by=betaY,byse =sdbetaGYC_IV)
        egger=mr_mvegger(object = inp) 
        ivw=mr_mvivw(object = inp)
        lasso=mr_mvlasso(object = inp,lambda=seq(from=0,to=4,by=0.01))
        median1=mr_median(object = inp1,iterations=100)
        
        hunter=data.frame(xbeta=betaX, xse=sdbetaGX_IV, ybeta=betaY,  yse=sdbetaGYC_IV,rs=rs_IV  )
        Hunter=hunt(dat=hunter,snp_col = "rs",xbeta_col = "xbeta",xse_col = "xse", ybeta_col = "ybeta",yse_col = "yse",xp_thresh = 0.001,init_pi = 0.6,init_sigmaIP = 1e-05,Bootstrapping = FALSE,M = 200,seed = 777,Plot = FALSE,show_adjustments = FALSE)
        DD=indexevent(seed = 2018,lambda = seq(0.25, 5, 0.25),B = 10, method = c("Hedges-Olkin"),prune = NULL,xbeta=betaX, xse=sdbetaGX_IV,ybeta=betaY,yse=sdbetaGYC_IV,weighted = T)
        bDP_egger[j,]=egger$Estimate
        bDP_ivw[j,]=ivw$Estimate
        bDP_lasso[j,]=lasso$Estimate
        bDP_median[j,]=median1$Estimate
        bDP_DD[j,]=DD$b
        bDP_SH[j,]=Hunter$b
      }
      b_egger=colMeans(bDP_egger)
      covb_egger=cov(bDP_egger)
      b_ivw=colMeans(bDP_ivw)
      covb_ivw=cov(bDP_ivw)
      b_lasso=colMeans(bDP_lasso)
      covb_lasso=cov(bDP_lasso) 
      b_median=colMeans(bDP_median)
      covb_median=cov(bDP_median) 
      b_DD=colMeans(bDP_DD)
      covb_DD=cov(bDP_DD) 
      b_SH=colMeans(bDP_SH)
      covb_SH=cov(bDP_SH) 
      
      
      
      
      
      
      
      
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
      
      
      for(j1 in 1:Nsnps){
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
      
      
      
      
      
      Want=cml_MA(N=Nsample,p=ncol(betaGX1),m=nrow(betaGX_IV),betax=betaGX_IV,betayc=betaGYC_IV,SIG1=SIG,threshold=0.001,maxit=500)
      
      b=Want$b
      covb=Want$Covb
      if(min(eigen(covb)$values)<0){
        covb=covb+diag((-1.1)*min(eigen(covb)$values),nrow=d,ncol = d)
      }
      valid=Want$validuse
      betaxvalid=Want$betaxvalid
      betaxIVvalid=Want$betaxIVvalid
      used=which(IV==TRUE)
      used=used[valid]
      Omega=Want$Omega
      print("function finish")
      sdbetaGYadj1=numeric()
      sdbetaGYadj2=numeric()
      varbeta=numeric()
      for(i in 1:1000){
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
      for(i in 1:Nsnps){
        covbetax=diag(sdbetaGX1[i,],nrow=d,ncol=d)%*%correlation1[2:nrow(correlation1),2:nrow(correlation1)]%*%diag(sdbetaGX1[i,],nrow=d,ncol=d)
        covyx=sdbetaGYC[i]*correlation1[1,2:ncol(correlation1)]*sdbetaGX1[i,]
        a=sdbetaGYC[i]^2
        b1=t(betaGX1[i,])%*%covb%*%betaGX1[i,]
        c1=t(b)%*%covbetax%*%b
        c2=t(b)%*%covbetax%*%b
        sdbetaGYadj1[i]=sqrt(a+b1+c1+sum(covb*covbetax))
        sdbetaGYadj2[i]=sqrt(a+b1+c2+sum(covb*covbetax)-2*sum(covyx*b))
      }
      
      
      betaGYadj=numeric()
      for(i in 1:nsnpstotal){
        betaGYadj[i]=betaGYC[i]-t(b)%*%betaGX1[i,]
      }
      
      sdbetaGYadj3=sqrt(varbeta)
      
      
      z1=abs(betaGYadj/sdbetaGYadj1)
      z2=abs(betaGYadj/sdbetaGYadj2)
      z3=abs(betaGYadj/sdbetaGYadj3)
      
      
      
      
      pvalueadj1=(1-pnorm(z1))*2
      pvalueadj2=(1-pnorm(z2))*2
      pvalueadj3=(1-pnorm(z3))*2
      
      
      
      
      
      
       
      
      M22total=list()
      covG=cov(Gmatrix[1:1000,1:1000])
      for(i in 1:1000){
        covGtemp=covG[-i,-i]
        betaGXtemp=matrix(betaGX[,1:1000],nrow=1,ncol=1000)[,-i]
        M22total[[i]]=betaUX%*%t(betaUX)*var(U)+varEX+t(betaGXtemp)%*%covGtemp%*%betaGXtemp
        print(i)
      }
      b_true=matrix(0,nrow = 1000,ncol = d)
      for(i in 1:Nsnps){
        b_true[i,]=-betaUY[1,1]*solve(M22total[[i]])%*%betaUX*var(U)
      }
      slope=matrix(0,ncol=8*d,nrow=Nsnps)
      slope_se=matrix(0,ncol=8*d,nrow=Nsnps)
      for(i in 1:nrow(slope)){
        slope[i,]=c(b_true[i,],b,b_egger,b_DD,b_SH,b_ivw,b_lasso,b_median)
        slope_se[i,]=c(0,sqrt(covb),sqrt(covb_egger),sqrt(covb_DD),sqrt(covb_SH),sqrt(covb_ivw),sqrt(covb_lasso),sqrt(covb_median))
      }
      
      
      name1=numeric()
      name2=numeric()
      name3=numeric()
      name4=numeric()
      name5=numeric()
      name6=numeric()
      name7=numeric()
      name8=numeric()
      for(i in 1:d){
        name1[i]=paste("true_slop",i,sep="")
        name2[i]=paste("cml_slop",i,sep="")
        name3[i]=paste("egger_slop",i,sep="")
        name4[i]=paste("DD_slop",i,sep="")
        name5[i]=paste("SH_slop",i,sep="")
        name6[i]=paste("ivw_slop",i,sep="")
        name7[i]=paste("lasso_slop",i,sep="")
        name8[i]=paste("median_slop",i,sep="")
      }
      name1se=numeric()
      name2se=numeric()
      name3se=numeric()
      name4se=numeric()
      name5se=numeric()
      name6se=numeric()
      name7se=numeric()
      name8se=numeric()
      for(i in 1:d){
        name1se[i]=paste("true_slop_se",i,sep="")
        name2se[i]=paste("cml_slop_se",i,sep="")
        name3se[i]=paste("egger_slop_se",i,sep="")
        name4se[i]=paste("DD_slop_se",i,sep="")
        name5se[i]=paste("SH_slop_se",i,sep="")
        name6se[i]=paste("ivw_slop_se",i,sep="")
        name7se[i]=paste("lasso_slop_se",i,sep="")
        name8se[i]=paste("median_slop_se",i,sep="")
      }
      slope=as.data.frame(slope)
      colnames(slope)=c(name1,name2,name3,name4,name5,name6,name7,name8)
      slope_se=as.data.frame(slope_se)
      colnames(slope_se)=c(name1se,name2se,name3se,name4se,name5se,name6se,name7se,name8se)
      
    
       result1=data.frame(rs=gwas_y$rs[1:1000],beta_before=gwas_y$beta[1:1000],se_before=gwas_y$s.d.[1:1000],
                         p_before=gwas_y$p[1:1000],beta_after=betaGYadj[1:1000],se_after1=sdbetaGYadj1[1:1000],
                         se_after2=sdbetaGYadj2[1:1000],se_after3=sdbetaGYadj3[1:1000],p_after1=pvalueadj1[1:1000],
                         p_after2=pvalueadj2[1:1000],p_after3=pvalueadj3[1:1000],beta_DD=betaGYadj_DD[1:1000],
                         se_DD1=SDbetaGYadj_DD1[1:1000],se_DD2=SDbetaGYadj_DD2[1:1000],p_DD1=pvalueadj_DD1[1:1000],
                         p_DD2=pvalueadj_DD2[1:1000],beta_egger=betaGYadj_egger[1:1000],se_egger1=SDbetaGYadj_egger1[1:1000],
                         se_egger2=SDbetaGYadj_egger2[1:1000],p_egger1=pvalueadj_egger1[1:1000],p_egger2=pvalueadj_egger2[1:1000],
                         beta_SH=betaGYadj_SH[1:1000],se_SH1=SDbetaGYadj_SH1[1:1000],se_SH2=SDbetaGYadj_SH2[1:1000],
                         p_SH1=pvalueadj_SH1[1:1000],p_SH2=pvalueadj_SH2[1:1000],beta_median=betaGYadj_median[1:1000],
                         se_median1=SDbetaGYadj_median1[1:1000],se_median2=SDbetaGYadj_median2[1:1000],p_median1=pvalueadj_median1[1:1000],
                         p_median2=pvalueadj_median2[1:1000],beta_ivw=betaGYadj_ivw[1:1000],se_ivw1=SDbetaGYadj_ivw1[1:1000],
                         se_ivw2=SDbetaGYadj_ivw2[1:1000],p_ivw1=pvalueadj_ivw1[1:1000],p_ivw2=pvalueadj_ivw2[1:1000],beta_lasso=betaGYadj_lasso[1:1000],
                         se_lasso1=SDbetaGYadj_lasso1[1:1000],se_lasso2=SDbetaGYadj_lasso2[1:1000],p_lasso1=pvalueadj_lasso1[1:1000],
                         p_lasso2=pvalueadj_lasso2[1:1000],beta_true=SNP_effect[1:1000]) 
      result=cbind(result1,slope,slope_se) 
    }
    write.table(output,file=outnames[t,t1],sep="\t",row.names=FALSE,quote = FALSE)  
  }
}

