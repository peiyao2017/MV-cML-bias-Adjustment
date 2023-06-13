


N=20
name=numeric()
betaXY=0

for(d in c(1,2,4,6,8)){
  for(rho in c(0,0.5,0.9,-0.5,-0.9)){
    if(10>0){
      filename=paste("D:/art sim/result/d=",d,"/d=",d,"betaXYis0rho",rho,".RData",sep="")
      setwd(paste("D:/art sim/","d=",d,sep=""))
for(i in 1:N){
  name[i]=paste("sim",rho,"_",i,"betaXYis",betaXY,".txt",sep="")
}
result=list()
for(i in 1:N){
  result[[i]]=read.table(file=name[i],header = TRUE,sep="\t")
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


all_snps[[1]]=data.frame(rs=rep(0,times=1000),beta_before=rep(0,times=1000),sd_before=rep(0,times=1000),p_before=rep(0,times=1000),beta_after=rep(0,times=1000),sd_after=rep(0,times=1000),p_after=rep(0,times=1000),beta_egger=rep(0,times=1000),sd_egger=rep(0,times=1000),p_egger=rep(0,times=1000),beta_true=rep(0,times=1000),bias=rep(0,times=1000),bias_true=rep(0,times=1000))

all_snps[[1]]=cbind(all_snps[[1]],matrix(0,ncol=3*d,nrow=nrow(all_snps[[1]])))



  for(j in  1:1000){
    all_snps[[1]][j,]=a[1+(j-1)*1000,]
  }
  data=all_snps[[1]][,(ncol(all_snps[[1]])-3*d+1):ncol(all_snps[[1]])]
  
 d_mean=colMeans(data)
 d_sd=numeric()
 for(i in 1:(3*d)){
   d_sd[i]=sd(data[,i])/sqrt(1000)
 }

 d_true=-d_mean[1:d]
 d_cml=d_mean[(d+1):(2*d)]
 d_egger=d_mean[(2*d+1):(3*d)]
 
 sd_true= d_sd[1:d]
 sd_cml= d_sd[(d+1):(2*d)]
 sd_egger= d_sd[(2*d+1):(3*d)]
 result=list(data.frame(d_true,sd_true,d_cml,sd_cml,d_egger,sd_egger))
 save(result,file=filename)
    }
  }
}


for(d in c(1,2,4,6,8)){
  for(rho in c(0,0.5,0.9,-0.5,-0.9)){
    if(10>0){
      filename=paste("D:/art sim/result/d=",d,"/d=",d,"betaXYnot0rho",rho,".RData",sep="")
      setwd(paste("D:/art sim/","d=",d,sep=""))
      for(i in 1:N){
        name[i]=paste("sim",rho,"_",i,"betaXYnot",betaXY,".txt",sep="")
      }
      result=list()
      for(i in 1:N){
        result[[i]]=read.table(file=name[i],header = TRUE,sep="\t")
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
      
      
      all_snps[[1]]=data.frame(rs=rep(0,times=1000),beta_before=rep(0,times=1000),sd_before=rep(0,times=1000),p_before=rep(0,times=1000),beta_after=rep(0,times=1000),sd_after=rep(0,times=1000),p_after=rep(0,times=1000),beta_egger=rep(0,times=1000),sd_egger=rep(0,times=1000),p_egger=rep(0,times=1000),beta_true=rep(0,times=1000),bias=rep(0,times=1000),bias_true=rep(0,times=1000))
      
      all_snps[[1]]=cbind(all_snps[[1]],matrix(0,ncol=3*d,nrow=nrow(all_snps[[1]])))
      
      
      
      for(j in  1:1000){
        all_snps[[1]][j,]=a[1+(j-1)*1000,]
      }
      data=all_snps[[1]][,(ncol(all_snps[[1]])-3*d+1):ncol(all_snps[[1]])]
      
      d_mean=colMeans(data)
      d_sd=numeric()
      for(i in 1:(3*d)){
        d_sd[i]=sd(data[,i])/sqrt(1000)
      }
      
      d_true=-d_mean[1:d]
      d_cml=d_mean[(d+1):(2*d)]
      d_egger=d_mean[(2*d+1):(3*d)]
      
      sd_true= d_sd[1:d]
      sd_cml= d_sd[(d+1):(2*d)]
      sd_egger= d_sd[(2*d+1):(3*d)]
      result=list(data.frame(d_true,sd_true,d_cml,sd_cml,d_egger,sd_egger))
      save(result,file=filename)
    }
  }
}

