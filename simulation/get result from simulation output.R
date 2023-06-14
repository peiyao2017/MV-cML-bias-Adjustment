


N=20
name=numeric()
betaXY=0

for(d in c(1,2,4,6,8)){
  for(rho in c(0,0.5,0.9,-0.5,-0.9)){
    if(10>0){
      setwd(paste("D:/art sim/","d=",d,sep=""))
plotname=numeric()
plotname[1]=paste("D:/art sim/result/d=",d,"/d=",d,"betaXYis0rho",rho,"_before.jpeg",sep="")
plotname[2]=paste("D:/art sim/result/d=",d,"/d=",d,"betaXYis0rho",rho,"_after.jpeg",sep="")
plotname[3]=paste("D:/art sim/result/d=",d,"/d=",d,"betaXYis0rho",rho,"_egger.jpeg",sep="")
filename=paste("D:/art sim/result/d=",d,"/d=",d,"betaXYis0rho",rho,".txt",sep="")
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

for(i in 1:1000){
  all_snps[[i]]=data.frame(rs=rep(0,times=1000),beta_before=rep(0,times=1000),sd_before=rep(0,times=1000),p_before=rep(0,times=1000),beta_after=rep(0,times=1000),sd_after=rep(0,times=1000),p_after=rep(0,times=1000),beta_egger=rep(0,times=1000),sd_egger=rep(0,times=1000),p_egger=rep(0,times=1000),beta_true=rep(0,times=1000),bias=rep(0,times=1000),bias_true=rep(0,times=1000))
  all_snps[[i]]=cbind(all_snps[[i]],matrix(0,ncol=3*d,nrow=nrow(all_snps[[i]])))
  }
for(i in 1:1000){
  for(j in  1:1000){
    all_snps[[i]][j,]=a[i+(j-1)*1000,]
  }
}
null_snps=list()
for(i in 1:850){
  null_snps[[i]]=all_snps[[i+150]]
}

all_s
all_snps_no_y=list()
for(i in 1:900){
  if(i<=50){
    all_snps_no_y[[i]]=all_snps[[i]]
  }
  if(i>50){
    all_snps_no_y[[i]]=all_snps[[(i-50)+150]]
  }
}
all_snps_H_no_y=list()
for(i in 1:50){
  all_snps_H_no_y[[i]]=all_snps[[i]]
}

all_snps_y=list()
for(i in 1:100){
  all_snps_y[[i]]=all_snps[[50+i]]
}
all_snps_H_and_y=list()
for(i in 1:50){
  all_snps_H_and_y[[i]]=all_snps[[100+i]]
}

a1=numeric()
vara1=numeric()
for(i in 1: length(all_snps_no_y)){
  a1[i]=mean(all_snps_no_y[[i]]$p_before<0.05)
  vara1[i]=a1[i]*(1-a1[i])/1000
}
a11=mean(a1)
sda11=sqrt(mean(vara1))
a2=numeric()
vara2=numeric()
for(i in 1: length(all_snps_H_no_y)){
  a2[i]=mean(all_snps_H_no_y[[i]]$p_before<0.05)
  vara2[i]=a2[i]*(1-a2[i])/1000
}
a21=mean(a2)
sda21=sqrt(mean(vara2))


a3=numeric()
for(i in 1:length(all_snps_no_y)){
  a3[i]=mean(all_snps_no_y[[i]]$p_before<0.05)
}
a31=max(a3)
sda31=sqrt(a31*(1-a31)/1000)
a4=numeric()
for(i in 1:1000){
  a=numeric()
  for(j in 1:50){
    a[j]=all_snps_H_no_y[[j]]$p_before[i]
  }
  a4[i]=(sum(a<0.05/50)>0)
}
a41=mean(a4)
sda41=sqrt(mean(a4*(1-a4)/1000))
a5=numeric()
vara5=numeric()
for(i in 1:length(all_snps_y)){
  a5[i]=mean(all_snps_y[[i]]$p_before<0.05)
  vara5[i]=(a5[i]*(1-a5[i])/1000)
}
a51=mean(a5)
sda51=sqrt(mean(vara5))
a6=numeric()
vara6=numeric()
for(i in 1:length(all_snps_H_and_y)){
  a6[i]=mean(all_snps_H_and_y[[i]]$p_before<0.05)
  vara6[i]=(a6[i]*(1-a6[i])/1000)
}
a61=mean(a6)
sda61=sqrt(mean(vara6))
power1=numeric()
power2=numeric()

for(i in 1:length(all_snps_y)){
  power1[i]=mean(all_snps_y[[i]]$p_before<0.05)
  power2[i]=mean(all_snps_y[[i]]$p_after<0.05)
}
change=power2-power1
if(sum(change>0)>0){
  a71=power1[which.max(change)]
}
if(sum(change>0)==0){
  a71=NA
}
sda71=sqrt(a71*(1-a71)/1000)
if(sum(change<0)>0){
  a81=power1[which.min(change)]
}
if(sum(change<0)==0){
  a81=NA
}
sda81=sqrt(a81*(1-a81)/1000)

a9=numeric()
vara9=numeric()
for(i in 1: length(null_snps)){
  a9[i]=mean(null_snps[[i]]$p_before<0.05)
  vara9[i]=a9[i]*(1-a9[i])/1000
}
a91=mean(a9)
sda91=sqrt(mean(vara9))





b1=numeric()
varb1=numeric()
for(i in 1: length(all_snps_no_y)){
  b1[i]=mean(all_snps_no_y[[i]]$p_after<0.05)
  varb1[i]=(1-b1[i])*b1[i]/1000
}
b11=mean(b1)
sdb11=sqrt(mean(varb1))
b2=numeric()
varb2=numeric()
for(i in 1: length(all_snps_H_no_y)){
  b2[i]=mean(all_snps_H_no_y[[i]]$p_after<0.05)
  varb2[i]=(1-b2[i])*b2[i]/1000
}
b21=mean(b2)
sdb21=sqrt(mean(varb2))


a=numeric()
b=numeric()
for(i in 1:length(all_snps_no_y)){
  a[i]=mean(all_snps_no_y[[i]]$p_before<0.05)
  b[i]=mean(all_snps_no_y[[i]]$p_after<0.05)
}
b31=b[which.max(a)]
sdb31=sqrt(b31*(1-b31)/1000)


b4=numeric()
for(i in 1:1000){
  a=numeric()
  for(j in 1:50){
    a[j]=all_snps_H_no_y[[j]]$p_after[i]
  }
  b4[i]=(sum(a<0.05/50)>0)
}
b41=mean(b4)
sdb41=sqrt(b41*(1-b41)/1000)

b5=numeric()
varb5=numeric()
for(i in 1:length(all_snps_y)){
  b5[i]=mean(all_snps_y[[i]]$p_after<0.05)
  varb5[i]= b5[i]*(1- b5[i])/1000
}
b51=mean(b5)
sdb51=sqrt(mean(varb5))

b6=numeric()
varb6=numeric()
for(i in 1:length(all_snps_H_and_y)){
  b6[i]=mean(all_snps_H_and_y[[i]]$p_after<0.05)
  varb6[i]= b6[i]*(1- b6[i])/1000
}
sdb61=sqrt(mean(varb6))
b61=mean(b6)
power1=numeric()
power2=numeric()



for(i in 1:length(all_snps_y)){
  power1[i]=mean(all_snps_y[[i]]$p_before<0.05)
  power2[i]=mean(all_snps_y[[i]]$p_after<0.05)
}
change=power2-power1
if(sum(change>0)>0){
  b71=power2[which.max(change)]
}
if(sum(change>0)==0){
  b71=NA
}
sdb71=(b71*(1-b71))/1000
if(sum(change<0)>0){
  b81=power2[which.min(change)]
}
if(sum(change<0)==0){
  b81=NA
}
sdb81=(b81*(1-b81))/1000

b9=numeric()
varb9=numeric()
for(i in 1: length(null_snps)){
  b9[i]=mean(null_snps[[i]]$p_after<0.05)
  varb9[i]=b9[i]*(1-b9[i])/1000
}
b91=mean(b9)
sdb91=sqrt(mean(varb9))









c1=numeric()
varc1=numeric()
for(i in 1: length(all_snps_no_y)){
  c1[i]=mean(all_snps_no_y[[i]]$p_egger<0.05)
  varc1[i]=(1-c1[i])*c1[i]/1000
}
c11=mean(c1)
sdc11=sqrt(mean(varc1))
c2=numeric()
varc2=numeric()
for(i in 1: length(all_snps_H_no_y)){
  c2[i]=mean(all_snps_H_no_y[[i]]$p_egger<0.05)
  varc2[i]=(1-c2[i])*c2[i]/1000
}
c21=mean(c2)
sdc21=sqrt(mean(varc2))


a=numeric()
c=numeric()
for(i in 1:length(all_snps_no_y)){
  a[i]=mean(all_snps_no_y[[i]]$p_before<0.05)
  c[i]=mean(all_snps_no_y[[i]]$p_egger<0.05)
}
c31=c[which.max(a)]
sdc31=sqrt(c31*(1-c31)/1000)


c4=numeric()
for(i in 1:1000){
  a=numeric()
  for(j in 1:50){
    a[j]=all_snps_H_no_y[[j]]$p_egger[i]
  }
  c4[i]=(sum(a<0.05/50)>0)
}
c41=mean(c4)
sdc41=sqrt(c41*(1-c41)/1000)

c5=numeric()
varc5=numeric()
for(i in 1:length(all_snps_y)){
  c5[i]=mean(all_snps_y[[i]]$p_egger<0.05)
  varc5[i]= c5[i]*(1- c5[i])/1000
}
c51=mean(c5)
sdc51=sqrt(mean(varc5))

c6=numeric()
varc6=numeric()
for(i in 1:length(all_snps_H_and_y)){
  c6[i]=mean(all_snps_H_and_y[[i]]$p_egger<0.05)
  varc6[i]= c6[i]*(1- c6[i])/1000
}
sdc61=sqrt(mean(varc6))
c61=mean(c6)
power1=numeric()
power2=numeric()



for(i in 1:length(all_snps_y)){
  power1[i]=mean(all_snps_y[[i]]$p_before<0.05)
  power2[i]=mean(all_snps_y[[i]]$p_egger<0.05)
}
change=power2-power1
if(sum(change>0)>0){
  c71=power2[which.max(change)]
}
if(sum(change>0)==0){
  c71=NA
}
sdc71=(c71*(1-c71))/1000
if(sum(change<0)>0){
  c81=power2[which.min(change)]
}
if(sum(change<0)==0){
  c81=NA
}
sdc81=(c81*(1-c81))/1000
c9=numeric()
varc9=numeric()
for(i in 1: length(null_snps)){
  c9[i]=mean(null_snps[[i]]$p_egger<0.05)
  varc9[i]=c9[i]*(1-c9[i])/1000
}
c91=mean(c9)
sdc91=sqrt(mean(varc9))








no=c(a11,a21,a31,a41,a51,a61,a71,a81,a91)
sdno=c(sda11,sda21,sda31,sda41,sda51,sda61,sda71,sda81,sda91)
yes=c(b11,b21,b31,b41,b51,b61,b71,b81,b91)
sdyes=c(sdb11,sdb21,sdb31,sdb41,sdb51,sdb61,sdb71,sdb81,sdb91)
egger=c(c11,c21,c31,c41,c51,c61,c71,c81,c91)
sdegger=c(sdc11,sdc21,sdc31,sdc41,sdc51,sdc61,sdc71,sdc81,sdc91)
result=data.frame(no=round(no,digits = 3),sdno=round(sdno,digits = 3),yes=round(yes,digits = 3),sdyes=round(sdyes,digits = 3),egger=round(egger,digits = 3),sdegger=round(sdegger,digits = 3))









beta_before=numeric()
beta_after=numeric()
beta_true=numeric()
beta_egger=numeric()
sd_egger=numeric()
sd_before=numeric()
sd_after=numeric()
x0_before=numeric()
y0_before=numeric()
x1_before=numeric()
y1_before=numeric()
x0_after=numeric()
y0_after=numeric()
x1_after=numeric()
y1_after=numeric()
x0_egger=numeric()
y0_egger=numeric()
x1_egger=numeric()
y1_egger=numeric()
for(i in 1:1000){
  beta_before[i]=mean(all_snps[[i]]$beta_before)
  beta_after[i]=mean(all_snps[[i]]$beta_after)
  sd_before[i]=mean(all_snps[[i]]$sd_before)
  sd_after[i]=mean(all_snps[[i]]$sd_after)
  beta_true[i]=mean(all_snps[[i]]$beta_true)
  beta_egger[i]=mean(all_snps[[i]]$beta_egger)
  sd_egger[i]=mean(all_snps[[i]]$sd_egger)
}
for(i in 1:1000){
  x0_before[i]=beta_true[i]
  y0_before[i]=beta_before[i]-sd_before[i]
  x1_before[i]=beta_true[i]
  y1_before[i]=beta_before[i]+sd_before[i]
  x0_after[i]=beta_true[i]
  y0_after[i]=beta_after[i]-sd_after[i]
  x1_after[i]=beta_true[i]
  y1_after[i]=beta_after[i]+sd_after[i]
  x0_egger[i]=beta_true[i]
  y0_egger[i]=beta_egger[i]-sd_egger[i]
  x1_egger[i]=beta_true[i]
  y1_egger[i]=beta_egger[i]+sd_egger[i]
}

jpeg(filename = plotname[1],width=500,height = 500,res=100)
  plot(beta_true,beta_before,xlim=c(min(beta_true),max(beta_true)),ylim=c(min(y0_before,y0_after,y0_egger),max(y1_before,y1_after,y1_egger)),main="SNP effects before adjustment",type="n",xlab = "true effect",ylab = "estimated effect")

legend("topleft",pch=c(15,16,17,18),col=c("black","darkorchid1","gold1","lightskyblue"),legend=c("H only","y only","H and Y","null SNPs"),cex=c(0.7,0.7,0.7,0.7))
abline(a=0,b=1)
for(i in 1:1000){
  segments(x0=x0_before[i],y0=y0_before[i],x1=x1_before[i],y1=y1_before[i],lwd=0.5,col="gray")
}
points(beta_true[1:50],beta_before[1:50],col=1,pch=15,cex=0.7)
points(beta_true[51:100],beta_before[51:100],col="darkorchid1",pch=16,cex=0.7)
points(beta_true[101:150],beta_before[101:150],col="gold1",pch=17,cex=0.7)
points(beta_true[151:1000],beta_before[151:1000],col="lightskyblue",pch=18,cex=0.7)
dev.off()









jpeg(filename = plotname[2],width=500,height = 500,res=100)

  plot(beta_true,beta_after,xlim=c(min(beta_true),max(beta_true)),ylim=c(min(y0_before,y0_after,y0_egger),max(y1_before,y1_after,y1_egger)),main="SNP effects after adjustment",type="n",xlab = "true effect",ylab = "estimated effect")

legend("topleft",pch=c(15,16,17,18),col=c("black","darkorchid1","gold1","lightskyblue"),legend=c("H only","y only","H and Y","null SNPs"),cex=c(0.7,0.7,0.7,0.7))
abline(a=0,b=1)
for(i in 1:1000){
  segments(x0=x0_after[i],y0=y0_after[i],x1=x1_after[i],y1=y1_after[i],lwd=0.5,col="gray")
}


points(beta_true[1:50],beta_after[1:50],col=1,pch=15,cex=0.7)
points(beta_true[51:100],beta_after[51:100],col="darkorchid1",pch=16,cex=0.7)
points(beta_true[101:150],beta_after[101:150],col="gold1",pch=17,cex=0.7)
points(beta_true[151:1000],beta_after[151:1000],col="lightskyblue",pch=18,cex=0.7)
dev.off()



jpeg(filename = plotname[3],width=500,height = 500,res=100)
 
  plot(beta_true,beta_egger,xlim=c(min(beta_true),max(beta_true)),ylim=c(min(y0_before,y0_egger,y0_egger),max(y1_before,y1_egger,y1_egger)),main="SNP effects MV-Egger adjustment",type="n",xlab = "true effect",ylab = "estimated effect")
 
legend("topleft",pch=c(15,16,17,18),col=c("black","darkorchid1","gold1","lightskyblue"),legend=c("H only","y only","H and Y","null SNPs"),cex=c(0.7,0.7,0.7,0.7))
abline(a=0,b=1)
for(i in 1:1000){
  segments(x0=x0_egger[i],y0=y0_egger[i],x1=x1_egger[i],y1=y1_egger[i],lwd=0.5,col="gray")
}


points(beta_true[1:50],beta_egger[1:50],col=1,pch=15,cex=0.7)
points(beta_true[51:100],beta_egger[51:100],col="darkorchid1",pch=16,cex=0.7)
points(beta_true[101:150],beta_egger[101:150],col="gold1",pch=17,cex=0.7)
points(beta_true[151:1000],beta_egger[151:1000],col="lightskyblue",pch=18,cex=0.7)
dev.off()
write.table(result,file = filename,sep="\t",col.names = TRUE,quote = FALSE,row.names = FALSE)
}
  }
}




for(d in c(1,2,4,6,8)){
  for(rho in c(0,0.5,0.9,-0.5,-0.9)){
    if(10>0){
      setwd(paste("D:/art sim/","d=",d,sep=""))
      plotname=numeric()
      plotname[1]=paste("D:/art sim/result/d=",d,"/d=",d,"betaXYnot0rho",rho,"_before.jpeg",sep="")
      plotname[2]=paste("D:/art sim/result/d=",d,"/d=",d,"betaXYnot0rho",rho,"_after.jpeg",sep="")
      plotname[3]=paste("D:/art sim/result/d=",d,"/d=",d,"betaXYnot0rho",rho,"_egger.jpeg",sep="")
      filename=paste("D:/art sim/result/d=",d,"/d=",d,"betaXYnot0rho",rho,".txt",sep="")
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
      
      for(i in 1:1000){
        all_snps[[i]]=data.frame(rs=rep(0,times=1000),beta_before=rep(0,times=1000),sd_before=rep(0,times=1000),p_before=rep(0,times=1000),beta_after=rep(0,times=1000),sd_after=rep(0,times=1000),p_after=rep(0,times=1000),beta_egger=rep(0,times=1000),sd_egger=rep(0,times=1000),p_egger=rep(0,times=1000),beta_true=rep(0,times=1000),bias=rep(0,times=1000),bias_true=rep(0,times=1000))
        all_snps[[i]]=cbind(all_snps[[i]],matrix(0,ncol=3*d,nrow=nrow(all_snps[[i]])))
        }
      for(i in 1:1000){
        for(j in  1:1000){
          all_snps[[i]][j,]=a[i+(j-1)*1000,]
        }
      }
      null_snps=list()
      for(i in 1:850){
        null_snps[[i]]=all_snps[[i+150]]
      }
      
      all_s
      all_snps_no_y=list()
      for(i in 1:900){
        if(i<=50){
          all_snps_no_y[[i]]=all_snps[[i]]
        }
        if(i>50){
          all_snps_no_y[[i]]=all_snps[[(i-50)+150]]
        }
      }
      all_snps_H_no_y=list()
      for(i in 1:50){
        all_snps_H_no_y[[i]]=all_snps[[i]]
      }
      
      all_snps_y=list()
      for(i in 1:100){
        all_snps_y[[i]]=all_snps[[50+i]]
      }
      all_snps_H_and_y=list()
      for(i in 1:50){
        all_snps_H_and_y[[i]]=all_snps[[100+i]]
      }
      
      a1=numeric()
      vara1=numeric()
      for(i in 1: length(all_snps_no_y)){
        a1[i]=mean(all_snps_no_y[[i]]$p_before<0.05)
        vara1[i]=a1[i]*(1-a1[i])/1000
      }
      a11=mean(a1)
      sda11=sqrt(mean(vara1))
      a2=numeric()
      vara2=numeric()
      for(i in 1: length(all_snps_H_no_y)){
        a2[i]=mean(all_snps_H_no_y[[i]]$p_before<0.05)
        vara2[i]=a2[i]*(1-a2[i])/1000
      }
      a21=mean(a2)
      sda21=sqrt(mean(vara2))
      
      
      a3=numeric()
      for(i in 1:length(all_snps_no_y)){
        a3[i]=mean(all_snps_no_y[[i]]$p_before<0.05)
      }
      a31=max(a3)
      sda31=sqrt(a31*(1-a31)/1000)
      a4=numeric()
      for(i in 1:1000){
        a=numeric()
        for(j in 1:50){
          a[j]=all_snps_H_no_y[[j]]$p_before[i]
        }
        a4[i]=(sum(a<0.05/50)>0)
      }
      a41=mean(a4)
      sda41=sqrt(mean(a4*(1-a4)/1000))
      a5=numeric()
      vara5=numeric()
      for(i in 1:length(all_snps_y)){
        a5[i]=mean(all_snps_y[[i]]$p_before<0.05)
        vara5[i]=(a5[i]*(1-a5[i])/1000)
      }
      a51=mean(a5)
      sda51=sqrt(mean(vara5))
      a6=numeric()
      vara6=numeric()
      for(i in 1:length(all_snps_H_and_y)){
        a6[i]=mean(all_snps_H_and_y[[i]]$p_before<0.05)
        vara6[i]=(a6[i]*(1-a6[i])/1000)
      }
      a61=mean(a6)
      sda61=sqrt(mean(vara6))
      power1=numeric()
      power2=numeric()
      
      for(i in 1:length(all_snps_y)){
        power1[i]=mean(all_snps_y[[i]]$p_before<0.05)
        power2[i]=mean(all_snps_y[[i]]$p_after<0.05)
      }
      change=power2-power1
      if(sum(change>0)>0){
        a71=power1[which.max(change)]
      }
      if(sum(change>0)==0){
        a71=NA
      }
      sda71=sqrt(a71*(1-a71)/1000)
      if(sum(change<0)>0){
        a81=power1[which.min(change)]
      }
      if(sum(change<0)==0){
        a81=NA
      }
      sda81=sqrt(a81*(1-a81)/1000)
      
      a9=numeric()
      vara9=numeric()
      for(i in 1: length(null_snps)){
        a9[i]=mean(null_snps[[i]]$p_before<0.05)
        vara9[i]=a9[i]*(1-a9[i])/1000
      }
      a91=mean(a9)
      sda91=sqrt(mean(vara9))
      
      
      
      
      
      b1=numeric()
      varb1=numeric()
      for(i in 1: length(all_snps_no_y)){
        b1[i]=mean(all_snps_no_y[[i]]$p_after<0.05)
        varb1[i]=(1-b1[i])*b1[i]/1000
      }
      b11=mean(b1)
      sdb11=sqrt(mean(varb1))
      b2=numeric()
      varb2=numeric()
      for(i in 1: length(all_snps_H_no_y)){
        b2[i]=mean(all_snps_H_no_y[[i]]$p_after<0.05)
        varb2[i]=(1-b2[i])*b2[i]/1000
      }
      b21=mean(b2)
      sdb21=sqrt(mean(varb2))
      
      
      a=numeric()
      b=numeric()
      for(i in 1:length(all_snps_no_y)){
        a[i]=mean(all_snps_no_y[[i]]$p_before<0.05)
        b[i]=mean(all_snps_no_y[[i]]$p_after<0.05)
      }
      b31=b[which.max(a)]
      sdb31=sqrt(b31*(1-b31)/1000)
      
      
      b4=numeric()
      for(i in 1:1000){
        a=numeric()
        for(j in 1:50){
          a[j]=all_snps_H_no_y[[j]]$p_after[i]
        }
        b4[i]=(sum(a<0.05/50)>0)
      }
      b41=mean(b4)
      sdb41=sqrt(b41*(1-b41)/1000)
      
      b5=numeric()
      varb5=numeric()
      for(i in 1:length(all_snps_y)){
        b5[i]=mean(all_snps_y[[i]]$p_after<0.05)
        varb5[i]= b5[i]*(1- b5[i])/1000
      }
      b51=mean(b5)
      sdb51=sqrt(mean(varb5))
      
      b6=numeric()
      varb6=numeric()
      for(i in 1:length(all_snps_H_and_y)){
        b6[i]=mean(all_snps_H_and_y[[i]]$p_after<0.05)
        varb6[i]= b6[i]*(1- b6[i])/1000
      }
      sdb61=sqrt(mean(varb6))
      b61=mean(b6)
      power1=numeric()
      power2=numeric()
      
      
      
      for(i in 1:length(all_snps_y)){
        power1[i]=mean(all_snps_y[[i]]$p_before<0.05)
        power2[i]=mean(all_snps_y[[i]]$p_after<0.05)
      }
      change=power2-power1
      if(sum(change>0)>0){
        b71=power2[which.max(change)]
      }
      if(sum(change>0)==0){
        b71=NA
      }
      sdb71=(b71*(1-b71))/1000
      if(sum(change<0)>0){
        b81=power2[which.min(change)]
      }
      if(sum(change<0)==0){
        b81=NA
      }
      sdb81=(b81*(1-b81))/1000
      
      b9=numeric()
      varb9=numeric()
      for(i in 1: length(null_snps)){
        b9[i]=mean(null_snps[[i]]$p_after<0.05)
        varb9[i]=b9[i]*(1-b9[i])/1000
      }
      b91=mean(b9)
      sdb91=sqrt(mean(varb9))
      
      
      
      
      
      
      
      
      
      c1=numeric()
      varc1=numeric()
      for(i in 1: length(all_snps_no_y)){
        c1[i]=mean(all_snps_no_y[[i]]$p_egger<0.05)
        varc1[i]=(1-c1[i])*c1[i]/1000
      }
      c11=mean(c1)
      sdc11=sqrt(mean(varc1))
      c2=numeric()
      varc2=numeric()
      for(i in 1: length(all_snps_H_no_y)){
        c2[i]=mean(all_snps_H_no_y[[i]]$p_egger<0.05)
        varc2[i]=(1-c2[i])*c2[i]/1000
      }
      c21=mean(c2)
      sdc21=sqrt(mean(varc2))
      
      
      a=numeric()
      c=numeric()
      for(i in 1:length(all_snps_no_y)){
        a[i]=mean(all_snps_no_y[[i]]$p_before<0.05)
        c[i]=mean(all_snps_no_y[[i]]$p_egger<0.05)
      }
      c31=c[which.max(a)]
      sdc31=sqrt(c31*(1-c31)/1000)
      
      
      c4=numeric()
      for(i in 1:1000){
        a=numeric()
        for(j in 1:50){
          a[j]=all_snps_H_no_y[[j]]$p_egger[i]
        }
        c4[i]=(sum(a<0.05/50)>0)
      }
      c41=mean(c4)
      sdc41=sqrt(c41*(1-c41)/1000)
      
      c5=numeric()
      varc5=numeric()
      for(i in 1:length(all_snps_y)){
        c5[i]=mean(all_snps_y[[i]]$p_egger<0.05)
        varc5[i]= c5[i]*(1- c5[i])/1000
      }
      c51=mean(c5)
      sdc51=sqrt(mean(varc5))
      
      c6=numeric()
      varc6=numeric()
      for(i in 1:length(all_snps_H_and_y)){
        c6[i]=mean(all_snps_H_and_y[[i]]$p_egger<0.05)
        varc6[i]= c6[i]*(1- c6[i])/1000
      }
      sdc61=sqrt(mean(varc6))
      c61=mean(c6)
      power1=numeric()
      power2=numeric()
      
      
      
      for(i in 1:length(all_snps_y)){
        power1[i]=mean(all_snps_y[[i]]$p_before<0.05)
        power2[i]=mean(all_snps_y[[i]]$p_egger<0.05)
      }
      change=power2-power1
      if(sum(change>0)>0){
        c71=power2[which.max(change)]
      }
      if(sum(change>0)==0){
        c71=NA
      }
      sdc71=(c71*(1-c71))/1000
      if(sum(change<0)>0){
        c81=power2[which.min(change)]
      }
      if(sum(change<0)==0){
        c81=NA
      }
      sdc81=(c81*(1-c81))/1000
      c9=numeric()
      varc9=numeric()
      for(i in 1: length(null_snps)){
        c9[i]=mean(null_snps[[i]]$p_egger<0.05)
        varc9[i]=c9[i]*(1-c9[i])/1000
      }
      c91=mean(c9)
      sdc91=sqrt(mean(varc9))
      
      
      
      
      
      
      
      
      no=c(a11,a21,a31,a41,a51,a61,a71,a81,a91)
      sdno=c(sda11,sda21,sda31,sda41,sda51,sda61,sda71,sda81,sda91)
      yes=c(b11,b21,b31,b41,b51,b61,b71,b81,b91)
      sdyes=c(sdb11,sdb21,sdb31,sdb41,sdb51,sdb61,sdb71,sdb81,sdb91)
      egger=c(c11,c21,c31,c41,c51,c61,c71,c81,c91)
      sdegger=c(sdc11,sdc21,sdc31,sdc41,sdc51,sdc61,sdc71,sdc81,sdc91)
      result=data.frame(no=round(no,digits = 3),sdno=round(sdno,digits = 3),yes=round(yes,digits = 3),sdyes=round(sdyes,digits = 3),egger=round(egger,digits = 3),sdegger=round(sdegger,digits = 3))
      
      
      
      
      
      
      
      
      
      beta_before=numeric()
      beta_after=numeric()
      beta_true=numeric()
      beta_egger=numeric()
      sd_egger=numeric()
      sd_before=numeric()
      sd_after=numeric()
      x0_before=numeric()
      y0_before=numeric()
      x1_before=numeric()
      y1_before=numeric()
      x0_after=numeric()
      y0_after=numeric()
      x1_after=numeric()
      y1_after=numeric()
      x0_egger=numeric()
      y0_egger=numeric()
      x1_egger=numeric()
      y1_egger=numeric()
      for(i in 1:1000){
        beta_before[i]=mean(all_snps[[i]]$beta_before)
        beta_after[i]=mean(all_snps[[i]]$beta_after)
        sd_before[i]=mean(all_snps[[i]]$sd_before)
        sd_after[i]=mean(all_snps[[i]]$sd_after)
        beta_true[i]=mean(all_snps[[i]]$beta_true)
        beta_egger[i]=mean(all_snps[[i]]$beta_egger)
        sd_egger[i]=mean(all_snps[[i]]$sd_egger)
      }
      for(i in 1:1000){
        x0_before[i]=beta_true[i]
        y0_before[i]=beta_before[i]-sd_before[i]
        x1_before[i]=beta_true[i]
        y1_before[i]=beta_before[i]+sd_before[i]
        x0_after[i]=beta_true[i]
        y0_after[i]=beta_after[i]-sd_after[i]
        x1_after[i]=beta_true[i]
        y1_after[i]=beta_after[i]+sd_after[i]
        x0_egger[i]=beta_true[i]
        y0_egger[i]=beta_egger[i]-sd_egger[i]
        x1_egger[i]=beta_true[i]
        y1_egger[i]=beta_egger[i]+sd_egger[i]
      }
      
      jpeg(filename = plotname[1],width=500,height = 500,res=100)
      plot(beta_true,beta_before,xlim=c(min(beta_true),max(beta_true)),ylim=c(min(y0_before,y0_after,y0_egger),max(y1_before,y1_after,y1_egger)),main="SNP effects before adjustment",type="n",xlab = "true effect",ylab = "estimated effect")
      
      legend("topleft",pch=c(15,16,17,18),col=c("black","darkorchid1","gold1","lightskyblue"),legend=c("H only","y only","H and Y","null SNPs"),cex=c(0.7,0.7,0.7,0.7))
      abline(a=0,b=1)
      for(i in 1:1000){
        segments(x0=x0_before[i],y0=y0_before[i],x1=x1_before[i],y1=y1_before[i],lwd=0.5,col="gray")
      }
      points(beta_true[1:50],beta_before[1:50],col=1,pch=15,cex=0.7)
      points(beta_true[51:100],beta_before[51:100],col="darkorchid1",pch=16,cex=0.7)
      points(beta_true[101:150],beta_before[101:150],col="gold1",pch=17,cex=0.7)
      points(beta_true[151:1000],beta_before[151:1000],col="lightskyblue",pch=18,cex=0.7)
      dev.off()
      
      
      
      
      
      
      
      
      
      jpeg(filename = plotname[2],width=500,height = 500,res=100)
      
      plot(beta_true,beta_after,xlim=c(min(beta_true),max(beta_true)),ylim=c(min(y0_before,y0_after,y0_egger),max(y1_before,y1_after,y1_egger)),main="SNP effects after adjustment",type="n",xlab = "true effect",ylab = "estimated effect")
      
      legend("topleft",pch=c(15,16,17,18),col=c("black","darkorchid1","gold1","lightskyblue"),legend=c("H only","y only","H and Y","null SNPs"),cex=c(0.7,0.7,0.7,0.7))
      abline(a=0,b=1)
      for(i in 1:1000){
        segments(x0=x0_after[i],y0=y0_after[i],x1=x1_after[i],y1=y1_after[i],lwd=0.5,col="gray")
      }
      
      
      points(beta_true[1:50],beta_after[1:50],col=1,pch=15,cex=0.7)
      points(beta_true[51:100],beta_after[51:100],col="darkorchid1",pch=16,cex=0.7)
      points(beta_true[101:150],beta_after[101:150],col="gold1",pch=17,cex=0.7)
      points(beta_true[151:1000],beta_after[151:1000],col="lightskyblue",pch=18,cex=0.7)
      dev.off()
      
      
      
      jpeg(filename = plotname[3],width=500,height = 500,res=100)
      
      plot(beta_true,beta_egger,xlim=c(min(beta_true),max(beta_true)),ylim=c(min(y0_before,y0_egger,y0_egger),max(y1_before,y1_egger,y1_egger)),main="SNP effects MV-Egger adjustment",type="n",xlab = "true effect",ylab = "estimated effect")
      
      legend("topleft",pch=c(15,16,17,18),col=c("black","darkorchid1","gold1","lightskyblue"),legend=c("H only","y only","H and Y","null SNPs"),cex=c(0.7,0.7,0.7,0.7))
      abline(a=0,b=1)
      for(i in 1:1000){
        segments(x0=x0_egger[i],y0=y0_egger[i],x1=x1_egger[i],y1=y1_egger[i],lwd=0.5,col="gray")
      }
      
      
      points(beta_true[1:50],beta_egger[1:50],col=1,pch=15,cex=0.7)
      points(beta_true[51:100],beta_egger[51:100],col="darkorchid1",pch=16,cex=0.7)
      points(beta_true[101:150],beta_egger[101:150],col="gold1",pch=17,cex=0.7)
      points(beta_true[151:1000],beta_egger[151:1000],col="lightskyblue",pch=18,cex=0.7)
      dev.off()
      write.table(result,file = filename,sep="\t",col.names = TRUE,quote = FALSE,row.names = FALSE)
    }
  }
}



