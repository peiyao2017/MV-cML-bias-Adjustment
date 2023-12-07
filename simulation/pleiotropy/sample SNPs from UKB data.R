library(BEDMatrix)
setwd("D:/art sim/New/real_snps_sim")
nIV=30
Nsnps=1000
Gmatrixtotal=BEDMatrix("D:/art sim/New/newGWAS/110K_QCed1.bed",simple_names = TRUE)
Indsnps=read.table("D:/art sim/New/newGWAS/110K_QCed0.001.bim",header = FALSE)
snps=read.table("D:/art sim/New/newGWAS/110K_QCed1.bim",header = FALSE)

IVID=sample(Indsnps$V2[Indsnps$V1==6],size=nIV,replace = FALSE)
ID=sample(snps$V2[snps$V1==6&!(snps$V2%in%IVID)],size=Nsnps-nIV,replace = FALSE)
invalid=0.3
vp=1-invalid
X_only=c(ID[1:(Nsnps*0.05-nIV*vp)],IVID[1:(nIV*vp)])
n1=(Nsnps*0.05-nIV*vp)+1
n2=(Nsnps*0.05-nIV*vp)+Nsnps*0.05
Y_only=ID[n1:n2]
n3=n2+1
n4=n2+Nsnps*0.05-nIV*(1-vp)
n5=nIV*vp+1
n6=nIV
X_and_Y=c(ID[n3:n4],IVID[n5:n6])
null=ID[(n4+1):length(ID)]

simID=c(X_only,Y_only,X_and_Y,null)
simIV=IVID
Gmatrix=Gmatrixtotal[,simID]
for(i in 1:ncol(Gmatrix)){
  Gmatrix[,i][is.na( Gmatrix[,i])]=mean(Gmatrix[,i],na.rm=TRUE)
}

Geno=list(Gmatrix,simID,simIV)
save(Geno,file = "Geno30.RDATA")



for(d in c(1)){
  for(rho in c(0,0.5,-0.5)){
    if(10>0){
      N=20
      setwd(paste("D:/art sim/New/real_snps_sim/50invalid/","d=",d,sep=""))
      plotname=numeric()
      name1_1=numeric()
      plotname[1]=paste("D:/art sim/New/real_snps_sim/result/50invalid/","d=",d,"/","d=",d,"betaXYnot0rho",rho,"_before.jpeg",sep="")
      plotname[2]=paste("D:/art sim/New/real_snps_sim/result/50invalid/","d=",d,"/","d=",d,"betaXYnot0rho",rho,"_after.jpeg",sep="")
      plotname[3]=paste("D:/art sim/New/real_snps_sim/result/50invalid/","d=",d,"/","d=",d,"betaXYnot0rho",rho,"_egger.jpeg",sep="")
      plotname[4]=paste("D:/art sim/New/real_snps_sim/result/50invalid/","d=",d,"/","d=",d,"betaXYnot0rho",rho,"_dd.jpeg",sep="")
      plotname[5]=paste("D:/art sim/New/real_snps_sim/result/50invalid/","d=",d,"/","d=",d,"betaXYnot0rho",rho,"_sh.jpeg",sep="")
      plotname[6]=paste("D:/art sim/New/real_snps_sim/result/50invalid/","d=",d,"/","d=",d,"betaXYnot0rho",rho,"_ivw.jpeg",sep="")
      plotname[7]=paste("D:/art sim/New/real_snps_sim/result/50invalid/","d=",d,"/","d=",d,"betaXYnot0rho",rho,"_lasso.jpeg",sep="")
      plotname[8]=paste("D:/art sim/New/real_snps_sim/result/50invalid/","d=",d,"/","d=",d,"betaXYnot0rho",rho,"_median.jpeg",sep="")
      filename1=paste("D:/art sim/New/real_snps_sim/result/50invalid/","d=",d,"/","d=",d,"betaXYnot0rho",rho,".txt",sep="")
      filename2=paste("D:/art sim/New/real_snps_sim/result/50invalid/","d=",d,"/","d=",d,"betaXYnot0rho",rho,"Slope.txt",sep="")
      filename3=paste("D:/art sim/New/real_snps_sim/result/50invalid/","d=",d,"/","d=",d,"betaXYnot0rho",rho,"point_effect_estimates.RData",sep="")
      
      for(i in 1:N){
        name1_1[i]=paste("sim",rho,"_",i,"betaXYnot",0,".txt",sep="")
      }
      result=list()
      for(i in 1:N){
        a1=read.table(file=name1_1[i],header = TRUE,sep="\t")
        result[[i]]=a1
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
        all_snps[[i]]=a[1:1000,]
      }
      for(i in 1:1000){
        for(j in  1:1000){
          all_snps[[i]][j,]=a[i+(j-1)*1000,]
        }
        print(i)
      }
      if(10>0){
        name6se=numeric()
        name7se=numeric()
        name8se=numeric()
        for(i in 1:d){
          name6se[i]=paste("ivw_slop_se",i,sep="")
          name7se[i]=paste("lasso_slop_se",i,sep="")
          name8se[i]=paste("median_slop_se",i,sep="")
        }
        for(i in 1:length(all_snps)){
          colnames(all_snps[[i]])=c(colnames(all_snps[[i]])[1:(ncol(all_snps[[i]])-3*d)],name6se,name7se,name8se)
        }
      }
      null_snps=list()
      for(i in 1:850){
        null_snps[[i]]=all_snps[[i+150]]
      }
      
      all_snps_no_y=list()
      for(i in c(1:900)){
        if(i<=50){
          all_snps_no_y[[i]]=all_snps[[i]]
        }
        if(i>50){
          all_snps_no_y[[i]]=all_snps[[i+100]]
        }
      }
      
      
      all_snps_H_no_y=list()
      for(i in 1:50){
        all_snps_H_no_y[[i]]=all_snps[[i]]
      }
      
      all_snps_y_only=list()
      for(i in 1:50){
        all_snps_y_only[[i]]=all_snps[[50+i]]
      }
      all_snps_H_and_y=list()
      for(i in 1:50){
        all_snps_H_and_y[[i]]=all_snps[[100+i]]
      }
      
      all_snps_y=list()
      for(i in 1:100){
        all_snps_y[[i]]=all_snps[[50+i]]
      }
      
      if(10>0){
        a1=numeric()
        vara1=numeric()
        for(i in 1: length(null_snps)){
          a1[i]=mean(null_snps[[i]]$p_before<0.05)
          vara1[i]=var(null_snps[[i]]$p_before<0.05)
        }
        a11=mean(a1)
        sda11=sqrt(mean(vara1))
        a2=numeric()
        vara2=numeric()
        for(i in 1: length(all_snps_H_no_y)){
          a2[i]=mean(all_snps_H_no_y[[i]]$p_before<0.05)
          vara2[i]=var(all_snps_H_no_y[[i]]$p_before<0.05)
        }
        a21=mean(a2)
        sda21=sqrt(mean(vara2))
        
        
        a3=numeric()
        vara3=numeric()
        for(i in 1:length(all_snps_no_y)){
          a3[i]=mean(all_snps_no_y[[i]]$p_before<0.05)
          vara3[i]=var(all_snps_no_y[[i]]$p_before<0.05)
        }
        a31=max(a3)
        sda31=sqrt(vara3[which.max(a3)])
        
        
        a4=numeric()
        for(i in 1:1000){
          a=numeric()
          for(j in 1:50){
            a[j]=all_snps_H_no_y[[j]]$p_before[i]
          }
          a4[i]=(sum(a<0.05/50)>0)
        }
        a41=mean(a4)
        sda41=sd(a4) 
        
        a5=numeric()
        vara5=numeric()
        for(i in 1:length(all_snps_y_only)){
          a5[i]=mean(all_snps_y_only[[i]]$p_before<0.05)
          vara5[i]=var(all_snps_y_only[[i]]$p_before<0.05)
        }
        a51=mean(a5)
        sda51=sqrt(mean(vara5))
        a6=numeric()
        vara6=numeric()
        for(i in 1:length(all_snps_H_and_y)){
          a6[i]=mean(all_snps_H_and_y[[i]]$p_before<0.05)
          vara6[i]=var(all_snps_H_and_y[[i]]$p_before<0.05)
        }
        a61=mean(a6)
        sda61=sqrt(mean(vara6))
        power1=numeric()
        power2=numeric()
        
        for(i in 1:length(all_snps_y)){
          power1[i]=mean(all_snps_y[[i]]$p_before<0.05)
          power2[i]=mean(all_snps_y[[i]]$p_after3<0.05)
        }
        change=power2-power1
        if(sum(change>0)>0){
          a71=power1[which.max(change)]
          sda71=sd(all_snps_y[[which.max(change)]]$p_before<0.05)
        }
        if(sum(change>0)==0){
          a71=NA
          sda71=NA
        }
        
        if(sum(change<0)>0){
          a81=power1[which.min(change)]
          sda81=sd(all_snps_y[[which.min(change)]]$p_before<0.05)
        }
        if(sum(change<0)==0){
          a81=NA
          sda81=NA
        }
        
      }
      
      if(10>0){
        b1=numeric()
        varb1=numeric()
        for(i in 1: length(null_snps)){
          b1[i]=mean(null_snps[[i]]$p_after1<0.05)
          varb1[i]=var(null_snps[[i]]$p_after1<0.05)
        }
        cml11=mean(b1)
        sdcml11=sqrt(mean(varb1))
        
        
        b2=numeric()
        varb2=numeric()
        for(i in 1: length(null_snps)){
          b2[i]=mean(null_snps[[i]]$p_after2<0.05)
          varb2[i]=var(null_snps[[i]]$p_after2<0.05)
        }
        cml12=mean(b2)
        sdcml12=sqrt(mean(varb2))
        
        b3=numeric()
        varb3=numeric()
        for(i in 1: length(null_snps)){
          b3[i]=mean(null_snps[[i]]$p_after3<0.05)
          varb3[i]=var(null_snps[[i]]$p_after3<0.05)
        }
        cml13=mean(b3)
        sdcml13=sqrt(mean(varb3))
        
        
        
        
        
        
        b4=numeric()
        varb4=numeric()
        for(i in 1: length(all_snps_H_no_y)){
          b4[i]=mean(all_snps_H_no_y[[i]]$p_after1<0.05)
          varb4[i]=var(all_snps_H_no_y[[i]]$p_after1<0.05)
        }
        cml21=mean(b4)
        sdcml21=sqrt(mean(varb4))
        
        
        b5=numeric()
        varb5=numeric()
        for(i in 1: length(all_snps_H_no_y)){
          b5[i]=mean(all_snps_H_no_y[[i]]$p_after2<0.05)
          varb5[i]=var(all_snps_H_no_y[[i]]$p_after2<0.05)
        }
        cml22=mean(b5)
        sdcml22=sqrt(mean(varb5))
        
        
        b6=numeric()
        varb6=numeric()
        for(i in 1: length(all_snps_H_no_y)){
          b6[i]=mean(all_snps_H_no_y[[i]]$p_after3<0.05)
          varb6[i]=var(all_snps_H_no_y[[i]]$p_after3<0.05)
        }
        cml23=mean(b6)
        sdcml23=sqrt(mean(varb6))
        
        
        
        a=numeric()
        b=numeric()
        for(i in 1:length(all_snps_no_y)){
          a[i]=mean(all_snps_no_y[[i]]$p_before<0.05)
          b[i]=mean(all_snps_no_y[[i]]$p_after1<0.05)
        }
        cml31=b[which.max(a)]
        sdcml31=sd(all_snps_no_y[[which.max(a)]]$p_after1<0.05)
        
        a=numeric()
        b=numeric()
        for(i in 1:length(all_snps_no_y)){
          a[i]=mean(all_snps_no_y[[i]]$p_before<0.05)
          b[i]=mean(all_snps_no_y[[i]]$p_after2<0.05)
        }
        cml32=b[which.max(a)]
        sdcml32=sd(all_snps_no_y[[which.max(a)]]$p_after2<0.05)
        
        
        
        a=numeric()
        b=numeric()
        for(i in 1:length(all_snps_no_y)){
          a[i]=mean(all_snps_no_y[[i]]$p_before<0.05)
          b[i]=mean(all_snps_no_y[[i]]$p_after3<0.05)
        }
        cml33=b[which.max(a)]
        sdcml33=sd(all_snps_no_y[[which.max(a)]]$p_after3<0.05)
        
        
        b7=numeric()
        for(i in 1:1000){
          a=numeric()
          for(j in 1:50){
            a[j]=all_snps_H_no_y[[j]]$p_after1[i]
          }
          b7[i]=(sum(a<0.05/50)>0)
        }
        cml41=mean(b7)
        sdcml41=sd(b7) 
        
        b8=numeric()
        for(i in 1:1000){
          a=numeric()
          for(j in 1:50){
            a[j]=all_snps_H_no_y[[j]]$p_after2[i]
          }
          b8[i]=(sum(a<0.05/50)>0)
        }
        cml42=mean(b8)
        sdcml42=sd(b8) 
        
        
        b9=numeric()
        for(i in 1:1000){
          a=numeric()
          for(j in 1:50){
            a[j]=all_snps_H_no_y[[j]]$p_after3[i]
          }
          b9[i]=(sum(a<0.05/50)>0)
        }
        cml43=mean(b9)
        sdcml43=sd(b9)
        
        
        
        
        b10=numeric()
        varb10=numeric()
        for(i in 1:length(all_snps_y_only)){
          b10[i]=mean(all_snps_y_only[[i]]$p_after1<0.05)
          varb10[i]= var(all_snps_y_only[[i]]$p_after1<0.05)
        }
        cml51=mean(b10)
        sdcml51=sqrt(mean(varb10))
        
        b11t=numeric()
        varb11t=numeric()
        for(i in 1:length(all_snps_y_only)){
          b11t[i]=mean(all_snps_y_only[[i]]$p_after2<0.05)
          varb11t[i]=var(all_snps_y_only[[i]]$p_after2<0.05)
        }
        cml52=mean(b11t)
        sdcml52=sqrt(mean(varb11t))
        
        b12t=numeric()
        varb12t=numeric()
        for(i in 1:length(all_snps_y_only)){
          b12t[i]=mean(all_snps_y_only[[i]]$p_after3<0.05)
          varb12t[i]=var(all_snps_y_only[[i]]$p_after3<0.05)
        }
        cml53=mean(b12t)
        sdcml53=sqrt(mean(varb12t))
        
        
        
        b13t=numeric()
        varb13t=numeric()
        for(i in 1:length(all_snps_H_and_y)){
          b13t[i]=mean(all_snps_H_and_y[[i]]$p_after1<0.05)
          varb13t[i]=var(all_snps_H_and_y[[i]]$p_after1<0.05)
        }
        cml61=mean(b13t)
        sdcml61=sqrt(mean(varb13t))
        
        
        b14t=numeric()
        varb14t=numeric()
        for(i in 1:length(all_snps_H_and_y)){
          b14t[i]=mean(all_snps_H_and_y[[i]]$p_after2<0.05)
          varb14t[i]=var(all_snps_H_and_y[[i]]$p_after2<0.05)
        }
        cml62=mean(b14t)
        sdcml62=sqrt(mean(varb14t))
        
        b15t=numeric()
        varb15t=numeric()
        for(i in 1:length(all_snps_H_and_y)){
          b15t[i]=mean(all_snps_H_and_y[[i]]$p_after3<0.05)
          varb15t[i]=var(all_snps_H_and_y[[i]]$p_after3<0.05)
        }
        cml63=mean(b15t)
        sdcml63=sqrt(mean(varb15t))
        
        
        
        power1=numeric()
        power21=numeric()
        power22=numeric()
        power23=numeric()
        
        
        for(i in 1:length(all_snps_y)){
          power1[i]=mean(all_snps_y[[i]]$p_before<0.05)
          power21[i]=mean(all_snps_y[[i]]$p_after1<0.05)
          power22[i]=mean(all_snps_y[[i]]$p_after2<0.05)
          power23[i]=mean(all_snps_y[[i]]$p_after3<0.05)
        }
        change1=power21-power1
        change2=power22-power1
        change3=power23-power1
        if(sum(change1>0)>0){
          cml71=power21[which.max(change1)]
          sdcml71=sd(all_snps_y[[which.max(change1)]]$p_after1<0.05)
        }
        if(sum(change1>0)==0){
          cml71=NA
          sdcml71=NA
        }
        if(sum(change2>0)>0){
          cml72=power22[which.max(change2)]
          sdcml72=sd(all_snps_y[[which.max(change2)]]$p_after2<0.05)
        }
        if(sum(change2>0)==0){
          cml72=NA
          sdcml72=NA
        }
        if(sum(change3>0)>0){
          cml73=power23[which.max(change3)]
          sdcml73=sd(all_snps_y[[which.max(change3)]]$p_after3<0.05)
        }
        if(sum(change3>0)==0){
          cml73=NA
          sdcml73=NA
        }
        
        
        
        if(sum(change1<0)>0){
          cml81=power21[which.min(change1)]
          sdcml81=sd(all_snps_y[[which.min(change1)]]$p_after1<0.05)
        }
        if(sum(change1<0)==0){
          cml81=NA
          sdcml81=NA
        }
        
        if(sum(change2<0)>0){
          cml82=power22[which.min(change2)]
          sdcml82=sd(all_snps_y[[which.min(change2)]]$p_after2<0.05)
        }
        if(sum(change2<0)==0){
          cml82=NA
          sdcml82=NA
        }
        
        
        if(sum(change3<0)>0){
          cml83=power23[which.min(change3)]
          sdcml83=sd(all_snps_y[[which.min(change3)]]$p_after3<0.05)
        }
        if(sum(change3<0)==0){
          cml83=NA
          sdcml83=NA
        }
      }
      
      if(10>0){
        b1=numeric()
        varb1=numeric()
        for(i in 1: length(null_snps)){
          b1[i]=mean(null_snps[[i]]$p_egger1<0.05)
          varb1[i]=var(null_snps[[i]]$p_egger1<0.05)
        }
        egger11=mean(b1)
        sdegger11=sqrt(mean(varb1))
        
        
        b2=numeric()
        varb2=numeric()
        for(i in 1: length(null_snps)){
          b2[i]=mean(null_snps[[i]]$p_egger2<0.05)
          varb2[i]=var(null_snps[[i]]$p_egger2<0.05)
        }
        egger12=mean(b2)
        sdegger12=sqrt(mean(varb2))
        
        
        
        
        b4=numeric()
        varb4=numeric()
        for(i in 1: length(all_snps_H_no_y)){
          b4[i]=mean(all_snps_H_no_y[[i]]$p_egger1<0.05)
          varb4[i]=var(all_snps_H_no_y[[i]]$p_egger1<0.05)
        }
        egger21=mean(b4)
        sdegger21=sqrt(mean(varb4))
        
        
        b5=numeric()
        varb5=numeric()
        for(i in 1: length(all_snps_H_no_y)){
          b5[i]=mean(all_snps_H_no_y[[i]]$p_egger2<0.05)
          varb5[i]=var(all_snps_H_no_y[[i]]$p_egger2<0.05)
        }
        egger22=mean(b5)
        sdegger22=sqrt(mean(varb5))
        
        
        
        
        
        a=numeric()
        b=numeric()
        for(i in 1:length(all_snps_no_y)){
          a[i]=mean(all_snps_no_y[[i]]$p_before<0.05)
          b[i]=mean(all_snps_no_y[[i]]$p_egger1<0.05)
        }
        egger31=b[which.max(a)]
        sdegger31=sd(all_snps_no_y[[which.max(a)]]$p_egger1<0.05)
        
        a=numeric()
        b=numeric()
        for(i in 1:length(all_snps_no_y)){
          a[i]=mean(all_snps_no_y[[i]]$p_before<0.05)
          b[i]=mean(all_snps_no_y[[i]]$p_egger2<0.05)
        }
        egger32=b[which.max(a)]
        sdegger32=sd(all_snps_no_y[[which.max(a)]]$p_egger2<0.05)
        
        
        
        
        
        b7=numeric()
        for(i in 1:1000){
          a=numeric()
          for(j in 1:50){
            a[j]=all_snps_H_no_y[[j]]$p_egger1[i]
          }
          b7[i]=(sum(a<0.05/50)>0)
        }
        egger41=mean(b7)
        sdegger41=sd(b7) 
        
        b8=numeric()
        for(i in 1:1000){
          a=numeric()
          for(j in 1:50){
            a[j]=all_snps_H_no_y[[j]]$p_egger2[i]
          }
          b8[i]=(sum(a<0.05/50)>0)
        }
        egger42=mean(b8)
        sdegger42=sd(b8) 
        
        
        
        
        
        
        b10=numeric()
        varb10=numeric()
        for(i in 1:length(all_snps_y_only)){
          b10[i]=mean(all_snps_y_only[[i]]$p_egger1<0.05)
          varb10[i]=var(all_snps_y_only[[i]]$p_egger1<0.05)
        }
        egger51=mean(b10)
        sdegger51=sqrt(mean(varb10))
        
        b11t=numeric()
        varb11t=numeric()
        for(i in 1:length(all_snps_y_only)){
          b11t[i]=mean(all_snps_y_only[[i]]$p_egger2<0.05)
          varb11t[i]=var(all_snps_y_only[[i]]$p_egger2<0.05)
        }
        egger52=mean(b11t)
        sdegger52=sqrt(mean(varb11t))
        
        
        
        
        
        b13t=numeric()
        varb13t=numeric()
        for(i in 1:length(all_snps_H_and_y)){
          b13t[i]=mean(all_snps_H_and_y[[i]]$p_egger1<0.05)
          varb13t[i]=var(all_snps_H_and_y[[i]]$p_egger1<0.05)
        }
        egger61=mean(b13t)
        sdegger61=sqrt(mean(varb13t))
        
        
        b14t=numeric()
        varb14t=numeric()
        for(i in 1:length(all_snps_H_and_y)){
          b14t[i]=mean(all_snps_H_and_y[[i]]$p_egger2<0.05)
          varb14t[i]=var(all_snps_H_and_y[[i]]$p_egger2<0.05)
        }
        egger62=mean(b14t)
        sdegger62=sqrt(mean(varb14t))
        
        
        
        
        power1=numeric()
        power21=numeric()
        power22=numeric()
        power_egger21=numeric()
        power_egger22=numeric()
        
        for(i in 1:length(all_snps_y)){
          power1[i]=mean(all_snps_y[[i]]$p_before<0.05)
          power21[i]=mean(all_snps_y[[i]]$p_after1<0.05)
          power22[i]=mean(all_snps_y[[i]]$p_after2<0.05)
          power_egger21[i]=mean(all_snps_y[[i]]$p_egger1<0.05)
          power_egger22[i]=mean(all_snps_y[[i]]$p_egger2<0.05)
        }
        
        change1=power21-power1
        change2=power22-power1
        
        if(sum(change1>0)>0){
          egger71=power_egger21[which.max(change1)]
          sdegger71=sd(all_snps_y[[which.max(change1)]]$p_egger1<0.05)
        }
        if(sum(change1>0)==0){
          egger71=NA
          sdegger71=NA
        }
        if(sum(change2>0)>0){
          egger72=power_egger22[which.max(change2)]
          sdegger72=sd(all_snps_y[[which.max(change2)]]$p_egger2<0.05)
        }
        if(sum(change2>0)==0){
          egger72=NA
          sdegger72=NA
        }
        
        
        if(sum(change1<0)>0){
          egger81=power_egger21[which.min(change1)]
          sdegger81=sd(all_snps_y[[which.min(change1)]]$p_egger1<0.05)
        }
        if(sum(change1<0)==0){
          egger81=NA
          sdegger81=NA
        }
        
        if(sum(change2<0)>0){
          egger82=power_egger22[which.min(change2)]
          sdegger82=sd(all_snps_y[[which.min(change2)]]$p_egger2<0.05)
        }
        if(sum(change2<0)==0){
          egger82=NA
          sdegger82=NA
        }
        
      }
      
      if(10>0){
        b1=numeric()
        varb1=numeric()
        for(i in 1: length(null_snps)){
          b1[i]=mean(null_snps[[i]]$p_DD1<0.05)
          varb1[i]=var(null_snps[[i]]$p_DD1<0.05)
        }
        DD11=mean(b1)
        sdDD11=sqrt(mean(varb1))
        
        
        b2=numeric()
        varb2=numeric()
        for(i in 1: length(null_snps)){
          b2[i]=mean(null_snps[[i]]$p_DD2<0.05)
          varb2[i]=var(null_snps[[i]]$p_DD2<0.05)
        }
        DD12=mean(b2)
        sdDD12=sqrt(mean(varb2))
        
        
        
        
        b4=numeric()
        varb4=numeric()
        for(i in 1: length(all_snps_H_no_y)){
          b4[i]=mean(all_snps_H_no_y[[i]]$p_DD1<0.05)
          varb4[i]=var(all_snps_H_no_y[[i]]$p_DD1<0.05)
        }
        DD21=mean(b4)
        sdDD21=sqrt(mean(varb4))
        
        
        b5=numeric()
        varb5=numeric()
        for(i in 1: length(all_snps_H_no_y)){
          b5[i]=mean(all_snps_H_no_y[[i]]$p_DD2<0.05)
          varb5[i]=var(all_snps_H_no_y[[i]]$p_DD2<0.05)
        }
        DD22=mean(b5)
        sdDD22=sqrt(mean(varb5))
        
        
        
        
        
        a=numeric()
        b=numeric()
        for(i in 1:length(all_snps_no_y)){
          a[i]=mean(all_snps_no_y[[i]]$p_before<0.05)
          b[i]=mean(all_snps_no_y[[i]]$p_DD1<0.05)
        }
        DD31=b[which.max(a)]
        sdDD31=sd(all_snps_no_y[[which.max(a)]]$p_DD1<0.05)
        
        a=numeric()
        b=numeric()
        for(i in 1:length(all_snps_no_y)){
          a[i]=mean(all_snps_no_y[[i]]$p_before<0.05)
          b[i]=mean(all_snps_no_y[[i]]$p_DD2<0.05)
        }
        DD32=b[which.max(a)]
        sdDD32=sd(all_snps_no_y[[which.max(a)]]$p_DD2<0.05)
        
        
        
        
        
        b7=numeric()
        for(i in 1:1000){
          a=numeric()
          for(j in 1:50){
            a[j]=all_snps_H_no_y[[j]]$p_DD1[i]
          }
          b7[i]=(sum(a<0.05/50)>0)
        }
        DD41=mean(b7)
        sdDD41=sd(b7) 
        
        b8=numeric()
        for(i in 1:1000){
          a=numeric()
          for(j in 1:50){
            a[j]=all_snps_H_no_y[[j]]$p_DD2[i]
          }
          b8[i]=(sum(a<0.05/50)>0)
        }
        DD42=mean(b8)
        sdDD42=sd(b8) 
        
        
        
        
        
        
        b10=numeric()
        varb10=numeric()
        for(i in 1:length(all_snps_y_only)){
          b10[i]=mean(all_snps_y_only[[i]]$p_DD1<0.05)
          varb10[i]=var(all_snps_y_only[[i]]$p_DD1<0.05)
        }
        DD51=mean(b10)
        sdDD51=sqrt(mean(varb10))
        
        b11t=numeric()
        varb11t=numeric()
        for(i in 1:length(all_snps_y_only)){
          b11t[i]=mean(all_snps_y_only[[i]]$p_DD2<0.05)
          varb11t[i]=var(all_snps_y_only[[i]]$p_DD2<0.05)
        }
        DD52=mean(b11t)
        sdDD52=sqrt(mean(varb11t))
        
        
        
        
        
        b13t=numeric()
        varb13t=numeric()
        for(i in 1:length(all_snps_H_and_y)){
          b13t[i]=mean(all_snps_H_and_y[[i]]$p_DD1<0.05)
          varb13t[i]=var(all_snps_H_and_y[[i]]$p_DD1<0.05)
        }
        DD61=mean(b13t)
        sdDD61=sqrt(mean(varb13t))
        
        
        b14t=numeric()
        varb14t=numeric()
        for(i in 1:length(all_snps_H_and_y)){
          b14t[i]=mean(all_snps_H_and_y[[i]]$p_DD2<0.05)
          varb14t[i]=var(all_snps_H_and_y[[i]]$p_DD2<0.05)
        }
        DD62=mean(b14t)
        sdDD62=sqrt(mean(varb14t))
        
        
        
        
        power1=numeric()
        power21=numeric()
        power22=numeric()
        power_DD21=numeric()
        power_DD22=numeric()
        
        for(i in 1:length(all_snps_y)){
          power1[i]=mean(all_snps_y[[i]]$p_before<0.05)
          power21[i]=mean(all_snps_y[[i]]$p_after1<0.05)
          power22[i]=mean(all_snps_y[[i]]$p_after2<0.05)
          power_DD21[i]=mean(all_snps_y[[i]]$p_DD1<0.05)
          power_DD22[i]=mean(all_snps_y[[i]]$p_DD2<0.05)
        }
        
        change1=power21-power1
        change2=power22-power1
        
        if(sum(change1>0)>0){
          DD71=power_DD21[which.max(change1)]
          sdDD71=sd(all_snps_y[[which.max(change1)]]$p_DD1<0.05)
        }
        if(sum(change1>0)==0){
          DD71=NA
          sdDD71=NA
        }
        if(sum(change2>0)>0){
          DD72=power_DD22[which.max(change2)]
          sdDD72=sd(all_snps_y[[which.max(change2)]]$p_DD2<0.05)
        }
        if(sum(change2>0)==0){
          DD72=NA
          sdDD72=NA
        }
        
        
        if(sum(change1<0)>0){
          DD81=power_DD21[which.min(change1)]
          sdDD81=sd(all_snps_y[[which.min(change1)]]$p_DD1<0.05)
        }
        if(sum(change1<0)==0){
          DD81=NA
          sdDD81=NA
        }
        
        if(sum(change2<0)>0){
          DD82=power_DD22[which.min(change2)]
          sdDD82=sd(all_snps_y[[which.min(change2)]]$p_DD2<0.05)
        }
        if(sum(change2<0)==0){
          DD82=NA
          sdDD82=NA
        }
        
      }
      
      if(10>0){
        b1=numeric()
        varb1=numeric()
        for(i in 1: length(null_snps)){
          b1[i]=mean(null_snps[[i]]$p_SH1<0.05)
          varb1[i]=var(null_snps[[i]]$p_SH1<0.05)
        }
        SH11=mean(b1)
        sdSH11=sqrt(mean(varb1))
        
        
        b2=numeric()
        varb2=numeric()
        for(i in 1: length(null_snps)){
          b2[i]=mean(null_snps[[i]]$p_SH2<0.05)
          varb2[i]=var(null_snps[[i]]$p_SH2<0.05)
        }
        SH12=mean(b2)
        sdSH12=sqrt(mean(varb2))
        
        
        
        
        b4=numeric()
        varb4=numeric()
        for(i in 1: length(all_snps_H_no_y)){
          b4[i]=mean(all_snps_H_no_y[[i]]$p_SH1<0.05)
          varb4[i]=var(all_snps_H_no_y[[i]]$p_SH1<0.05)
        }
        SH21=mean(b4)
        sdSH21=sqrt(mean(varb4))
        
        
        b5=numeric()
        varb5=numeric()
        for(i in 1: length(all_snps_H_no_y)){
          b5[i]=mean(all_snps_H_no_y[[i]]$p_SH2<0.05)
          varb5[i]=var(all_snps_H_no_y[[i]]$p_SH2<0.05)
        }
        SH22=mean(b5)
        sdSH22=sqrt(mean(varb5))
        
        
        
        
        
        a=numeric()
        b=numeric()
        for(i in 1:length(all_snps_no_y)){
          a[i]=mean(all_snps_no_y[[i]]$p_before<0.05)
          b[i]=mean(all_snps_no_y[[i]]$p_SH1<0.05)
        }
        SH31=b[which.max(a)]
        sdSH31=sd(all_snps_no_y[[which.max(a)]]$p_SH1<0.05)
        
        a=numeric()
        b=numeric()
        for(i in 1:length(all_snps_no_y)){
          a[i]=mean(all_snps_no_y[[i]]$p_before<0.05)
          b[i]=mean(all_snps_no_y[[i]]$p_SH2<0.05)
        }
        SH32=b[which.max(a)]
        sdSH32=sd(all_snps_no_y[[which.max(a)]]$p_SH2<0.05)
        
        
        
        
        
        b7=numeric()
        for(i in 1:1000){
          a=numeric()
          for(j in 1:50){
            a[j]=all_snps_H_no_y[[j]]$p_SH1[i]
          }
          b7[i]=(sum(a<0.05/50)>0)
        }
        SH41=mean(b7)
        sdSH41=sd(b7) 
        
        b8=numeric()
        for(i in 1:1000){
          a=numeric()
          for(j in 1:50){
            a[j]=all_snps_H_no_y[[j]]$p_SH2[i]
          }
          b8[i]=(sum(a<0.05/50)>0)
        }
        SH42=mean(b8)
        sdSH42=sd(b8) 
        
        
        
        
        
        
        b10=numeric()
        varb10=numeric()
        for(i in 1:length(all_snps_y_only)){
          b10[i]=mean(all_snps_y_only[[i]]$p_SH1<0.05)
          varb10[i]=var(all_snps_y_only[[i]]$p_SH1<0.05)
        }
        SH51=mean(b10)
        sdSH51=sqrt(mean(varb10))
        
        b11t=numeric()
        varb11t=numeric()
        for(i in 1:length(all_snps_y_only)){
          b11t[i]=mean(all_snps_y_only[[i]]$p_SH2<0.05)
          varb11t[i]=var(all_snps_y_only[[i]]$p_SH2<0.05)
        }
        SH52=mean(b11t)
        sdSH52=sqrt(mean(varb11t))
        
        
        
        
        
        b13t=numeric()
        varb13t=numeric()
        for(i in 1:length(all_snps_H_and_y)){
          b13t[i]=mean(all_snps_H_and_y[[i]]$p_SH1<0.05)
          varb13t[i]=var(all_snps_H_and_y[[i]]$p_SH1<0.05)
        }
        SH61=mean(b13t)
        sdSH61=sqrt(mean(varb13t))
        
        
        b14t=numeric()
        varb14t=numeric()
        for(i in 1:length(all_snps_H_and_y)){
          b14t[i]=mean(all_snps_H_and_y[[i]]$p_SH2<0.05)
          varb14t[i]=var(all_snps_H_and_y[[i]]$p_SH2<0.05)
        }
        SH62=mean(b14t)
        sdSH62=sqrt(mean(varb14t))
        
        
        
        
        power1=numeric()
        power21=numeric()
        power22=numeric()
        power_SH21=numeric()
        power_SH22=numeric()
        
        for(i in 1:length(all_snps_y)){
          power1[i]=mean(all_snps_y[[i]]$p_before<0.05)
          power21[i]=mean(all_snps_y[[i]]$p_after1<0.05)
          power22[i]=mean(all_snps_y[[i]]$p_after2<0.05)
          power_SH21[i]=mean(all_snps_y[[i]]$p_SH1<0.05)
          power_SH22[i]=mean(all_snps_y[[i]]$p_SH2<0.05)
        }
        
        change1=power21-power1
        change2=power22-power1
        
        if(sum(change1>0)>0){
          SH71=power_SH21[which.max(change1)]
          sdSH71=sd(all_snps_y[[which.max(change1)]]$p_SH1<0.05)
        }
        if(sum(change1>0)==0){
          SH71=NA
          sdSH71=NA
        }
        if(sum(change2>0)>0){
          SH72=power_SH22[which.max(change2)]
          sdSH72=sd(all_snps_y[[which.max(change2)]]$p_SH2<0.05)
        }
        if(sum(change2>0)==0){
          SH72=NA
          sdSH72=NA
        }
        
        
        if(sum(change1<0)>0){
          SH81=power_SH21[which.min(change1)]
          sdSH81=sd(all_snps_y[[which.min(change1)]]$p_SH1<0.05)
        }
        if(sum(change1<0)==0){
          SH81=NA
          sdSH81=NA
        }
        
        if(sum(change2<0)>0){
          SH82=power_SH22[which.min(change2)]
          sdSH82=sd(all_snps_y[[which.min(change2)]]$p_SH2<0.05)
        }
        if(sum(change2<0)==0){
          SH82=NA
          sdSH82=NA
        }
        
      }
      
      if(10>0){
        b1=numeric()
        varb1=numeric()
        for(i in 1: length(null_snps)){
          b1[i]=mean(null_snps[[i]]$p_median1<0.05)
          varb1[i]=var(null_snps[[i]]$p_median1<0.05)
        }
        median11=mean(b1)
        sdmedian11=sqrt(mean(varb1))
        
        
        b2=numeric()
        varb2=numeric()
        for(i in 1: length(null_snps)){
          b2[i]=mean(null_snps[[i]]$p_median2<0.05)
          varb2[i]=var(null_snps[[i]]$p_median2<0.05)
        }
        median12=mean(b2)
        sdmedian12=sqrt(mean(varb2))
        
        
        
        
        b4=numeric()
        varb4=numeric()
        for(i in 1: length(all_snps_H_no_y)){
          b4[i]=mean(all_snps_H_no_y[[i]]$p_median1<0.05)
          varb4[i]=var(all_snps_H_no_y[[i]]$p_median1<0.05)
        }
        median21=mean(b4)
        sdmedian21=sqrt(mean(varb4))
        
        
        b5=numeric()
        varb5=numeric()
        for(i in 1: length(all_snps_H_no_y)){
          b5[i]=mean(all_snps_H_no_y[[i]]$p_median2<0.05)
          varb5[i]=var(all_snps_H_no_y[[i]]$p_median2<0.05)
        }
        median22=mean(b5)
        sdmedian22=sqrt(mean(varb5))
        
        
        
        
        
        a=numeric()
        b=numeric()
        for(i in 1:length(all_snps_no_y)){
          a[i]=mean(all_snps_no_y[[i]]$p_before<0.05)
          b[i]=mean(all_snps_no_y[[i]]$p_median1<0.05)
        }
        median31=b[which.max(a)]
        sdmedian31=sd(all_snps_no_y[[which.max(a)]]$p_median1<0.05)
        
        a=numeric()
        b=numeric()
        for(i in 1:length(all_snps_no_y)){
          a[i]=mean(all_snps_no_y[[i]]$p_before<0.05)
          b[i]=mean(all_snps_no_y[[i]]$p_median2<0.05)
        }
        median32=b[which.max(a)]
        sdmedian32=sd(all_snps_no_y[[which.max(a)]]$p_median2<0.05)
        
        
        
        
        
        b7=numeric()
        for(i in 1:1000){
          a=numeric()
          for(j in 1:50){
            a[j]=all_snps_H_no_y[[j]]$p_median1[i]
          }
          b7[i]=(sum(a<0.05/50)>0)
        }
        median41=mean(b7)
        sdmedian41=sd(b7) 
        
        b8=numeric()
        for(i in 1:1000){
          a=numeric()
          for(j in 1:50){
            a[j]=all_snps_H_no_y[[j]]$p_median2[i]
          }
          b8[i]=(sum(a<0.05/50)>0)
        }
        median42=mean(b8)
        sdmedian42=sd(b8) 
        
        
        
        
        
        
        b10=numeric()
        varb10=numeric()
        for(i in 1:length(all_snps_y_only)){
          b10[i]=mean(all_snps_y_only[[i]]$p_median1<0.05)
          varb10[i]=var(all_snps_y_only[[i]]$p_median1<0.05)
        }
        median51=mean(b10)
        sdmedian51=sqrt(mean(varb10))
        
        b11t=numeric()
        varb11t=numeric()
        for(i in 1:length(all_snps_y_only)){
          b11t[i]=mean(all_snps_y_only[[i]]$p_median2<0.05)
          varb11t[i]=var(all_snps_y_only[[i]]$p_median2<0.05)
        }
        median52=mean(b11t)
        sdmedian52=sqrt(mean(varb11t))
        
        
        
        
        
        b13t=numeric()
        varb13t=numeric()
        for(i in 1:length(all_snps_H_and_y)){
          b13t[i]=mean(all_snps_H_and_y[[i]]$p_median1<0.05)
          varb13t[i]=var(all_snps_H_and_y[[i]]$p_median1<0.05)
        }
        median61=mean(b13t)
        sdmedian61=sqrt(mean(varb13t))
        
        
        b14t=numeric()
        varb14t=numeric()
        for(i in 1:length(all_snps_H_and_y)){
          b14t[i]=mean(all_snps_H_and_y[[i]]$p_median2<0.05)
          varb14t[i]=var(all_snps_H_and_y[[i]]$p_median2<0.05)
        }
        median62=mean(b14t)
        sdmedian62=sqrt(mean(varb14t))
        
        
        
        
        power1=numeric()
        power21=numeric()
        power22=numeric()
        power_median21=numeric()
        power_median22=numeric()
        
        for(i in 1:length(all_snps_y)){
          power1[i]=mean(all_snps_y[[i]]$p_before<0.05)
          power21[i]=mean(all_snps_y[[i]]$p_after1<0.05)
          power22[i]=mean(all_snps_y[[i]]$p_after2<0.05)
          power_median21[i]=mean(all_snps_y[[i]]$p_median1<0.05)
          power_median22[i]=mean(all_snps_y[[i]]$p_median2<0.05)
        }
        
        change1=power21-power1
        change2=power22-power1
        
        if(sum(change1>0)>0){
          median71=power_median21[which.max(change1)]
          sdmedian71=sd(all_snps_y[[which.max(change1)]]$p_median1<0.05)
        }
        if(sum(change1>0)==0){
          median71=NA
          sdmedian71=NA
        }
        if(sum(change2>0)>0){
          median72=power_median22[which.max(change2)]
          sdmedian72=sd(all_snps_y[[which.max(change2)]]$p_median2<0.05)
        }
        if(sum(change2>0)==0){
          median72=NA
          sdmedian72=NA
        }
        
        
        if(sum(change1<0)>0){
          median81=power_median21[which.min(change1)]
          sdmedian81=sd(all_snps_y[[which.min(change1)]]$p_median1<0.05)
        }
        if(sum(change1<0)==0){
          median81=NA
          sdmedian81=NA
        }
        
        if(sum(change2<0)>0){
          median82=power_median22[which.min(change2)]
          sdmedian82=sd(all_snps_y[[which.min(change2)]]$p_median2<0.05)
        }
        if(sum(change2<0)==0){
          median82=NA
          sdmedian82=NA
        }
        
      }
      
      if(10>0){
        b1=numeric()
        varb1=numeric()
        for(i in 1: length(null_snps)){
          b1[i]=mean(null_snps[[i]]$p_ivw1<0.05)
          varb1[i]=var(null_snps[[i]]$p_ivw1<0.05)
        }
        ivw11=mean(b1)
        sdivw11=sqrt(mean(varb1))
        
        
        b2=numeric()
        varb2=numeric()
        for(i in 1: length(null_snps)){
          b2[i]=mean(null_snps[[i]]$p_ivw2<0.05)
          varb2[i]=var(null_snps[[i]]$p_ivw2<0.05)
        }
        ivw12=mean(b2)
        sdivw12=sqrt(mean(varb2))
        
        
        
        
        b4=numeric()
        varb4=numeric()
        for(i in 1: length(all_snps_H_no_y)){
          b4[i]=mean(all_snps_H_no_y[[i]]$p_ivw1<0.05)
          varb4[i]=var(all_snps_H_no_y[[i]]$p_ivw1<0.05)
        }
        ivw21=mean(b4)
        sdivw21=sqrt(mean(varb4))
        
        
        b5=numeric()
        varb5=numeric()
        for(i in 1: length(all_snps_H_no_y)){
          b5[i]=mean(all_snps_H_no_y[[i]]$p_ivw2<0.05)
          varb5[i]=var(all_snps_H_no_y[[i]]$p_ivw2<0.05)
        }
        ivw22=mean(b5)
        sdivw22=sqrt(mean(varb5))
        
        
        
        
        
        a=numeric()
        b=numeric()
        for(i in 1:length(all_snps_no_y)){
          a[i]=mean(all_snps_no_y[[i]]$p_before<0.05)
          b[i]=mean(all_snps_no_y[[i]]$p_ivw1<0.05)
        }
        ivw31=b[which.max(a)]
        sdivw31=sd(all_snps_no_y[[which.max(a)]]$p_ivw1<0.05)
        
        a=numeric()
        b=numeric()
        for(i in 1:length(all_snps_no_y)){
          a[i]=mean(all_snps_no_y[[i]]$p_before<0.05)
          b[i]=mean(all_snps_no_y[[i]]$p_ivw2<0.05)
        }
        ivw32=b[which.max(a)]
        sdivw32=sd(all_snps_no_y[[which.max(a)]]$p_ivw2<0.05)
        
        
        
        
        
        b7=numeric()
        for(i in 1:1000){
          a=numeric()
          for(j in 1:50){
            a[j]=all_snps_H_no_y[[j]]$p_ivw1[i]
          }
          b7[i]=(sum(a<0.05/50)>0)
        }
        ivw41=mean(b7)
        sdivw41=sd(b7) 
        
        b8=numeric()
        for(i in 1:1000){
          a=numeric()
          for(j in 1:50){
            a[j]=all_snps_H_no_y[[j]]$p_ivw2[i]
          }
          b8[i]=(sum(a<0.05/50)>0)
        }
        ivw42=mean(b8)
        sdivw42=sd(b8) 
        
        
        
        
        
        
        b10=numeric()
        varb10=numeric()
        for(i in 1:length(all_snps_y_only)){
          b10[i]=mean(all_snps_y_only[[i]]$p_ivw1<0.05)
          varb10[i]=var(all_snps_y_only[[i]]$p_ivw1<0.05)
        }
        ivw51=mean(b10)
        sdivw51=sqrt(mean(varb10))
        
        b11t=numeric()
        varb11t=numeric()
        for(i in 1:length(all_snps_y_only)){
          b11t[i]=mean(all_snps_y_only[[i]]$p_ivw2<0.05)
          varb11t[i]=var(all_snps_y_only[[i]]$p_ivw2<0.05)
        }
        ivw52=mean(b11t)
        sdivw52=sqrt(mean(varb11t))
        
        
        
        
        
        b13t=numeric()
        varb13t=numeric()
        for(i in 1:length(all_snps_H_and_y)){
          b13t[i]=mean(all_snps_H_and_y[[i]]$p_ivw1<0.05)
          varb13t[i]=var(all_snps_H_and_y[[i]]$p_ivw1<0.05)
        }
        ivw61=mean(b13t)
        sdivw61=sqrt(mean(varb13t))
        
        
        b14t=numeric()
        varb14t=numeric()
        for(i in 1:length(all_snps_H_and_y)){
          b14t[i]=mean(all_snps_H_and_y[[i]]$p_ivw2<0.05)
          varb14t[i]=var(all_snps_H_and_y[[i]]$p_ivw2<0.05)
        }
        ivw62=mean(b14t)
        sdivw62=sqrt(mean(varb14t))
        
        
        
        
        power1=numeric()
        power21=numeric()
        power22=numeric()
        power_ivw21=numeric()
        power_ivw22=numeric()
        
        for(i in 1:length(all_snps_y)){
          power1[i]=mean(all_snps_y[[i]]$p_before<0.05)
          power21[i]=mean(all_snps_y[[i]]$p_after1<0.05)
          power22[i]=mean(all_snps_y[[i]]$p_after2<0.05)
          power_ivw21[i]=mean(all_snps_y[[i]]$p_ivw1<0.05)
          power_ivw22[i]=mean(all_snps_y[[i]]$p_ivw2<0.05)
        }
        
        change1=power21-power1
        change2=power22-power1
        
        if(sum(change1>0)>0){
          ivw71=power_ivw21[which.max(change1)]
          sdivw71=sd(all_snps_y[[which.max(change1)]]$p_ivw1<0.05)
        }
        if(sum(change1>0)==0){
          ivw71=NA
          sdivw71=NA
        }
        if(sum(change2>0)>0){
          ivw72=power_ivw22[which.max(change2)]
          sdivw72=sd(all_snps_y[[which.max(change2)]]$p_ivw2<0.05)
        }
        if(sum(change2>0)==0){
          ivw72=NA
          sdivw72=NA
        }
        
        
        if(sum(change1<0)>0){
          ivw81=power_ivw21[which.min(change1)]
          sdivw81=sd(all_snps_y[[which.min(change1)]]$p_ivw1<0.05)
        }
        if(sum(change1<0)==0){
          ivw81=NA
          sdivw81=NA
        }
        
        if(sum(change2<0)>0){
          ivw82=power_ivw22[which.min(change2)]
          sdivw82=sd(all_snps_y[[which.min(change2)]]$p_ivw2<0.05)
        }
        if(sum(change2<0)==0){
          ivw82=NA
          sdivw82=NA
        }
        
      }
      
      if(10>0){
        b1=numeric()
        varb1=numeric()
        for(i in 1: length(null_snps)){
          b1[i]=mean(null_snps[[i]]$p_lasso1<0.05)
          varb1[i]=var(null_snps[[i]]$p_lasso1<0.05)
        }
        lasso11=mean(b1)
        sdlasso11=sqrt(mean(varb1))
        
        
        b2=numeric()
        varb2=numeric()
        for(i in 1: length(null_snps)){
          b2[i]=mean(null_snps[[i]]$p_lasso2<0.05)
          varb2[i]=var(null_snps[[i]]$p_lasso2<0.05)
        }
        lasso12=mean(b2)
        sdlasso12=sqrt(mean(varb2))
        
        
        
        
        b4=numeric()
        varb4=numeric()
        for(i in 1: length(all_snps_H_no_y)){
          b4[i]=mean(all_snps_H_no_y[[i]]$p_lasso1<0.05)
          varb4[i]=var(all_snps_H_no_y[[i]]$p_lasso1<0.05)
        }
        lasso21=mean(b4)
        sdlasso21=sqrt(mean(varb4))
        
        
        b5=numeric()
        varb5=numeric()
        for(i in 1: length(all_snps_H_no_y)){
          b5[i]=mean(all_snps_H_no_y[[i]]$p_lasso2<0.05)
          varb5[i]=var(all_snps_H_no_y[[i]]$p_lasso2<0.05)
        }
        lasso22=mean(b5)
        sdlasso22=sqrt(mean(varb5))
        
        
        
        
        
        a=numeric()
        b=numeric()
        for(i in 1:length(all_snps_no_y)){
          a[i]=mean(all_snps_no_y[[i]]$p_before<0.05)
          b[i]=mean(all_snps_no_y[[i]]$p_lasso1<0.05)
        }
        lasso31=b[which.max(a)]
        sdlasso31=sd(all_snps_no_y[[which.max(a)]]$p_lasso1<0.05)
        
        a=numeric()
        b=numeric()
        for(i in 1:length(all_snps_no_y)){
          a[i]=mean(all_snps_no_y[[i]]$p_before<0.05)
          b[i]=mean(all_snps_no_y[[i]]$p_lasso2<0.05)
        }
        lasso32=b[which.max(a)]
        sdlasso32=sd(all_snps_no_y[[which.max(a)]]$p_lasso2<0.05)
        
        
        
        
        
        b7=numeric()
        for(i in 1:1000){
          a=numeric()
          for(j in 1:50){
            a[j]=all_snps_H_no_y[[j]]$p_lasso1[i]
          }
          b7[i]=(sum(a<0.05/50)>0)
        }
        lasso41=mean(b7)
        sdlasso41=sd(b7) 
        
        b8=numeric()
        for(i in 1:1000){
          a=numeric()
          for(j in 1:50){
            a[j]=all_snps_H_no_y[[j]]$p_lasso2[i]
          }
          b8[i]=(sum(a<0.05/50)>0)
        }
        lasso42=mean(b8)
        sdlasso42=sd(b8) 
        
        
        
        
        
        
        b10=numeric()
        varb10=numeric()
        for(i in 1:length(all_snps_y_only)){
          b10[i]=mean(all_snps_y_only[[i]]$p_lasso1<0.05)
          varb10[i]=var(all_snps_y_only[[i]]$p_lasso1<0.05)
        }
        lasso51=mean(b10)
        sdlasso51=sqrt(mean(varb10))
        
        b11t=numeric()
        varb11t=numeric()
        for(i in 1:length(all_snps_y_only)){
          b11t[i]=mean(all_snps_y_only[[i]]$p_lasso2<0.05)
          varb11t[i]=var(all_snps_y_only[[i]]$p_lasso2<0.05)
        }
        lasso52=mean(b11t)
        sdlasso52=sqrt(mean(varb11t))
        
        
        
        
        
        b13t=numeric()
        varb13t=numeric()
        for(i in 1:length(all_snps_H_and_y)){
          b13t[i]=mean(all_snps_H_and_y[[i]]$p_lasso1<0.05)
          varb13t[i]=var(all_snps_H_and_y[[i]]$p_lasso1<0.05)
        }
        lasso61=mean(b13t)
        sdlasso61=sqrt(mean(varb13t))
        
        
        b14t=numeric()
        varb14t=numeric()
        for(i in 1:length(all_snps_H_and_y)){
          b14t[i]=mean(all_snps_H_and_y[[i]]$p_lasso2<0.05)
          varb14t[i]=var(all_snps_H_and_y[[i]]$p_lasso2<0.05)
        }
        lasso62=mean(b14t)
        sdlasso62=sqrt(mean(varb14t))
        
        
        
        
        power1=numeric()
        power21=numeric()
        power22=numeric()
        power_lasso21=numeric()
        power_lasso22=numeric()
        
        for(i in 1:length(all_snps_y)){
          power1[i]=mean(all_snps_y[[i]]$p_before<0.05)
          power21[i]=mean(all_snps_y[[i]]$p_after1<0.05)
          power22[i]=mean(all_snps_y[[i]]$p_after2<0.05)
          power_lasso21[i]=mean(all_snps_y[[i]]$p_lasso1<0.05)
          power_lasso22[i]=mean(all_snps_y[[i]]$p_lasso2<0.05)
        }
        
        change1=power21-power1
        change2=power22-power1
        
        if(sum(change1>0)>0){
          lasso71=power_lasso21[which.max(change1)]
          sdlasso71=sd(all_snps_y[[which.max(change1)]]$p_lasso1<0.05)
        }
        if(sum(change1>0)==0){
          lasso71=NA
          sdlasso71=NA
        }
        if(sum(change2>0)>0){
          lasso72=power_lasso22[which.max(change2)]
          sdlasso72=sd(all_snps_y[[which.max(change2)]]$p_lasso2<0.05)
        }
        if(sum(change2>0)==0){
          lasso72=NA
          sdlasso72=NA
        }
        
        
        if(sum(change1<0)>0){
          lasso81=power_lasso21[which.min(change1)]
          sdlasso81=sd(all_snps_y[[which.min(change1)]]$p_lasso1<0.05)
        }
        if(sum(change1<0)==0){
          lasso81=NA
          sdlasso81=NA
        }
        
        if(sum(change2<0)>0){
          lasso82=power_lasso22[which.min(change2)]
          sdlasso82=sd(all_snps_y[[which.min(change2)]]$p_lasso2<0.05)
        }
        if(sum(change2<0)==0){
          lasso82=NA
          sdlasso82=NA
        }
        
      }
      
      
      
      if(10>0){
        no=round(c(a11,a21,a31,a41,a51,a61,a71,a81 ),3)
        sdno=round(c(sda11,sda21,sda31,sda41,sda51,sda61,sda71,sda81 ),3)
        cml1=round(c(cml11,cml21,cml31,cml41,cml51,cml61,cml71,cml81 ),3)
        sdcml1=round(c(sdcml11,sdcml21,sdcml31,sdcml41,sdcml51,sdcml61,sdcml71,sdcml81 ),3)
        cml2=round(c(cml12,cml22,cml32,cml42,cml52,cml62,cml72,cml82 ),3)
        sdcml2=round(c(sdcml12,sdcml22,sdcml32,sdcml42,sdcml52,sdcml62,sdcml72,sdcml82 ),3)
        cml3=round(c(cml13,cml23,cml33,cml43,cml53,cml63,cml73,cml83 ),3)
        sdcml3=round(c(sdcml13,sdcml23,sdcml33,sdcml43,sdcml53,sdcml63,sdcml73,sdcml83 ),3)
        
        egger1=round(c(egger11,egger21,egger31,egger41,egger51,egger61,egger71,egger81 ),3)
        sdegger1=round(c(sdegger11,sdegger21,sdegger31,sdegger41,sdegger51,sdegger61,sdegger71,sdegger81 ),3)
        egger2=round(c(egger12,egger22,egger32,egger42,egger52,egger62,egger72,egger82 ),3)
        sdegger2=round(c(sdegger12,sdegger22,sdegger32,sdegger42,sdegger52,sdegger62,sdegger72,sdegger82 ),3)
        
        DD1=round(c(DD11,DD21,DD31,DD41,DD51,DD61,DD71,DD81 ),3)
        sdDD1=round(c(sdDD11,sdDD21,sdDD31,sdDD41,sdDD51,sdDD61,sdDD71,sdDD81 ),3)
        DD2=round(c(DD12,DD22,DD32,DD42,DD52,DD62,DD72,DD82 ),3)
        sdDD2=round(c(sdDD12,sdDD22,sdDD32,sdDD42,sdDD52,sdDD62,sdDD72,sdDD82 ),3)
        
        SH1=round(c(SH11,SH21,SH31,SH41,SH51,SH61,SH71,SH81 ),3)
        sdSH1=round(c(sdSH11,sdSH21,sdSH31,sdSH41,sdSH51,sdSH61,sdSH71,sdSH81 ),3)
        SH2=round(c(SH12,SH22,SH32,SH42,SH52,SH62,SH72,SH82 ),3)
        sdSH2=round(c(sdSH12,sdSH22,sdSH32,sdSH42,sdSH52,sdSH62,sdSH72,sdSH82 ),3)
        
        median1=round(c(median11,median21,median31,median41,median51,median61,median71,median81 ),3)
        sdmedian1=round(c(sdmedian11,sdmedian21,sdmedian31,sdmedian41,sdmedian51,sdmedian61,sdmedian71,sdmedian81 ),3)
        median2=round(c(median12,median22,median32,median42,median52,median62,median72,median82 ),3)
        sdmedian2=round(c(sdmedian12,sdmedian22,sdmedian32,sdmedian42,sdmedian52,sdmedian62,sdmedian72,sdmedian82 ),3)
        
        ivw1=round(c(ivw11,ivw21,ivw31,ivw41,ivw51,ivw61,ivw71,ivw81 ),3)
        sdivw1=round(c(sdivw11,sdivw21,sdivw31,sdivw41,sdivw51,sdivw61,sdivw71,sdivw81 ),3)
        ivw2=round(c(ivw12,ivw22,ivw32,ivw42,ivw52,ivw62,ivw72,ivw82 ),3)
        sdivw2=round(c(sdivw12,sdivw22,sdivw32,sdivw42,sdivw52,sdivw62,sdivw72,sdivw82 ),3)
        
        
        lasso1=round(c(lasso11,lasso21,lasso31,lasso41,lasso51,lasso61,lasso71,lasso81 ),3)
        sdlasso1=round(c(sdlasso11,sdlasso21,sdlasso31,sdlasso41,sdlasso51,sdlasso61,sdlasso71,sdlasso81 ),3)
        lasso2=round(c(lasso12,lasso22,lasso32,lasso42,lasso52,lasso62,lasso72,lasso82 ),3)
        sdlasso2=round(c(sdlasso12,sdlasso22,sdlasso32,sdlasso42,sdlasso52,sdlasso62,sdlasso72,sdlasso82 ),3)
        
      }
      
      
      
      
      
      
      
      result1=data.frame(no,sdno,cml1,sdcml1,cml2,sdcml2,cml3,sdcml3,egger1,sdegger1,egger2,sdegger2,DD1,sdDD1,DD2,sdDD2,SH1,sdSH1,SH2,sdSH2,median1,sdmedian1,median2,sdmedian2,ivw1,sdivw1,ivw2,sdivw2,lasso1,sdlasso1,lasso2,sdlasso2)
      slope_estimated=data.frame(b_cml=mean(all_snps[[1]]$cml_slop1),b_cml_sd=sd(all_snps[[1]]$cml_slop1),b_cml_se=mean(all_snps[[1]]$cml_slop_se1),b_egger=mean(all_snps[[1]]$egger_slop1),b_egger_sd=sd(all_snps[[1]]$egger_slop1),b_egger_se=mean(all_snps[[1]]$egger_slop_se1),b_DD=mean(all_snps[[1]]$DD_slop1),b_DD_sd=sd(all_snps[[1]]$DD_slop1),b_DD_se=mean(all_snps[[1]]$DD_slop_se1),b_SH=mean(all_snps[[1]]$SH_slop1),b_SH_sd=sd(all_snps[[1]]$SH_slop1),b_SH_se=mean(all_snps[[1]]$SH_slop_se1),b_median=mean(all_snps[[1]]$median_slop1),b_median_sd=sd(all_snps[[1]]$median_slop1),b_median_se=mean(all_snps[[1]]$median_slop_se1),b_lasso=mean(all_snps[[1]]$lasso_slop1),b_lasso_sd=sd(all_snps[[1]]$lasso_slop1),b_lasso_se=mean(all_snps[[1]]$lasso_slop_se1),b_ivw=mean(all_snps[[1]]$ivw_slop1),b_ivw_sd=sd(all_snps[[1]]$ivw_slop1),b_ivw_se=mean(all_snps[[1]]$ivw_slop_se1),b_true=mean(all_snps[[1]]$true_slop1))
      slope_estimated=round(slope_estimated,digits=3)
      
      H_only=data.frame(beta_true=rep(0,50),beta_before=rep(0,50),sd_before=rep(0,50),se_before=rep(0,50),error_before=rep(0,50),
                        beta_cml=rep(0,50),sd_cml=rep(0,50),se_cml1=rep(0,50),se_cml2=rep(0,50),se_cml3=rep(0,50),error_cml1=rep(0,50),error_cml2=rep(0,50),error_cml3=rep(0,50),
                        beta_egger=rep(0,50),sd_egger=rep(0,50),se_egger1=rep(0,50),se_egger2=rep(0,50),error_egger1=rep(0,50),error_egger2=rep(0,50),
                        beta_DD=rep(0,50),sd_DD=rep(0,50),se_DD1=rep(0,50),se_DD2=rep(0,50),error_DD1=rep(0,50),error_DD2=rep(0,50),
                        beta_SH=rep(0,50),sd_SH=rep(0,50),se_SH1=rep(0,50),se_SH2=rep(0,50),error_SH1=rep(0,50),error_SH2=rep(0,50),
                        beta_ivw=rep(0,50),sd_ivw=rep(0,50),se_ivw1=rep(0,50),se_ivw2=rep(0,50),error_ivw1=rep(0,50),error_ivw2=rep(0,50),
                        beta_lasso=rep(0,50),sd_lasso=rep(0,50),se_lasso1=rep(0,50),se_lasso2=rep(0,50),error_lasso1=rep(0,50),error_lasso2=rep(0,50))
      
      
      
      for(i1 in 1:50){
        H_only$beta_true[i1]=mean(all_snps[[i1]]$beta_true)
        H_only$beta_before[i1]=mean(all_snps[[i1]]$beta_before)
        H_only$sd_before[i1]=sd(all_snps[[i1]]$beta_before)
        H_only$se_before[i1]=mean(all_snps[[i1]]$se_before)
        H_only$error_before[i1]=mean(all_snps[[i1]]$p_before<0.05)
        H_only$beta_cml[i1]=mean(all_snps[[i1]]$beta_after)
        H_only$sd_cml[i1]=sd(all_snps[[i1]]$beta_after)
        H_only$se_cml1[i1]=mean(all_snps[[i1]]$se_after1)
        H_only$se_cml2[i1]=mean(all_snps[[i1]]$se_after2)
        H_only$se_cml3[i1]=mean(all_snps[[i1]]$se_after3)
        H_only$error_cml1[i1]=mean(all_snps[[i1]]$p_after1<0.05)
        H_only$error_cml2[i1]=mean(all_snps[[i1]]$p_after2<0.05)
        H_only$error_cml3[i1]=mean(all_snps[[i1]]$p_after3<0.05)
        H_only$beta_egger[i1]=mean(all_snps[[i1]]$beta_egger)
        H_only$sd_egger[i1]=sd(all_snps[[i1]]$beta_egger)
        H_only$se_egger1[i1]=mean(all_snps[[i1]]$se_egger1)
        H_only$se_egger2[i1]=mean(all_snps[[i1]]$se_egger2)
        H_only$error_egger1[i1]=mean(all_snps[[i1]]$p_egger1<0.05)
        H_only$error_egger2[i1]=mean(all_snps[[i1]]$p_egger2<0.05)
        H_only$beta_DD[i1]=mean(all_snps[[i1]]$beta_DD)
        H_only$sd_DD[i1]=sd(all_snps[[i1]]$beta_DD)
        H_only$se_DD1[i1]=mean(all_snps[[i1]]$se_DD1)
        H_only$se_DD2[i1]=mean(all_snps[[i1]]$se_DD2)
        H_only$error_DD1[i1]=mean(all_snps[[i1]]$p_DD1<0.05)
        H_only$error_DD2[i1]=mean(all_snps[[i1]]$p_DD2<0.05)
        H_only$beta_SH[i1]=mean(all_snps[[i1]]$beta_SH)
        H_only$sd_SH[i1]=sd(all_snps[[i1]]$beta_SH)
        H_only$se_SH1[i1]=mean(all_snps[[i1]]$se_SH1)
        H_only$se_SH2[i1]=mean(all_snps[[i1]]$se_SH2)
        H_only$error_SH1[i1]=mean(all_snps[[i1]]$p_SH1<0.05)
        H_only$error_SH2[i1]=mean(all_snps[[i1]]$p_SH2<0.05)
        H_only$beta_lasso[i1]=mean(all_snps[[i1]]$beta_lasso)
        H_only$sd_lasso[i1]=sd(all_snps[[i1]]$beta_lasso)
        H_only$se_lasso1[i1]=mean(all_snps[[i1]]$se_lasso1)
        H_only$se_lasso2[i1]=mean(all_snps[[i1]]$se_lasso2)
        H_only$error_lasso1[i1]=mean(all_snps[[i1]]$p_lasso1<0.05)
        H_only$error_lasso2[i1]=mean(all_snps[[i1]]$p_lasso2<0.05)
        H_only$beta_median[i1]=mean(all_snps[[i1]]$beta_median)
        H_only$sd_median[i1]=sd(all_snps[[i1]]$beta_median)
        H_only$se_median1[i1]=mean(all_snps[[i1]]$se_median1)
        H_only$se_median2[i1]=mean(all_snps[[i1]]$se_median2)
        H_only$error_median1[i1]=mean(all_snps[[i1]]$p_median1<0.05)
        H_only$error_median2[i1]=mean(all_snps[[i1]]$p_median2<0.05)
        H_only$beta_ivw[i1]=mean(all_snps[[i1]]$beta_ivw)
        H_only$sd_ivw[i1]=sd(all_snps[[i1]]$beta_ivw)
        H_only$se_ivw1[i1]=mean(all_snps[[i1]]$se_ivw1)
        H_only$se_ivw2[i1]=mean(all_snps[[i1]]$se_ivw2)
        H_only$error_ivw1[i1]=mean(all_snps[[i1]]$p_ivw1<0.05)
        H_only$error_ivw2[i1]=mean(all_snps[[i1]]$p_ivw2<0.05)
      }
      
      H_and_Y=data.frame(beta_true=rep(0,50),beta_before=rep(0,50),sd_before=rep(0,50),se_before=rep(0,50),power_before=rep(0,50),
                         beta_cml=rep(0,50),sd_cml=rep(0,50),se_cml1=rep(0,50),se_cml2=rep(0,50),se_cml3=rep(0,50),power_cml1=rep(0,50),power_cml2=rep(0,50),power_cml3=rep(0,50),
                         beta_egger=rep(0,50),sd_egger=rep(0,50),se_egger1=rep(0,50),se_egger2=rep(0,50),power_egger1=rep(0,50),power_egger2=rep(0,50),
                         beta_DD=rep(0,50),sd_DD=rep(0,50),se_DD1=rep(0,50),se_DD2=rep(0,50),power_DD1=rep(0,50),power_DD2=rep(0,50),
                         beta_SH=rep(0,50),sd_SH=rep(0,50),se_SH1=rep(0,50),se_SH2=rep(0,50),power_SH1=rep(0,50),power_SH2=rep(0,50),
                         beta_ivw=rep(0,50),sd_ivw=rep(0,50),se_ivw1=rep(0,50),se_ivw2=rep(0,50),power_ivw1=rep(0,50),power_ivw2=rep(0,50),
                         beta_lasso=rep(0,50),sd_lasso=rep(0,50),se_lasso1=rep(0,50),se_lasso2=rep(0,50),power_lasso1=rep(0,50),power_lasso2=rep(0,50))
      
      
      for(i1 in 101:150){
        H_and_Y$beta_true[i1-100]=mean(all_snps[[i1]]$beta_true)
        H_and_Y$beta_before[i1-100]=mean(all_snps[[i1]]$beta_before)
        H_and_Y$sd_before[i1-100]=sd(all_snps[[i1]]$beta_before)
        H_and_Y$se_before[i1-100]=mean(all_snps[[i1]]$se_before)
        H_and_Y$power_before[i1-100]=mean(all_snps[[i1]]$p_before<0.05)
        H_and_Y$beta_cml[i1-100]=mean(all_snps[[i1]]$beta_after)
        H_and_Y$sd_cml[i1-100]=sd(all_snps[[i1]]$beta_after)
        H_and_Y$se_cml1[i1-100]=mean(all_snps[[i1]]$se_after1)
        H_and_Y$se_cml2[i1-100]=mean(all_snps[[i1]]$se_after2)
        H_and_Y$se_cml3[i1-100]=mean(all_snps[[i1]]$se_after3)
        H_and_Y$power_cml1[i1-100]=mean(all_snps[[i1]]$p_after1<0.05)
        H_and_Y$power_cml2[i1-100]=mean(all_snps[[i1]]$p_after2<0.05)
        H_and_Y$power_cml3[i1-100]=mean(all_snps[[i1]]$p_after3<0.05)
        H_and_Y$beta_egger[i1-100]=mean(all_snps[[i1]]$beta_egger)
        H_and_Y$sd_egger[i1-100]=sd(all_snps[[i1]]$beta_egger)
        H_and_Y$se_egger1[i1-100]=mean(all_snps[[i1]]$se_egger1)
        H_and_Y$se_egger2[i1-100]=mean(all_snps[[i1]]$se_egger2)
        H_and_Y$power_egger1[i1-100]=mean(all_snps[[i1]]$p_egger1<0.05)
        H_and_Y$power_egger2[i1-100]=mean(all_snps[[i1]]$p_egger2<0.05)
        H_and_Y$beta_DD[i1-100]=mean(all_snps[[i1]]$beta_DD)
        H_and_Y$sd_DD[i1-100]=sd(all_snps[[i1]]$beta_DD)
        H_and_Y$se_DD1[i1-100]=mean(all_snps[[i1]]$se_DD1)
        H_and_Y$se_DD2[i1-100]=mean(all_snps[[i1]]$se_DD2)
        H_and_Y$power_DD1[i1-100]=mean(all_snps[[i1]]$p_DD1<0.05)
        H_and_Y$power_DD2[i1-100]=mean(all_snps[[i1]]$p_DD2<0.05)
        H_and_Y$beta_SH[i1-100]=mean(all_snps[[i1]]$beta_SH)
        H_and_Y$sd_SH[i1-100]=sd(all_snps[[i1]]$beta_SH)
        H_and_Y$se_SH1[i1-100]=mean(all_snps[[i1]]$se_SH1)
        H_and_Y$se_SH2[i1-100]=mean(all_snps[[i1]]$se_SH2)
        H_and_Y$power_SH1[i1-100]=mean(all_snps[[i1]]$p_SH1<0.05)
        H_and_Y$power_SH2[i1-100]=mean(all_snps[[i1]]$p_SH2<0.05)
        H_and_Y$beta_lasso[i1-100]=mean(all_snps[[i1]]$beta_lasso)
        H_and_Y$sd_lasso[i1-100]=sd(all_snps[[i1]]$beta_lasso)
        H_and_Y$se_lasso1[i1-100]=mean(all_snps[[i1]]$se_lasso1)
        H_and_Y$se_lasso2[i1-100]=mean(all_snps[[i1]]$se_lasso2)
        H_and_Y$power_lasso1[i1-100]=mean(all_snps[[i1]]$p_lasso1<0.05)
        H_and_Y$power_lasso2[i1-100]=mean(all_snps[[i1]]$p_lasso2<0.05)
        H_and_Y$beta_median[i1-100]=mean(all_snps[[i1]]$beta_median)
        H_and_Y$sd_median[i1-100]=sd(all_snps[[i1]]$beta_median)
        H_and_Y$se_median1[i1-100]=mean(all_snps[[i1]]$se_median1)
        H_and_Y$se_median2[i1-100]=mean(all_snps[[i1]]$se_median2)
        H_and_Y$power_median1[i1-100]=mean(all_snps[[i1]]$p_median1<0.05)
        H_and_Y$power_median2[i1-100]=mean(all_snps[[i1]]$p_median2<0.05)
        H_and_Y$beta_ivw[i1-100]=mean(all_snps[[i1]]$beta_ivw)
        H_and_Y$sd_ivw[i1-100]=sd(all_snps[[i1]]$beta_ivw)
        H_and_Y$se_ivw1[i1-100]=mean(all_snps[[i1]]$se_ivw1)
        H_and_Y$se_ivw2[i1-100]=mean(all_snps[[i1]]$se_ivw2)
        H_and_Y$power_ivw1[i1-100]=mean(all_snps[[i1]]$p_ivw1<0.05)
        H_and_Y$power_ivw2[i1-100]=mean(all_snps[[i1]]$p_ivw2<0.05)
      }
      
      
      
      if(10>0){
        beta_true=numeric()
        beta_before=numeric()
        beta_after=numeric()
        beta_egger=numeric()
        beta_DD=numeric()
        beta_SH=numeric()
        beta_ivw=numeric()
        beta_median=numeric()
        beta_lasso=numeric()
        
        
        sd_before=numeric()
        sd_after=numeric()
        sd_egger=numeric()
        sd_DD=numeric()
        sd_SH=numeric()
        sd_ivw=numeric()
        sd_median=numeric()
        sd_lasso=numeric()
        
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
        x0_DD=numeric()
        y0_DD=numeric()
        x1_DD=numeric()
        y1_DD=numeric()
        x0_SH=numeric()
        y0_SH=numeric()
        x1_SH=numeric()
        y1_SH=numeric()
        x0_median=numeric()
        y0_median=numeric()
        x1_median=numeric()
        y1_median=numeric()
        x0_lasso=numeric()
        y0_lasso=numeric()
        x1_lasso=numeric()
        y1_lasso=numeric()
        x0_ivw=numeric()
        y0_ivw=numeric()
        x1_ivw=numeric()
        y1_ivw=numeric()
        
      }
      
      
      
      for(i in 1:1000){
        beta_before[i]=mean(all_snps[[i]]$beta_before)
        beta_after[i]=mean(all_snps[[i]]$beta_after)
        sd_before[i]=mean(all_snps[[i]]$se_before)
        sd_after[i]=mean(all_snps[[i]]$se_after3)
        beta_true[i]=mean(all_snps[[i]]$beta_true)
        beta_egger[i]=mean(all_snps[[i]]$beta_egger)
        sd_egger[i]=mean(all_snps[[i]]$se_egger2)
        beta_DD[i]=mean(all_snps[[i]]$beta_DD)
        sd_DD[i]=mean(all_snps[[i]]$se_DD1)
        beta_SH[i]=mean(all_snps[[i]]$beta_SH)
        sd_SH[i]=mean(all_snps[[i]]$se_SH1)
        beta_median[i]=mean(all_snps[[i]]$beta_median)
        sd_median[i]=mean(all_snps[[i]]$se_median1)
        beta_lasso[i]=mean(all_snps[[i]]$beta_lasso)
        sd_lasso[i]=mean(all_snps[[i]]$se_lasso1)
        beta_ivw[i]=mean(all_snps[[i]]$beta_ivw)
        sd_ivw[i]=mean(all_snps[[i]]$se_ivw1)
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
        x0_DD[i]=beta_true[i]
        y0_DD[i]=beta_DD[i]-sd_DD[i]
        x1_DD[i]=beta_true[i]
        y1_DD[i]=beta_DD[i]+sd_DD[i]
        x0_SH[i]=beta_true[i]
        y0_SH[i]=beta_SH[i]-sd_SH[i]
        x1_SH[i]=beta_true[i]
        y1_SH[i]=beta_SH[i]+sd_SH[i]
        x0_ivw[i]=beta_true[i]
        y0_ivw[i]=beta_ivw[i]-sd_ivw[i]
        x1_ivw[i]=beta_true[i]
        y1_ivw[i]=beta_ivw[i]+sd_ivw[i]
        x0_median[i]=beta_true[i]
        y0_median[i]=beta_median[i]-sd_median[i]
        x1_median[i]=beta_true[i]
        y1_median[i]=beta_median[i]+sd_median[i]
        x0_lasso[i]=beta_true[i]
        y0_lasso[i]=beta_lasso[i]-sd_lasso[i]
        x1_lasso[i]=beta_true[i]
        y1_lasso[i]=beta_lasso[i]+sd_lasso[i]
      }
      if(10>0){
        jpeg(filename = plotname[1],width=500,height = 500,res=100)
        plot(beta_true,beta_before,xlim=c(min(beta_true),max(beta_true)),ylim=c(min(y0_before,y0_after,y0_egger,y0_DD,y0_SH,y0_ivw,y0_lasso,y0_median),max(y1_before,y1_after,y1_egger,y1_DD,y1_SH,y1_ivw,y1_lasso,y1_median)),main="SNP effects before bias correction",type="n",xlab = "true effect",ylab = "estimated effect")
        
        legend("topleft",pch=c(15,17),col=c("black","darkorchid1"),legend=c("H only","H and Y"),cex=c(0.7,0.7))
        abline(a=0,b=1)
        for(i in 1:50){
          segments(x0=x0_before[i],y0=y0_before[i],x1=x1_before[i],y1=y1_before[i],lwd=0.5,col="gray")
        }
        for(i in 101:150){
          segments(x0=x0_before[i],y0=y0_before[i],x1=x1_before[i],y1=y1_before[i],lwd=0.5,col="gray")
        }
        points(beta_true[1:50],beta_before[1:50],col=1,pch=15,cex=0.7)
        points(beta_true[101:150],beta_before[101:150],col="darkorchid1",pch=17,cex=0.7)
        dev.off()
        
      }
      if(10>0){
        jpeg(filename = plotname[2],width=500,height = 500,res=100)
        
        plot(beta_true,beta_after,xlim=c(min(beta_true),max(beta_true)),ylim=c(min(y0_before,y0_after,y0_egger,y0_DD,y0_SH,y0_ivw,y0_lasso,y0_median),max(y1_before,y1_after,y1_egger,y1_DD,y1_SH,y1_ivw,y1_lasso,y1_median)),main="SNP effects after MVMR-cML bias correction",type="n",xlab = "true effect",ylab = "estimated effect")
        
        legend("topleft",pch=c(15,17),col=c("black","darkorchid1"),legend=c("H only","H and Y"),cex=c(0.7,0.7))
        abline(a=0,b=1)
        for(i in 1:50){
          segments(x0=x0_after[i],y0=y0_after[i],x1=x1_after[i],y1=y1_after[i],lwd=0.5,col="gray")
        }
        for(i in 101:150){
          segments(x0=x0_after[i],y0=y0_after[i],x1=x1_after[i],y1=y1_after[i],lwd=0.5,col="gray")
        }
        
        
        points(beta_true[1:50],beta_after[1:50],col=1,pch=15,cex=0.7)
        points(beta_true[101:150],beta_after[101:150],col="darkorchid1",pch=17,cex=0.7)
        dev.off()
      }
      if(10>0){
        jpeg(filename = plotname[3],width=500,height = 500,res=100)
        
        plot(beta_true,beta_egger,xlim=c(min(beta_true),max(beta_true)),ylim=c(min(y0_before,y0_after,y0_egger,y0_DD,y0_SH,y0_ivw,y0_lasso,y0_median),max(y1_before,y1_after,y1_egger,y1_DD,y1_SH,y1_ivw,y1_lasso,y1_median)),main="SNP effects after MV-Egger bias correction",type="n",xlab = "true effect",ylab = "estimated effect")
        
        legend("topleft",pch=c(15,17),col=c("black","darkorchid1"),legend=c("H only","H and Y"),cex=c(0.7,0.7))
        abline(a=0,b=1)
        for(i in 1:50){
          segments(x0=x0_egger[i],y0=y0_egger[i],x1=x1_egger[i],y1=y1_egger[i],lwd=0.5,col="gray")
        }
        for(i in 101:150){
          segments(x0=x0_egger[i],y0=y0_egger[i],x1=x1_egger[i],y1=y1_egger[i],lwd=0.5,col="gray")
        }
        
        
        points(beta_true[1:50],beta_egger[1:50],col=1,pch=15,cex=0.7)
        points(beta_true[101:150],beta_egger[101:150],col="darkorchid1",pch=17,cex=0.7)
        dev.off()
        
      }
      if(10>0){
        jpeg(filename = plotname[4],width=500,height = 500,res=100)
        
        plot(beta_true,beta_DD,xlim=c(min(beta_true),max(beta_true)),ylim=c(min(y0_before,y0_after,y0_egger,y0_DD,y0_SH,y0_ivw,y0_lasso,y0_median),max(y1_before,y1_after,y1_egger,y1_DD,y1_SH,y1_ivw,y1_lasso,y1_median)),main="SNP effects after DHO bias correction",type="n",xlab = "true effect",ylab = "estimated effect")
        
        legend("topleft",pch=c(15,17),col=c("black","darkorchid1"),legend=c("H only","H and Y"),cex=c(0.7,0.7))
        abline(a=0,b=1)
        for(i in 1:50){
          segments(x0=x0_DD[i],y0=y0_DD[i],x1=x1_DD[i],y1=y1_DD[i],lwd=0.5,col="gray")
        }
        for(i in 101:150){
          segments(x0=x0_DD[i],y0=y0_DD[i],x1=x1_DD[i],y1=y1_DD[i],lwd=0.5,col="gray")
        }
        
        
        points(beta_true[1:50],beta_DD[1:50],col=1,pch=15,cex=0.7)
        points(beta_true[101:150],beta_DD[101:150],col="darkorchid1",pch=17,cex=0.7)
        dev.off()
      }
      if(10>0){
        jpeg(filename = plotname[5],width=500,height = 500,res=100)
        
        plot(beta_true,beta_SH,xlim=c(min(beta_true),max(beta_true)),ylim=c(min(y0_before,y0_after,y0_egger,y0_DD,y0_SH,y0_ivw,y0_lasso,y0_median),max(y1_before,y1_after,y1_egger,y1_DD,y1_SH,y1_ivw,y1_lasso,y1_median)),main="SNP effects after SH bias correction",type="n",xlab = "true effect",ylab = "estimated effect")
        
        legend("topleft",pch=c(15,17),col=c("black","darkorchid1"),legend=c("H only","H and Y"),cex=c(0.7,0.7))
        abline(a=0,b=1)
        for(i in 1:50){
          segments(x0=x0_SH[i],y0=y0_SH[i],x1=x1_SH[i],y1=y1_SH[i],lwd=0.5,col="gray")
        }
        for(i in 101:150){
          segments(x0=x0_SH[i],y0=y0_SH[i],x1=x1_SH[i],y1=y1_SH[i],lwd=0.5,col="gray")
        }
        
        
        points(beta_true[1:50],beta_SH[1:50],col=1,pch=15,cex=0.7)
        points(beta_true[101:150],beta_SH[101:150],col="darkorchid1",pch=17,cex=0.7)
        dev.off()
      }
      if(10>0){
        jpeg(filename = plotname[6],width=500,height = 500,res=100)
        
        plot(beta_true,beta_ivw,xlim=c(min(beta_true),max(beta_true)),ylim=c(min(y0_before,y0_after,y0_egger,y0_DD,y0_ivw,y0_ivw,y0_lasso,y0_median),max(y1_before,y1_after,y1_egger,y1_DD,y1_ivw,y1_ivw,y1_lasso,y1_median)),main="SNP effects after MVMR-IVW bias correction",type="n",xlab = "true effect",ylab = "estimated effect")
        
        legend("topleft",pch=c(15,17),col=c("black","darkorchid1"),legend=c("H only","H and Y"),cex=c(0.7,0.7))
        abline(a=0,b=1)
        for(i in 1:50){
          segments(x0=x0_ivw[i],y0=y0_ivw[i],x1=x1_ivw[i],y1=y1_ivw[i],lwd=0.5,col="gray")
        }
        for(i in 101:150){
          segments(x0=x0_ivw[i],y0=y0_ivw[i],x1=x1_ivw[i],y1=y1_ivw[i],lwd=0.5,col="gray")
        }
        
        
        points(beta_true[1:50],beta_ivw[1:50],col=1,pch=15,cex=0.7)
        points(beta_true[101:150],beta_ivw[101:150],col="darkorchid1",pch=17,cex=0.7)
        dev.off()
      }
      if(10>0){
        jpeg(filename = plotname[7],width=500,height = 500,res=100)
        
        plot(beta_true,beta_lasso,xlim=c(min(beta_true),max(beta_true)),ylim=c(min(y0_before,y0_after,y0_egger,y0_DD,y0_lasso,y0_lasso,y0_lasso,y0_lasso),max(y1_before,y1_after,y1_egger,y1_DD,y1_lasso,y1_lasso,y1_lasso,y1_lasso)),main="SNP effects after MVMR-LASSO bias correction",type="n",xlab = "true effect",ylab = "estimated effect")
        
        legend("topleft",pch=c(15,17),col=c("black","darkorchid1"),legend=c("H only","H and Y"),cex=c(0.7,0.7))
        abline(a=0,b=1)
        for(i in 1:50){
          segments(x0=x0_lasso[i],y0=y0_lasso[i],x1=x1_lasso[i],y1=y1_lasso[i],lwd=0.5,col="gray")
        }
        for(i in 101:150){
          segments(x0=x0_lasso[i],y0=y0_lasso[i],x1=x1_lasso[i],y1=y1_lasso[i],lwd=0.5,col="gray")
        }
        
        
        points(beta_true[1:50],beta_lasso[1:50],col=1,pch=15,cex=0.7)
        points(beta_true[101:150],beta_lasso[101:150],col="darkorchid1",pch=17,cex=0.7)
        dev.off()
        
      }
      if(10>0){
        jpeg(filename = plotname[8],width=500,height = 500,res=100)
        
        plot(beta_true,beta_median,xlim=c(min(beta_true),max(beta_true)),ylim=c(min(y0_before,y0_after,y0_egger,y0_DD,y0_median,y0_median,y0_median,y0_median),max(y1_before,y1_after,y1_egger,y1_DD,y1_median,y1_median,y1_median,y1_median)),main="SNP effects after MVMR-Median bias correction",type="n",xlab = "true effect",ylab = "estimated effect")
        
        legend("topleft",pch=c(15,17),col=c("black","darkorchid1"),legend=c("H only","H and Y"),cex=c(0.7,0.7))
        abline(a=0,b=1)
        for(i in 1:50){
          segments(x0=x0_median[i],y0=y0_median[i],x1=x1_median[i],y1=y1_median[i],lwd=0.5,col="gray")
        }
        for(i in 101:150){
          segments(x0=x0_median[i],y0=y0_median[i],x1=x1_median[i],y1=y1_median[i],lwd=0.5,col="gray")
        }
        
        
        points(beta_true[1:50],beta_median[1:50],col=1,pch=15,cex=0.7)
        points(beta_true[101:150],beta_median[101:150],col="darkorchid1",pch=17,cex=0.7)
        dev.off()
        
      }
      
      
      
      point_est=list(H_only,H_and_Y)
      write.table(result1,file = filename1,sep="\t",col.names = TRUE,quote = FALSE,row.names = FALSE)
      write.table(slope_estimated,file = filename2,sep="\t",col.names = TRUE,quote = FALSE,row.names = FALSE)
      save(point_est,file=filename3)
      
    }
  }
  print(d)
}


for(d in c(1)){
  for(rho in c(0,0.5,-0.5)){
    if(10>0){
      N=20
      setwd(paste("D:/art sim/New/real_snps_sim/30invalid/","d=",d,sep=""))
      plotname=numeric()
      name1_1=numeric()
      plotname[1]=paste("D:/art sim/New/real_snps_sim/result/30invalid/","d=",d,"/","d=",d,"betaXYnot0rho",rho,"_before.jpeg",sep="")
      plotname[2]=paste("D:/art sim/New/real_snps_sim/result/30invalid/","d=",d,"/","d=",d,"betaXYnot0rho",rho,"_after.jpeg",sep="")
      plotname[3]=paste("D:/art sim/New/real_snps_sim/result/30invalid/","d=",d,"/","d=",d,"betaXYnot0rho",rho,"_egger.jpeg",sep="")
      plotname[4]=paste("D:/art sim/New/real_snps_sim/result/30invalid/","d=",d,"/","d=",d,"betaXYnot0rho",rho,"_dd.jpeg",sep="")
      plotname[5]=paste("D:/art sim/New/real_snps_sim/result/30invalid/","d=",d,"/","d=",d,"betaXYnot0rho",rho,"_sh.jpeg",sep="")
      plotname[6]=paste("D:/art sim/New/real_snps_sim/result/30invalid/","d=",d,"/","d=",d,"betaXYnot0rho",rho,"_ivw.jpeg",sep="")
      plotname[7]=paste("D:/art sim/New/real_snps_sim/result/30invalid/","d=",d,"/","d=",d,"betaXYnot0rho",rho,"_lasso.jpeg",sep="")
      plotname[8]=paste("D:/art sim/New/real_snps_sim/result/30invalid/","d=",d,"/","d=",d,"betaXYnot0rho",rho,"_median.jpeg",sep="")
      filename1=paste("D:/art sim/New/real_snps_sim/result/30invalid/","d=",d,"/","d=",d,"betaXYnot0rho",rho,".txt",sep="")
      filename2=paste("D:/art sim/New/real_snps_sim/result/30invalid/","d=",d,"/","d=",d,"betaXYnot0rho",rho,"Slope.txt",sep="")
      filename3=paste("D:/art sim/New/real_snps_sim/result/30invalid/","d=",d,"/","d=",d,"betaXYnot0rho",rho,"point_effect_estimates.RData",sep="")
      
      for(i in 1:N){
        name1_1[i]=paste("sim",rho,"_",i,"betaXYnot",0,".txt",sep="")
      }
      result=list()
      for(i in 1:N){
        a1=read.table(file=name1_1[i],header = TRUE,sep="\t")
        result[[i]]=a1
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
        all_snps[[i]]=a[1:1000,]
      }
      for(i in 1:1000){
        for(j in  1:1000){
          all_snps[[i]][j,]=a[i+(j-1)*1000,]
        }
        print(i)
      }
      if(10>0){
        name6se=numeric()
        name7se=numeric()
        name8se=numeric()
        for(i in 1:d){
          name6se[i]=paste("ivw_slop_se",i,sep="")
          name7se[i]=paste("lasso_slop_se",i,sep="")
          name8se[i]=paste("median_slop_se",i,sep="")
        }
        for(i in 1:length(all_snps)){
          colnames(all_snps[[i]])=c(colnames(all_snps[[i]])[1:(ncol(all_snps[[i]])-3*d)],name6se,name7se,name8se)
        }
      }
      null_snps=list()
      for(i in 1:850){
        null_snps[[i]]=all_snps[[i+150]]
      }
      
      all_snps_no_y=list()
      for(i in c(1:900)){
        if(i<=50){
          all_snps_no_y[[i]]=all_snps[[i]]
        }
        if(i>50){
          all_snps_no_y[[i]]=all_snps[[i+100]]
        }
      }
      
      
      all_snps_H_no_y=list()
      for(i in 1:50){
        all_snps_H_no_y[[i]]=all_snps[[i]]
      }
      
      all_snps_y_only=list()
      for(i in 1:50){
        all_snps_y_only[[i]]=all_snps[[50+i]]
      }
      all_snps_H_and_y=list()
      for(i in 1:50){
        all_snps_H_and_y[[i]]=all_snps[[100+i]]
      }
      
      all_snps_y=list()
      for(i in 1:100){
        all_snps_y[[i]]=all_snps[[50+i]]
      }
      
      if(10>0){
        a1=numeric()
        vara1=numeric()
        for(i in 1: length(null_snps)){
          a1[i]=mean(null_snps[[i]]$p_before<0.05)
          vara1[i]=var(null_snps[[i]]$p_before<0.05)
        }
        a11=mean(a1)
        sda11=sqrt(mean(vara1))
        a2=numeric()
        vara2=numeric()
        for(i in 1: length(all_snps_H_no_y)){
          a2[i]=mean(all_snps_H_no_y[[i]]$p_before<0.05)
          vara2[i]=var(all_snps_H_no_y[[i]]$p_before<0.05)
        }
        a21=mean(a2)
        sda21=sqrt(mean(vara2))
        
        
        a3=numeric()
        vara3=numeric()
        for(i in 1:length(all_snps_no_y)){
          a3[i]=mean(all_snps_no_y[[i]]$p_before<0.05)
          vara3[i]=var(all_snps_no_y[[i]]$p_before<0.05)
        }
        a31=max(a3)
        sda31=sqrt(vara3[which.max(a3)])
        
        
        a4=numeric()
        for(i in 1:1000){
          a=numeric()
          for(j in 1:50){
            a[j]=all_snps_H_no_y[[j]]$p_before[i]
          }
          a4[i]=(sum(a<0.05/50)>0)
        }
        a41=mean(a4)
        sda41=sd(a4) 
        
        a5=numeric()
        vara5=numeric()
        for(i in 1:length(all_snps_y_only)){
          a5[i]=mean(all_snps_y_only[[i]]$p_before<0.05)
          vara5[i]=var(all_snps_y_only[[i]]$p_before<0.05)
        }
        a51=mean(a5)
        sda51=sqrt(mean(vara5))
        a6=numeric()
        vara6=numeric()
        for(i in 1:length(all_snps_H_and_y)){
          a6[i]=mean(all_snps_H_and_y[[i]]$p_before<0.05)
          vara6[i]=var(all_snps_H_and_y[[i]]$p_before<0.05)
        }
        a61=mean(a6)
        sda61=sqrt(mean(vara6))
        power1=numeric()
        power2=numeric()
        
        for(i in 1:length(all_snps_y)){
          power1[i]=mean(all_snps_y[[i]]$p_before<0.05)
          power2[i]=mean(all_snps_y[[i]]$p_after3<0.05)
        }
        change=power2-power1
        if(sum(change>0)>0){
          a71=power1[which.max(change)]
          sda71=sd(all_snps_y[[which.max(change)]]$p_before<0.05)
        }
        if(sum(change>0)==0){
          a71=NA
          sda71=NA
        }
        
        if(sum(change<0)>0){
          a81=power1[which.min(change)]
          sda81=sd(all_snps_y[[which.min(change)]]$p_before<0.05)
        }
        if(sum(change<0)==0){
          a81=NA
          sda81=NA
        }
        
      }
      
      if(10>0){
        b1=numeric()
        varb1=numeric()
        for(i in 1: length(null_snps)){
          b1[i]=mean(null_snps[[i]]$p_after1<0.05)
          varb1[i]=var(null_snps[[i]]$p_after1<0.05)
        }
        cml11=mean(b1)
        sdcml11=sqrt(mean(varb1))
        
        
        b2=numeric()
        varb2=numeric()
        for(i in 1: length(null_snps)){
          b2[i]=mean(null_snps[[i]]$p_after2<0.05)
          varb2[i]=var(null_snps[[i]]$p_after2<0.05)
        }
        cml12=mean(b2)
        sdcml12=sqrt(mean(varb2))
        
        b3=numeric()
        varb3=numeric()
        for(i in 1: length(null_snps)){
          b3[i]=mean(null_snps[[i]]$p_after3<0.05)
          varb3[i]=var(null_snps[[i]]$p_after3<0.05)
        }
        cml13=mean(b3)
        sdcml13=sqrt(mean(varb3))
        
        
        
        
        
        
        b4=numeric()
        varb4=numeric()
        for(i in 1: length(all_snps_H_no_y)){
          b4[i]=mean(all_snps_H_no_y[[i]]$p_after1<0.05)
          varb4[i]=var(all_snps_H_no_y[[i]]$p_after1<0.05)
        }
        cml21=mean(b4)
        sdcml21=sqrt(mean(varb4))
        
        
        b5=numeric()
        varb5=numeric()
        for(i in 1: length(all_snps_H_no_y)){
          b5[i]=mean(all_snps_H_no_y[[i]]$p_after2<0.05)
          varb5[i]=var(all_snps_H_no_y[[i]]$p_after2<0.05)
        }
        cml22=mean(b5)
        sdcml22=sqrt(mean(varb5))
        
        
        b6=numeric()
        varb6=numeric()
        for(i in 1: length(all_snps_H_no_y)){
          b6[i]=mean(all_snps_H_no_y[[i]]$p_after3<0.05)
          varb6[i]=var(all_snps_H_no_y[[i]]$p_after3<0.05)
        }
        cml23=mean(b6)
        sdcml23=sqrt(mean(varb6))
        
        
        
        a=numeric()
        b=numeric()
        for(i in 1:length(all_snps_no_y)){
          a[i]=mean(all_snps_no_y[[i]]$p_before<0.05)
          b[i]=mean(all_snps_no_y[[i]]$p_after1<0.05)
        }
        cml31=b[which.max(a)]
        sdcml31=sd(all_snps_no_y[[which.max(a)]]$p_after1<0.05)
        
        a=numeric()
        b=numeric()
        for(i in 1:length(all_snps_no_y)){
          a[i]=mean(all_snps_no_y[[i]]$p_before<0.05)
          b[i]=mean(all_snps_no_y[[i]]$p_after2<0.05)
        }
        cml32=b[which.max(a)]
        sdcml32=sd(all_snps_no_y[[which.max(a)]]$p_after2<0.05)
        
        
        
        a=numeric()
        b=numeric()
        for(i in 1:length(all_snps_no_y)){
          a[i]=mean(all_snps_no_y[[i]]$p_before<0.05)
          b[i]=mean(all_snps_no_y[[i]]$p_after3<0.05)
        }
        cml33=b[which.max(a)]
        sdcml33=sd(all_snps_no_y[[which.max(a)]]$p_after3<0.05)
        
        
        b7=numeric()
        for(i in 1:1000){
          a=numeric()
          for(j in 1:50){
            a[j]=all_snps_H_no_y[[j]]$p_after1[i]
          }
          b7[i]=(sum(a<0.05/50)>0)
        }
        cml41=mean(b7)
        sdcml41=sd(b7) 
        
        b8=numeric()
        for(i in 1:1000){
          a=numeric()
          for(j in 1:50){
            a[j]=all_snps_H_no_y[[j]]$p_after2[i]
          }
          b8[i]=(sum(a<0.05/50)>0)
        }
        cml42=mean(b8)
        sdcml42=sd(b8) 
        
        
        b9=numeric()
        for(i in 1:1000){
          a=numeric()
          for(j in 1:50){
            a[j]=all_snps_H_no_y[[j]]$p_after3[i]
          }
          b9[i]=(sum(a<0.05/50)>0)
        }
        cml43=mean(b9)
        sdcml43=sd(b9)
        
        
        
        
        b10=numeric()
        varb10=numeric()
        for(i in 1:length(all_snps_y_only)){
          b10[i]=mean(all_snps_y_only[[i]]$p_after1<0.05)
          varb10[i]= var(all_snps_y_only[[i]]$p_after1<0.05)
        }
        cml51=mean(b10)
        sdcml51=sqrt(mean(varb10))
        
        b11t=numeric()
        varb11t=numeric()
        for(i in 1:length(all_snps_y_only)){
          b11t[i]=mean(all_snps_y_only[[i]]$p_after2<0.05)
          varb11t[i]=var(all_snps_y_only[[i]]$p_after2<0.05)
        }
        cml52=mean(b11t)
        sdcml52=sqrt(mean(varb11t))
        
        b12t=numeric()
        varb12t=numeric()
        for(i in 1:length(all_snps_y_only)){
          b12t[i]=mean(all_snps_y_only[[i]]$p_after3<0.05)
          varb12t[i]=var(all_snps_y_only[[i]]$p_after3<0.05)
        }
        cml53=mean(b12t)
        sdcml53=sqrt(mean(varb12t))
        
        
        
        b13t=numeric()
        varb13t=numeric()
        for(i in 1:length(all_snps_H_and_y)){
          b13t[i]=mean(all_snps_H_and_y[[i]]$p_after1<0.05)
          varb13t[i]=var(all_snps_H_and_y[[i]]$p_after1<0.05)
        }
        cml61=mean(b13t)
        sdcml61=sqrt(mean(varb13t))
        
        
        b14t=numeric()
        varb14t=numeric()
        for(i in 1:length(all_snps_H_and_y)){
          b14t[i]=mean(all_snps_H_and_y[[i]]$p_after2<0.05)
          varb14t[i]=var(all_snps_H_and_y[[i]]$p_after2<0.05)
        }
        cml62=mean(b14t)
        sdcml62=sqrt(mean(varb14t))
        
        b15t=numeric()
        varb15t=numeric()
        for(i in 1:length(all_snps_H_and_y)){
          b15t[i]=mean(all_snps_H_and_y[[i]]$p_after3<0.05)
          varb15t[i]=var(all_snps_H_and_y[[i]]$p_after3<0.05)
        }
        cml63=mean(b15t)
        sdcml63=sqrt(mean(varb15t))
        
        
        
        power1=numeric()
        power21=numeric()
        power22=numeric()
        power23=numeric()
        
        
        for(i in 1:length(all_snps_y)){
          power1[i]=mean(all_snps_y[[i]]$p_before<0.05)
          power21[i]=mean(all_snps_y[[i]]$p_after1<0.05)
          power22[i]=mean(all_snps_y[[i]]$p_after2<0.05)
          power23[i]=mean(all_snps_y[[i]]$p_after3<0.05)
        }
        change1=power21-power1
        change2=power22-power1
        change3=power23-power1
        if(sum(change1>0)>0){
          cml71=power21[which.max(change1)]
          sdcml71=sd(all_snps_y[[which.max(change1)]]$p_after1<0.05)
        }
        if(sum(change1>0)==0){
          cml71=NA
          sdcml71=NA
        }
        if(sum(change2>0)>0){
          cml72=power22[which.max(change2)]
          sdcml72=sd(all_snps_y[[which.max(change2)]]$p_after2<0.05)
        }
        if(sum(change2>0)==0){
          cml72=NA
          sdcml72=NA
        }
        if(sum(change3>0)>0){
          cml73=power23[which.max(change3)]
          sdcml73=sd(all_snps_y[[which.max(change3)]]$p_after3<0.05)
        }
        if(sum(change3>0)==0){
          cml73=NA
          sdcml73=NA
        }
        
        
        
        if(sum(change1<0)>0){
          cml81=power21[which.min(change1)]
          sdcml81=sd(all_snps_y[[which.min(change1)]]$p_after1<0.05)
        }
        if(sum(change1<0)==0){
          cml81=NA
          sdcml81=NA
        }
        
        if(sum(change2<0)>0){
          cml82=power22[which.min(change2)]
          sdcml82=sd(all_snps_y[[which.min(change2)]]$p_after2<0.05)
        }
        if(sum(change2<0)==0){
          cml82=NA
          sdcml82=NA
        }
        
        
        if(sum(change3<0)>0){
          cml83=power23[which.min(change3)]
          sdcml83=sd(all_snps_y[[which.min(change3)]]$p_after3<0.05)
        }
        if(sum(change3<0)==0){
          cml83=NA
          sdcml83=NA
        }
      }
      
      if(10>0){
        b1=numeric()
        varb1=numeric()
        for(i in 1: length(null_snps)){
          b1[i]=mean(null_snps[[i]]$p_egger1<0.05)
          varb1[i]=var(null_snps[[i]]$p_egger1<0.05)
        }
        egger11=mean(b1)
        sdegger11=sqrt(mean(varb1))
        
        
        b2=numeric()
        varb2=numeric()
        for(i in 1: length(null_snps)){
          b2[i]=mean(null_snps[[i]]$p_egger2<0.05)
          varb2[i]=var(null_snps[[i]]$p_egger2<0.05)
        }
        egger12=mean(b2)
        sdegger12=sqrt(mean(varb2))
        
        
        
        
        b4=numeric()
        varb4=numeric()
        for(i in 1: length(all_snps_H_no_y)){
          b4[i]=mean(all_snps_H_no_y[[i]]$p_egger1<0.05)
          varb4[i]=var(all_snps_H_no_y[[i]]$p_egger1<0.05)
        }
        egger21=mean(b4)
        sdegger21=sqrt(mean(varb4))
        
        
        b5=numeric()
        varb5=numeric()
        for(i in 1: length(all_snps_H_no_y)){
          b5[i]=mean(all_snps_H_no_y[[i]]$p_egger2<0.05)
          varb5[i]=var(all_snps_H_no_y[[i]]$p_egger2<0.05)
        }
        egger22=mean(b5)
        sdegger22=sqrt(mean(varb5))
        
        
        
        
        
        a=numeric()
        b=numeric()
        for(i in 1:length(all_snps_no_y)){
          a[i]=mean(all_snps_no_y[[i]]$p_before<0.05)
          b[i]=mean(all_snps_no_y[[i]]$p_egger1<0.05)
        }
        egger31=b[which.max(a)]
        sdegger31=sd(all_snps_no_y[[which.max(a)]]$p_egger1<0.05)
        
        a=numeric()
        b=numeric()
        for(i in 1:length(all_snps_no_y)){
          a[i]=mean(all_snps_no_y[[i]]$p_before<0.05)
          b[i]=mean(all_snps_no_y[[i]]$p_egger2<0.05)
        }
        egger32=b[which.max(a)]
        sdegger32=sd(all_snps_no_y[[which.max(a)]]$p_egger2<0.05)
        
        
        
        
        
        b7=numeric()
        for(i in 1:1000){
          a=numeric()
          for(j in 1:50){
            a[j]=all_snps_H_no_y[[j]]$p_egger1[i]
          }
          b7[i]=(sum(a<0.05/50)>0)
        }
        egger41=mean(b7)
        sdegger41=sd(b7) 
        
        b8=numeric()
        for(i in 1:1000){
          a=numeric()
          for(j in 1:50){
            a[j]=all_snps_H_no_y[[j]]$p_egger2[i]
          }
          b8[i]=(sum(a<0.05/50)>0)
        }
        egger42=mean(b8)
        sdegger42=sd(b8) 
        
        
        
        
        
        
        b10=numeric()
        varb10=numeric()
        for(i in 1:length(all_snps_y_only)){
          b10[i]=mean(all_snps_y_only[[i]]$p_egger1<0.05)
          varb10[i]=var(all_snps_y_only[[i]]$p_egger1<0.05)
        }
        egger51=mean(b10)
        sdegger51=sqrt(mean(varb10))
        
        b11t=numeric()
        varb11t=numeric()
        for(i in 1:length(all_snps_y_only)){
          b11t[i]=mean(all_snps_y_only[[i]]$p_egger2<0.05)
          varb11t[i]=var(all_snps_y_only[[i]]$p_egger2<0.05)
        }
        egger52=mean(b11t)
        sdegger52=sqrt(mean(varb11t))
        
        
        
        
        
        b13t=numeric()
        varb13t=numeric()
        for(i in 1:length(all_snps_H_and_y)){
          b13t[i]=mean(all_snps_H_and_y[[i]]$p_egger1<0.05)
          varb13t[i]=var(all_snps_H_and_y[[i]]$p_egger1<0.05)
        }
        egger61=mean(b13t)
        sdegger61=sqrt(mean(varb13t))
        
        
        b14t=numeric()
        varb14t=numeric()
        for(i in 1:length(all_snps_H_and_y)){
          b14t[i]=mean(all_snps_H_and_y[[i]]$p_egger2<0.05)
          varb14t[i]=var(all_snps_H_and_y[[i]]$p_egger2<0.05)
        }
        egger62=mean(b14t)
        sdegger62=sqrt(mean(varb14t))
        
        
        
        
        power1=numeric()
        power21=numeric()
        power22=numeric()
        power_egger21=numeric()
        power_egger22=numeric()
        
        for(i in 1:length(all_snps_y)){
          power1[i]=mean(all_snps_y[[i]]$p_before<0.05)
          power21[i]=mean(all_snps_y[[i]]$p_after1<0.05)
          power22[i]=mean(all_snps_y[[i]]$p_after2<0.05)
          power_egger21[i]=mean(all_snps_y[[i]]$p_egger1<0.05)
          power_egger22[i]=mean(all_snps_y[[i]]$p_egger2<0.05)
        }
        
        change1=power21-power1
        change2=power22-power1
        
        if(sum(change1>0)>0){
          egger71=power_egger21[which.max(change1)]
          sdegger71=sd(all_snps_y[[which.max(change1)]]$p_egger1<0.05)
        }
        if(sum(change1>0)==0){
          egger71=NA
          sdegger71=NA
        }
        if(sum(change2>0)>0){
          egger72=power_egger22[which.max(change2)]
          sdegger72=sd(all_snps_y[[which.max(change2)]]$p_egger2<0.05)
        }
        if(sum(change2>0)==0){
          egger72=NA
          sdegger72=NA
        }
        
        
        if(sum(change1<0)>0){
          egger81=power_egger21[which.min(change1)]
          sdegger81=sd(all_snps_y[[which.min(change1)]]$p_egger1<0.05)
        }
        if(sum(change1<0)==0){
          egger81=NA
          sdegger81=NA
        }
        
        if(sum(change2<0)>0){
          egger82=power_egger22[which.min(change2)]
          sdegger82=sd(all_snps_y[[which.min(change2)]]$p_egger2<0.05)
        }
        if(sum(change2<0)==0){
          egger82=NA
          sdegger82=NA
        }
        
      }
      
      if(10>0){
        b1=numeric()
        varb1=numeric()
        for(i in 1: length(null_snps)){
          b1[i]=mean(null_snps[[i]]$p_DD1<0.05)
          varb1[i]=var(null_snps[[i]]$p_DD1<0.05)
        }
        DD11=mean(b1)
        sdDD11=sqrt(mean(varb1))
        
        
        b2=numeric()
        varb2=numeric()
        for(i in 1: length(null_snps)){
          b2[i]=mean(null_snps[[i]]$p_DD2<0.05)
          varb2[i]=var(null_snps[[i]]$p_DD2<0.05)
        }
        DD12=mean(b2)
        sdDD12=sqrt(mean(varb2))
        
        
        
        
        b4=numeric()
        varb4=numeric()
        for(i in 1: length(all_snps_H_no_y)){
          b4[i]=mean(all_snps_H_no_y[[i]]$p_DD1<0.05)
          varb4[i]=var(all_snps_H_no_y[[i]]$p_DD1<0.05)
        }
        DD21=mean(b4)
        sdDD21=sqrt(mean(varb4))
        
        
        b5=numeric()
        varb5=numeric()
        for(i in 1: length(all_snps_H_no_y)){
          b5[i]=mean(all_snps_H_no_y[[i]]$p_DD2<0.05)
          varb5[i]=var(all_snps_H_no_y[[i]]$p_DD2<0.05)
        }
        DD22=mean(b5)
        sdDD22=sqrt(mean(varb5))
        
        
        
        
        
        a=numeric()
        b=numeric()
        for(i in 1:length(all_snps_no_y)){
          a[i]=mean(all_snps_no_y[[i]]$p_before<0.05)
          b[i]=mean(all_snps_no_y[[i]]$p_DD1<0.05)
        }
        DD31=b[which.max(a)]
        sdDD31=sd(all_snps_no_y[[which.max(a)]]$p_DD1<0.05)
        
        a=numeric()
        b=numeric()
        for(i in 1:length(all_snps_no_y)){
          a[i]=mean(all_snps_no_y[[i]]$p_before<0.05)
          b[i]=mean(all_snps_no_y[[i]]$p_DD2<0.05)
        }
        DD32=b[which.max(a)]
        sdDD32=sd(all_snps_no_y[[which.max(a)]]$p_DD2<0.05)
        
        
        
        
        
        b7=numeric()
        for(i in 1:1000){
          a=numeric()
          for(j in 1:50){
            a[j]=all_snps_H_no_y[[j]]$p_DD1[i]
          }
          b7[i]=(sum(a<0.05/50)>0)
        }
        DD41=mean(b7)
        sdDD41=sd(b7) 
        
        b8=numeric()
        for(i in 1:1000){
          a=numeric()
          for(j in 1:50){
            a[j]=all_snps_H_no_y[[j]]$p_DD2[i]
          }
          b8[i]=(sum(a<0.05/50)>0)
        }
        DD42=mean(b8)
        sdDD42=sd(b8) 
        
        
        
        
        
        
        b10=numeric()
        varb10=numeric()
        for(i in 1:length(all_snps_y_only)){
          b10[i]=mean(all_snps_y_only[[i]]$p_DD1<0.05)
          varb10[i]=var(all_snps_y_only[[i]]$p_DD1<0.05)
        }
        DD51=mean(b10)
        sdDD51=sqrt(mean(varb10))
        
        b11t=numeric()
        varb11t=numeric()
        for(i in 1:length(all_snps_y_only)){
          b11t[i]=mean(all_snps_y_only[[i]]$p_DD2<0.05)
          varb11t[i]=var(all_snps_y_only[[i]]$p_DD2<0.05)
        }
        DD52=mean(b11t)
        sdDD52=sqrt(mean(varb11t))
        
        
        
        
        
        b13t=numeric()
        varb13t=numeric()
        for(i in 1:length(all_snps_H_and_y)){
          b13t[i]=mean(all_snps_H_and_y[[i]]$p_DD1<0.05)
          varb13t[i]=var(all_snps_H_and_y[[i]]$p_DD1<0.05)
        }
        DD61=mean(b13t)
        sdDD61=sqrt(mean(varb13t))
        
        
        b14t=numeric()
        varb14t=numeric()
        for(i in 1:length(all_snps_H_and_y)){
          b14t[i]=mean(all_snps_H_and_y[[i]]$p_DD2<0.05)
          varb14t[i]=var(all_snps_H_and_y[[i]]$p_DD2<0.05)
        }
        DD62=mean(b14t)
        sdDD62=sqrt(mean(varb14t))
        
        
        
        
        power1=numeric()
        power21=numeric()
        power22=numeric()
        power_DD21=numeric()
        power_DD22=numeric()
        
        for(i in 1:length(all_snps_y)){
          power1[i]=mean(all_snps_y[[i]]$p_before<0.05)
          power21[i]=mean(all_snps_y[[i]]$p_after1<0.05)
          power22[i]=mean(all_snps_y[[i]]$p_after2<0.05)
          power_DD21[i]=mean(all_snps_y[[i]]$p_DD1<0.05)
          power_DD22[i]=mean(all_snps_y[[i]]$p_DD2<0.05)
        }
        
        change1=power21-power1
        change2=power22-power1
        
        if(sum(change1>0)>0){
          DD71=power_DD21[which.max(change1)]
          sdDD71=sd(all_snps_y[[which.max(change1)]]$p_DD1<0.05)
        }
        if(sum(change1>0)==0){
          DD71=NA
          sdDD71=NA
        }
        if(sum(change2>0)>0){
          DD72=power_DD22[which.max(change2)]
          sdDD72=sd(all_snps_y[[which.max(change2)]]$p_DD2<0.05)
        }
        if(sum(change2>0)==0){
          DD72=NA
          sdDD72=NA
        }
        
        
        if(sum(change1<0)>0){
          DD81=power_DD21[which.min(change1)]
          sdDD81=sd(all_snps_y[[which.min(change1)]]$p_DD1<0.05)
        }
        if(sum(change1<0)==0){
          DD81=NA
          sdDD81=NA
        }
        
        if(sum(change2<0)>0){
          DD82=power_DD22[which.min(change2)]
          sdDD82=sd(all_snps_y[[which.min(change2)]]$p_DD2<0.05)
        }
        if(sum(change2<0)==0){
          DD82=NA
          sdDD82=NA
        }
        
      }
      
      if(10>0){
        b1=numeric()
        varb1=numeric()
        for(i in 1: length(null_snps)){
          b1[i]=mean(null_snps[[i]]$p_SH1<0.05)
          varb1[i]=var(null_snps[[i]]$p_SH1<0.05)
        }
        SH11=mean(b1)
        sdSH11=sqrt(mean(varb1))
        
        
        b2=numeric()
        varb2=numeric()
        for(i in 1: length(null_snps)){
          b2[i]=mean(null_snps[[i]]$p_SH2<0.05)
          varb2[i]=var(null_snps[[i]]$p_SH2<0.05)
        }
        SH12=mean(b2)
        sdSH12=sqrt(mean(varb2))
        
        
        
        
        b4=numeric()
        varb4=numeric()
        for(i in 1: length(all_snps_H_no_y)){
          b4[i]=mean(all_snps_H_no_y[[i]]$p_SH1<0.05)
          varb4[i]=var(all_snps_H_no_y[[i]]$p_SH1<0.05)
        }
        SH21=mean(b4)
        sdSH21=sqrt(mean(varb4))
        
        
        b5=numeric()
        varb5=numeric()
        for(i in 1: length(all_snps_H_no_y)){
          b5[i]=mean(all_snps_H_no_y[[i]]$p_SH2<0.05)
          varb5[i]=var(all_snps_H_no_y[[i]]$p_SH2<0.05)
        }
        SH22=mean(b5)
        sdSH22=sqrt(mean(varb5))
        
        
        
        
        
        a=numeric()
        b=numeric()
        for(i in 1:length(all_snps_no_y)){
          a[i]=mean(all_snps_no_y[[i]]$p_before<0.05)
          b[i]=mean(all_snps_no_y[[i]]$p_SH1<0.05)
        }
        SH31=b[which.max(a)]
        sdSH31=sd(all_snps_no_y[[which.max(a)]]$p_SH1<0.05)
        
        a=numeric()
        b=numeric()
        for(i in 1:length(all_snps_no_y)){
          a[i]=mean(all_snps_no_y[[i]]$p_before<0.05)
          b[i]=mean(all_snps_no_y[[i]]$p_SH2<0.05)
        }
        SH32=b[which.max(a)]
        sdSH32=sd(all_snps_no_y[[which.max(a)]]$p_SH2<0.05)
        
        
        
        
        
        b7=numeric()
        for(i in 1:1000){
          a=numeric()
          for(j in 1:50){
            a[j]=all_snps_H_no_y[[j]]$p_SH1[i]
          }
          b7[i]=(sum(a<0.05/50)>0)
        }
        SH41=mean(b7)
        sdSH41=sd(b7) 
        
        b8=numeric()
        for(i in 1:1000){
          a=numeric()
          for(j in 1:50){
            a[j]=all_snps_H_no_y[[j]]$p_SH2[i]
          }
          b8[i]=(sum(a<0.05/50)>0)
        }
        SH42=mean(b8)
        sdSH42=sd(b8) 
        
        
        
        
        
        
        b10=numeric()
        varb10=numeric()
        for(i in 1:length(all_snps_y_only)){
          b10[i]=mean(all_snps_y_only[[i]]$p_SH1<0.05)
          varb10[i]=var(all_snps_y_only[[i]]$p_SH1<0.05)
        }
        SH51=mean(b10)
        sdSH51=sqrt(mean(varb10))
        
        b11t=numeric()
        varb11t=numeric()
        for(i in 1:length(all_snps_y_only)){
          b11t[i]=mean(all_snps_y_only[[i]]$p_SH2<0.05)
          varb11t[i]=var(all_snps_y_only[[i]]$p_SH2<0.05)
        }
        SH52=mean(b11t)
        sdSH52=sqrt(mean(varb11t))
        
        
        
        
        
        b13t=numeric()
        varb13t=numeric()
        for(i in 1:length(all_snps_H_and_y)){
          b13t[i]=mean(all_snps_H_and_y[[i]]$p_SH1<0.05)
          varb13t[i]=var(all_snps_H_and_y[[i]]$p_SH1<0.05)
        }
        SH61=mean(b13t)
        sdSH61=sqrt(mean(varb13t))
        
        
        b14t=numeric()
        varb14t=numeric()
        for(i in 1:length(all_snps_H_and_y)){
          b14t[i]=mean(all_snps_H_and_y[[i]]$p_SH2<0.05)
          varb14t[i]=var(all_snps_H_and_y[[i]]$p_SH2<0.05)
        }
        SH62=mean(b14t)
        sdSH62=sqrt(mean(varb14t))
        
        
        
        
        power1=numeric()
        power21=numeric()
        power22=numeric()
        power_SH21=numeric()
        power_SH22=numeric()
        
        for(i in 1:length(all_snps_y)){
          power1[i]=mean(all_snps_y[[i]]$p_before<0.05)
          power21[i]=mean(all_snps_y[[i]]$p_after1<0.05)
          power22[i]=mean(all_snps_y[[i]]$p_after2<0.05)
          power_SH21[i]=mean(all_snps_y[[i]]$p_SH1<0.05)
          power_SH22[i]=mean(all_snps_y[[i]]$p_SH2<0.05)
        }
        
        change1=power21-power1
        change2=power22-power1
        
        if(sum(change1>0)>0){
          SH71=power_SH21[which.max(change1)]
          sdSH71=sd(all_snps_y[[which.max(change1)]]$p_SH1<0.05)
        }
        if(sum(change1>0)==0){
          SH71=NA
          sdSH71=NA
        }
        if(sum(change2>0)>0){
          SH72=power_SH22[which.max(change2)]
          sdSH72=sd(all_snps_y[[which.max(change2)]]$p_SH2<0.05)
        }
        if(sum(change2>0)==0){
          SH72=NA
          sdSH72=NA
        }
        
        
        if(sum(change1<0)>0){
          SH81=power_SH21[which.min(change1)]
          sdSH81=sd(all_snps_y[[which.min(change1)]]$p_SH1<0.05)
        }
        if(sum(change1<0)==0){
          SH81=NA
          sdSH81=NA
        }
        
        if(sum(change2<0)>0){
          SH82=power_SH22[which.min(change2)]
          sdSH82=sd(all_snps_y[[which.min(change2)]]$p_SH2<0.05)
        }
        if(sum(change2<0)==0){
          SH82=NA
          sdSH82=NA
        }
        
      }
      
      if(10>0){
        b1=numeric()
        varb1=numeric()
        for(i in 1: length(null_snps)){
          b1[i]=mean(null_snps[[i]]$p_median1<0.05)
          varb1[i]=var(null_snps[[i]]$p_median1<0.05)
        }
        median11=mean(b1)
        sdmedian11=sqrt(mean(varb1))
        
        
        b2=numeric()
        varb2=numeric()
        for(i in 1: length(null_snps)){
          b2[i]=mean(null_snps[[i]]$p_median2<0.05)
          varb2[i]=var(null_snps[[i]]$p_median2<0.05)
        }
        median12=mean(b2)
        sdmedian12=sqrt(mean(varb2))
        
        
        
        
        b4=numeric()
        varb4=numeric()
        for(i in 1: length(all_snps_H_no_y)){
          b4[i]=mean(all_snps_H_no_y[[i]]$p_median1<0.05)
          varb4[i]=var(all_snps_H_no_y[[i]]$p_median1<0.05)
        }
        median21=mean(b4)
        sdmedian21=sqrt(mean(varb4))
        
        
        b5=numeric()
        varb5=numeric()
        for(i in 1: length(all_snps_H_no_y)){
          b5[i]=mean(all_snps_H_no_y[[i]]$p_median2<0.05)
          varb5[i]=var(all_snps_H_no_y[[i]]$p_median2<0.05)
        }
        median22=mean(b5)
        sdmedian22=sqrt(mean(varb5))
        
        
        
        
        
        a=numeric()
        b=numeric()
        for(i in 1:length(all_snps_no_y)){
          a[i]=mean(all_snps_no_y[[i]]$p_before<0.05)
          b[i]=mean(all_snps_no_y[[i]]$p_median1<0.05)
        }
        median31=b[which.max(a)]
        sdmedian31=sd(all_snps_no_y[[which.max(a)]]$p_median1<0.05)
        
        a=numeric()
        b=numeric()
        for(i in 1:length(all_snps_no_y)){
          a[i]=mean(all_snps_no_y[[i]]$p_before<0.05)
          b[i]=mean(all_snps_no_y[[i]]$p_median2<0.05)
        }
        median32=b[which.max(a)]
        sdmedian32=sd(all_snps_no_y[[which.max(a)]]$p_median2<0.05)
        
        
        
        
        
        b7=numeric()
        for(i in 1:1000){
          a=numeric()
          for(j in 1:50){
            a[j]=all_snps_H_no_y[[j]]$p_median1[i]
          }
          b7[i]=(sum(a<0.05/50)>0)
        }
        median41=mean(b7)
        sdmedian41=sd(b7) 
        
        b8=numeric()
        for(i in 1:1000){
          a=numeric()
          for(j in 1:50){
            a[j]=all_snps_H_no_y[[j]]$p_median2[i]
          }
          b8[i]=(sum(a<0.05/50)>0)
        }
        median42=mean(b8)
        sdmedian42=sd(b8) 
        
        
        
        
        
        
        b10=numeric()
        varb10=numeric()
        for(i in 1:length(all_snps_y_only)){
          b10[i]=mean(all_snps_y_only[[i]]$p_median1<0.05)
          varb10[i]=var(all_snps_y_only[[i]]$p_median1<0.05)
        }
        median51=mean(b10)
        sdmedian51=sqrt(mean(varb10))
        
        b11t=numeric()
        varb11t=numeric()
        for(i in 1:length(all_snps_y_only)){
          b11t[i]=mean(all_snps_y_only[[i]]$p_median2<0.05)
          varb11t[i]=var(all_snps_y_only[[i]]$p_median2<0.05)
        }
        median52=mean(b11t)
        sdmedian52=sqrt(mean(varb11t))
        
        
        
        
        
        b13t=numeric()
        varb13t=numeric()
        for(i in 1:length(all_snps_H_and_y)){
          b13t[i]=mean(all_snps_H_and_y[[i]]$p_median1<0.05)
          varb13t[i]=var(all_snps_H_and_y[[i]]$p_median1<0.05)
        }
        median61=mean(b13t)
        sdmedian61=sqrt(mean(varb13t))
        
        
        b14t=numeric()
        varb14t=numeric()
        for(i in 1:length(all_snps_H_and_y)){
          b14t[i]=mean(all_snps_H_and_y[[i]]$p_median2<0.05)
          varb14t[i]=var(all_snps_H_and_y[[i]]$p_median2<0.05)
        }
        median62=mean(b14t)
        sdmedian62=sqrt(mean(varb14t))
        
        
        
        
        power1=numeric()
        power21=numeric()
        power22=numeric()
        power_median21=numeric()
        power_median22=numeric()
        
        for(i in 1:length(all_snps_y)){
          power1[i]=mean(all_snps_y[[i]]$p_before<0.05)
          power21[i]=mean(all_snps_y[[i]]$p_after1<0.05)
          power22[i]=mean(all_snps_y[[i]]$p_after2<0.05)
          power_median21[i]=mean(all_snps_y[[i]]$p_median1<0.05)
          power_median22[i]=mean(all_snps_y[[i]]$p_median2<0.05)
        }
        
        change1=power21-power1
        change2=power22-power1
        
        if(sum(change1>0)>0){
          median71=power_median21[which.max(change1)]
          sdmedian71=sd(all_snps_y[[which.max(change1)]]$p_median1<0.05)
        }
        if(sum(change1>0)==0){
          median71=NA
          sdmedian71=NA
        }
        if(sum(change2>0)>0){
          median72=power_median22[which.max(change2)]
          sdmedian72=sd(all_snps_y[[which.max(change2)]]$p_median2<0.05)
        }
        if(sum(change2>0)==0){
          median72=NA
          sdmedian72=NA
        }
        
        
        if(sum(change1<0)>0){
          median81=power_median21[which.min(change1)]
          sdmedian81=sd(all_snps_y[[which.min(change1)]]$p_median1<0.05)
        }
        if(sum(change1<0)==0){
          median81=NA
          sdmedian81=NA
        }
        
        if(sum(change2<0)>0){
          median82=power_median22[which.min(change2)]
          sdmedian82=sd(all_snps_y[[which.min(change2)]]$p_median2<0.05)
        }
        if(sum(change2<0)==0){
          median82=NA
          sdmedian82=NA
        }
        
      }
      
      if(10>0){
        b1=numeric()
        varb1=numeric()
        for(i in 1: length(null_snps)){
          b1[i]=mean(null_snps[[i]]$p_ivw1<0.05)
          varb1[i]=var(null_snps[[i]]$p_ivw1<0.05)
        }
        ivw11=mean(b1)
        sdivw11=sqrt(mean(varb1))
        
        
        b2=numeric()
        varb2=numeric()
        for(i in 1: length(null_snps)){
          b2[i]=mean(null_snps[[i]]$p_ivw2<0.05)
          varb2[i]=var(null_snps[[i]]$p_ivw2<0.05)
        }
        ivw12=mean(b2)
        sdivw12=sqrt(mean(varb2))
        
        
        
        
        b4=numeric()
        varb4=numeric()
        for(i in 1: length(all_snps_H_no_y)){
          b4[i]=mean(all_snps_H_no_y[[i]]$p_ivw1<0.05)
          varb4[i]=var(all_snps_H_no_y[[i]]$p_ivw1<0.05)
        }
        ivw21=mean(b4)
        sdivw21=sqrt(mean(varb4))
        
        
        b5=numeric()
        varb5=numeric()
        for(i in 1: length(all_snps_H_no_y)){
          b5[i]=mean(all_snps_H_no_y[[i]]$p_ivw2<0.05)
          varb5[i]=var(all_snps_H_no_y[[i]]$p_ivw2<0.05)
        }
        ivw22=mean(b5)
        sdivw22=sqrt(mean(varb5))
        
        
        
        
        
        a=numeric()
        b=numeric()
        for(i in 1:length(all_snps_no_y)){
          a[i]=mean(all_snps_no_y[[i]]$p_before<0.05)
          b[i]=mean(all_snps_no_y[[i]]$p_ivw1<0.05)
        }
        ivw31=b[which.max(a)]
        sdivw31=sd(all_snps_no_y[[which.max(a)]]$p_ivw1<0.05)
        
        a=numeric()
        b=numeric()
        for(i in 1:length(all_snps_no_y)){
          a[i]=mean(all_snps_no_y[[i]]$p_before<0.05)
          b[i]=mean(all_snps_no_y[[i]]$p_ivw2<0.05)
        }
        ivw32=b[which.max(a)]
        sdivw32=sd(all_snps_no_y[[which.max(a)]]$p_ivw2<0.05)
        
        
        
        
        
        b7=numeric()
        for(i in 1:1000){
          a=numeric()
          for(j in 1:50){
            a[j]=all_snps_H_no_y[[j]]$p_ivw1[i]
          }
          b7[i]=(sum(a<0.05/50)>0)
        }
        ivw41=mean(b7)
        sdivw41=sd(b7) 
        
        b8=numeric()
        for(i in 1:1000){
          a=numeric()
          for(j in 1:50){
            a[j]=all_snps_H_no_y[[j]]$p_ivw2[i]
          }
          b8[i]=(sum(a<0.05/50)>0)
        }
        ivw42=mean(b8)
        sdivw42=sd(b8) 
        
        
        
        
        
        
        b10=numeric()
        varb10=numeric()
        for(i in 1:length(all_snps_y_only)){
          b10[i]=mean(all_snps_y_only[[i]]$p_ivw1<0.05)
          varb10[i]=var(all_snps_y_only[[i]]$p_ivw1<0.05)
        }
        ivw51=mean(b10)
        sdivw51=sqrt(mean(varb10))
        
        b11t=numeric()
        varb11t=numeric()
        for(i in 1:length(all_snps_y_only)){
          b11t[i]=mean(all_snps_y_only[[i]]$p_ivw2<0.05)
          varb11t[i]=var(all_snps_y_only[[i]]$p_ivw2<0.05)
        }
        ivw52=mean(b11t)
        sdivw52=sqrt(mean(varb11t))
        
        
        
        
        
        b13t=numeric()
        varb13t=numeric()
        for(i in 1:length(all_snps_H_and_y)){
          b13t[i]=mean(all_snps_H_and_y[[i]]$p_ivw1<0.05)
          varb13t[i]=var(all_snps_H_and_y[[i]]$p_ivw1<0.05)
        }
        ivw61=mean(b13t)
        sdivw61=sqrt(mean(varb13t))
        
        
        b14t=numeric()
        varb14t=numeric()
        for(i in 1:length(all_snps_H_and_y)){
          b14t[i]=mean(all_snps_H_and_y[[i]]$p_ivw2<0.05)
          varb14t[i]=var(all_snps_H_and_y[[i]]$p_ivw2<0.05)
        }
        ivw62=mean(b14t)
        sdivw62=sqrt(mean(varb14t))
        
        
        
        
        power1=numeric()
        power21=numeric()
        power22=numeric()
        power_ivw21=numeric()
        power_ivw22=numeric()
        
        for(i in 1:length(all_snps_y)){
          power1[i]=mean(all_snps_y[[i]]$p_before<0.05)
          power21[i]=mean(all_snps_y[[i]]$p_after1<0.05)
          power22[i]=mean(all_snps_y[[i]]$p_after2<0.05)
          power_ivw21[i]=mean(all_snps_y[[i]]$p_ivw1<0.05)
          power_ivw22[i]=mean(all_snps_y[[i]]$p_ivw2<0.05)
        }
        
        change1=power21-power1
        change2=power22-power1
        
        if(sum(change1>0)>0){
          ivw71=power_ivw21[which.max(change1)]
          sdivw71=sd(all_snps_y[[which.max(change1)]]$p_ivw1<0.05)
        }
        if(sum(change1>0)==0){
          ivw71=NA
          sdivw71=NA
        }
        if(sum(change2>0)>0){
          ivw72=power_ivw22[which.max(change2)]
          sdivw72=sd(all_snps_y[[which.max(change2)]]$p_ivw2<0.05)
        }
        if(sum(change2>0)==0){
          ivw72=NA
          sdivw72=NA
        }
        
        
        if(sum(change1<0)>0){
          ivw81=power_ivw21[which.min(change1)]
          sdivw81=sd(all_snps_y[[which.min(change1)]]$p_ivw1<0.05)
        }
        if(sum(change1<0)==0){
          ivw81=NA
          sdivw81=NA
        }
        
        if(sum(change2<0)>0){
          ivw82=power_ivw22[which.min(change2)]
          sdivw82=sd(all_snps_y[[which.min(change2)]]$p_ivw2<0.05)
        }
        if(sum(change2<0)==0){
          ivw82=NA
          sdivw82=NA
        }
        
      }
      
      if(10>0){
        b1=numeric()
        varb1=numeric()
        for(i in 1: length(null_snps)){
          b1[i]=mean(null_snps[[i]]$p_lasso1<0.05)
          varb1[i]=var(null_snps[[i]]$p_lasso1<0.05)
        }
        lasso11=mean(b1)
        sdlasso11=sqrt(mean(varb1))
        
        
        b2=numeric()
        varb2=numeric()
        for(i in 1: length(null_snps)){
          b2[i]=mean(null_snps[[i]]$p_lasso2<0.05)
          varb2[i]=var(null_snps[[i]]$p_lasso2<0.05)
        }
        lasso12=mean(b2)
        sdlasso12=sqrt(mean(varb2))
        
        
        
        
        b4=numeric()
        varb4=numeric()
        for(i in 1: length(all_snps_H_no_y)){
          b4[i]=mean(all_snps_H_no_y[[i]]$p_lasso1<0.05)
          varb4[i]=var(all_snps_H_no_y[[i]]$p_lasso1<0.05)
        }
        lasso21=mean(b4)
        sdlasso21=sqrt(mean(varb4))
        
        
        b5=numeric()
        varb5=numeric()
        for(i in 1: length(all_snps_H_no_y)){
          b5[i]=mean(all_snps_H_no_y[[i]]$p_lasso2<0.05)
          varb5[i]=var(all_snps_H_no_y[[i]]$p_lasso2<0.05)
        }
        lasso22=mean(b5)
        sdlasso22=sqrt(mean(varb5))
        
        
        
        
        
        a=numeric()
        b=numeric()
        for(i in 1:length(all_snps_no_y)){
          a[i]=mean(all_snps_no_y[[i]]$p_before<0.05)
          b[i]=mean(all_snps_no_y[[i]]$p_lasso1<0.05)
        }
        lasso31=b[which.max(a)]
        sdlasso31=sd(all_snps_no_y[[which.max(a)]]$p_lasso1<0.05)
        
        a=numeric()
        b=numeric()
        for(i in 1:length(all_snps_no_y)){
          a[i]=mean(all_snps_no_y[[i]]$p_before<0.05)
          b[i]=mean(all_snps_no_y[[i]]$p_lasso2<0.05)
        }
        lasso32=b[which.max(a)]
        sdlasso32=sd(all_snps_no_y[[which.max(a)]]$p_lasso2<0.05)
        
        
        
        
        
        b7=numeric()
        for(i in 1:1000){
          a=numeric()
          for(j in 1:50){
            a[j]=all_snps_H_no_y[[j]]$p_lasso1[i]
          }
          b7[i]=(sum(a<0.05/50)>0)
        }
        lasso41=mean(b7)
        sdlasso41=sd(b7) 
        
        b8=numeric()
        for(i in 1:1000){
          a=numeric()
          for(j in 1:50){
            a[j]=all_snps_H_no_y[[j]]$p_lasso2[i]
          }
          b8[i]=(sum(a<0.05/50)>0)
        }
        lasso42=mean(b8)
        sdlasso42=sd(b8) 
        
        
        
        
        
        
        b10=numeric()
        varb10=numeric()
        for(i in 1:length(all_snps_y_only)){
          b10[i]=mean(all_snps_y_only[[i]]$p_lasso1<0.05)
          varb10[i]=var(all_snps_y_only[[i]]$p_lasso1<0.05)
        }
        lasso51=mean(b10)
        sdlasso51=sqrt(mean(varb10))
        
        b11t=numeric()
        varb11t=numeric()
        for(i in 1:length(all_snps_y_only)){
          b11t[i]=mean(all_snps_y_only[[i]]$p_lasso2<0.05)
          varb11t[i]=var(all_snps_y_only[[i]]$p_lasso2<0.05)
        }
        lasso52=mean(b11t)
        sdlasso52=sqrt(mean(varb11t))
        
        
        
        
        
        b13t=numeric()
        varb13t=numeric()
        for(i in 1:length(all_snps_H_and_y)){
          b13t[i]=mean(all_snps_H_and_y[[i]]$p_lasso1<0.05)
          varb13t[i]=var(all_snps_H_and_y[[i]]$p_lasso1<0.05)
        }
        lasso61=mean(b13t)
        sdlasso61=sqrt(mean(varb13t))
        
        
        b14t=numeric()
        varb14t=numeric()
        for(i in 1:length(all_snps_H_and_y)){
          b14t[i]=mean(all_snps_H_and_y[[i]]$p_lasso2<0.05)
          varb14t[i]=var(all_snps_H_and_y[[i]]$p_lasso2<0.05)
        }
        lasso62=mean(b14t)
        sdlasso62=sqrt(mean(varb14t))
        
        
        
        
        power1=numeric()
        power21=numeric()
        power22=numeric()
        power_lasso21=numeric()
        power_lasso22=numeric()
        
        for(i in 1:length(all_snps_y)){
          power1[i]=mean(all_snps_y[[i]]$p_before<0.05)
          power21[i]=mean(all_snps_y[[i]]$p_after1<0.05)
          power22[i]=mean(all_snps_y[[i]]$p_after2<0.05)
          power_lasso21[i]=mean(all_snps_y[[i]]$p_lasso1<0.05)
          power_lasso22[i]=mean(all_snps_y[[i]]$p_lasso2<0.05)
        }
        
        change1=power21-power1
        change2=power22-power1
        
        if(sum(change1>0)>0){
          lasso71=power_lasso21[which.max(change1)]
          sdlasso71=sd(all_snps_y[[which.max(change1)]]$p_lasso1<0.05)
        }
        if(sum(change1>0)==0){
          lasso71=NA
          sdlasso71=NA
        }
        if(sum(change2>0)>0){
          lasso72=power_lasso22[which.max(change2)]
          sdlasso72=sd(all_snps_y[[which.max(change2)]]$p_lasso2<0.05)
        }
        if(sum(change2>0)==0){
          lasso72=NA
          sdlasso72=NA
        }
        
        
        if(sum(change1<0)>0){
          lasso81=power_lasso21[which.min(change1)]
          sdlasso81=sd(all_snps_y[[which.min(change1)]]$p_lasso1<0.05)
        }
        if(sum(change1<0)==0){
          lasso81=NA
          sdlasso81=NA
        }
        
        if(sum(change2<0)>0){
          lasso82=power_lasso22[which.min(change2)]
          sdlasso82=sd(all_snps_y[[which.min(change2)]]$p_lasso2<0.05)
        }
        if(sum(change2<0)==0){
          lasso82=NA
          sdlasso82=NA
        }
        
      }
      
      
      
      if(10>0){
        no=round(c(a11,a21,a31,a41,a51,a61,a71,a81 ),3)
        sdno=round(c(sda11,sda21,sda31,sda41,sda51,sda61,sda71,sda81 ),3)
        cml1=round(c(cml11,cml21,cml31,cml41,cml51,cml61,cml71,cml81 ),3)
        sdcml1=round(c(sdcml11,sdcml21,sdcml31,sdcml41,sdcml51,sdcml61,sdcml71,sdcml81 ),3)
        cml2=round(c(cml12,cml22,cml32,cml42,cml52,cml62,cml72,cml82 ),3)
        sdcml2=round(c(sdcml12,sdcml22,sdcml32,sdcml42,sdcml52,sdcml62,sdcml72,sdcml82 ),3)
        cml3=round(c(cml13,cml23,cml33,cml43,cml53,cml63,cml73,cml83 ),3)
        sdcml3=round(c(sdcml13,sdcml23,sdcml33,sdcml43,sdcml53,sdcml63,sdcml73,sdcml83 ),3)
        
        egger1=round(c(egger11,egger21,egger31,egger41,egger51,egger61,egger71,egger81 ),3)
        sdegger1=round(c(sdegger11,sdegger21,sdegger31,sdegger41,sdegger51,sdegger61,sdegger71,sdegger81 ),3)
        egger2=round(c(egger12,egger22,egger32,egger42,egger52,egger62,egger72,egger82 ),3)
        sdegger2=round(c(sdegger12,sdegger22,sdegger32,sdegger42,sdegger52,sdegger62,sdegger72,sdegger82 ),3)
        
        DD1=round(c(DD11,DD21,DD31,DD41,DD51,DD61,DD71,DD81 ),3)
        sdDD1=round(c(sdDD11,sdDD21,sdDD31,sdDD41,sdDD51,sdDD61,sdDD71,sdDD81 ),3)
        DD2=round(c(DD12,DD22,DD32,DD42,DD52,DD62,DD72,DD82 ),3)
        sdDD2=round(c(sdDD12,sdDD22,sdDD32,sdDD42,sdDD52,sdDD62,sdDD72,sdDD82 ),3)
        
        SH1=round(c(SH11,SH21,SH31,SH41,SH51,SH61,SH71,SH81 ),3)
        sdSH1=round(c(sdSH11,sdSH21,sdSH31,sdSH41,sdSH51,sdSH61,sdSH71,sdSH81 ),3)
        SH2=round(c(SH12,SH22,SH32,SH42,SH52,SH62,SH72,SH82 ),3)
        sdSH2=round(c(sdSH12,sdSH22,sdSH32,sdSH42,sdSH52,sdSH62,sdSH72,sdSH82 ),3)
        
        median1=round(c(median11,median21,median31,median41,median51,median61,median71,median81 ),3)
        sdmedian1=round(c(sdmedian11,sdmedian21,sdmedian31,sdmedian41,sdmedian51,sdmedian61,sdmedian71,sdmedian81 ),3)
        median2=round(c(median12,median22,median32,median42,median52,median62,median72,median82 ),3)
        sdmedian2=round(c(sdmedian12,sdmedian22,sdmedian32,sdmedian42,sdmedian52,sdmedian62,sdmedian72,sdmedian82 ),3)
        
        ivw1=round(c(ivw11,ivw21,ivw31,ivw41,ivw51,ivw61,ivw71,ivw81 ),3)
        sdivw1=round(c(sdivw11,sdivw21,sdivw31,sdivw41,sdivw51,sdivw61,sdivw71,sdivw81 ),3)
        ivw2=round(c(ivw12,ivw22,ivw32,ivw42,ivw52,ivw62,ivw72,ivw82 ),3)
        sdivw2=round(c(sdivw12,sdivw22,sdivw32,sdivw42,sdivw52,sdivw62,sdivw72,sdivw82 ),3)
        
        
        lasso1=round(c(lasso11,lasso21,lasso31,lasso41,lasso51,lasso61,lasso71,lasso81 ),3)
        sdlasso1=round(c(sdlasso11,sdlasso21,sdlasso31,sdlasso41,sdlasso51,sdlasso61,sdlasso71,sdlasso81 ),3)
        lasso2=round(c(lasso12,lasso22,lasso32,lasso42,lasso52,lasso62,lasso72,lasso82 ),3)
        sdlasso2=round(c(sdlasso12,sdlasso22,sdlasso32,sdlasso42,sdlasso52,sdlasso62,sdlasso72,sdlasso82 ),3)
        
      }
      
      
      
      
      
      
      
      result1=data.frame(no,sdno,cml1,sdcml1,cml2,sdcml2,cml3,sdcml3,egger1,sdegger1,egger2,sdegger2,DD1,sdDD1,DD2,sdDD2,SH1,sdSH1,SH2,sdSH2,median1,sdmedian1,median2,sdmedian2,ivw1,sdivw1,ivw2,sdivw2,lasso1,sdlasso1,lasso2,sdlasso2)
      slope_estimated=data.frame(b_cml=mean(all_snps[[1]]$cml_slop1),b_cml_sd=sd(all_snps[[1]]$cml_slop1),b_cml_se=mean(all_snps[[1]]$cml_slop_se1),b_egger=mean(all_snps[[1]]$egger_slop1),b_egger_sd=sd(all_snps[[1]]$egger_slop1),b_egger_se=mean(all_snps[[1]]$egger_slop_se1),b_DD=mean(all_snps[[1]]$DD_slop1),b_DD_sd=sd(all_snps[[1]]$DD_slop1),b_DD_se=mean(all_snps[[1]]$DD_slop_se1),b_SH=mean(all_snps[[1]]$SH_slop1),b_SH_sd=sd(all_snps[[1]]$SH_slop1),b_SH_se=mean(all_snps[[1]]$SH_slop_se1),b_median=mean(all_snps[[1]]$median_slop1),b_median_sd=sd(all_snps[[1]]$median_slop1),b_median_se=mean(all_snps[[1]]$median_slop_se1),b_lasso=mean(all_snps[[1]]$lasso_slop1),b_lasso_sd=sd(all_snps[[1]]$lasso_slop1),b_lasso_se=mean(all_snps[[1]]$lasso_slop_se1),b_ivw=mean(all_snps[[1]]$ivw_slop1),b_ivw_sd=sd(all_snps[[1]]$ivw_slop1),b_ivw_se=mean(all_snps[[1]]$ivw_slop_se1),b_true=mean(all_snps[[1]]$true_slop1))
      slope_estimated=round(slope_estimated,digits=3)
      
      H_only=data.frame(beta_true=rep(0,50),beta_before=rep(0,50),sd_before=rep(0,50),se_before=rep(0,50),error_before=rep(0,50),
                        beta_cml=rep(0,50),sd_cml=rep(0,50),se_cml1=rep(0,50),se_cml2=rep(0,50),se_cml3=rep(0,50),error_cml1=rep(0,50),error_cml2=rep(0,50),error_cml3=rep(0,50),
                        beta_egger=rep(0,50),sd_egger=rep(0,50),se_egger1=rep(0,50),se_egger2=rep(0,50),error_egger1=rep(0,50),error_egger2=rep(0,50),
                        beta_DD=rep(0,50),sd_DD=rep(0,50),se_DD1=rep(0,50),se_DD2=rep(0,50),error_DD1=rep(0,50),error_DD2=rep(0,50),
                        beta_SH=rep(0,50),sd_SH=rep(0,50),se_SH1=rep(0,50),se_SH2=rep(0,50),error_SH1=rep(0,50),error_SH2=rep(0,50),
                        beta_ivw=rep(0,50),sd_ivw=rep(0,50),se_ivw1=rep(0,50),se_ivw2=rep(0,50),error_ivw1=rep(0,50),error_ivw2=rep(0,50),
                        beta_lasso=rep(0,50),sd_lasso=rep(0,50),se_lasso1=rep(0,50),se_lasso2=rep(0,50),error_lasso1=rep(0,50),error_lasso2=rep(0,50))
      
      
      
      for(i1 in 1:50){
        H_only$beta_true[i1]=mean(all_snps[[i1]]$beta_true)
        H_only$beta_before[i1]=mean(all_snps[[i1]]$beta_before)
        H_only$sd_before[i1]=sd(all_snps[[i1]]$beta_before)
        H_only$se_before[i1]=mean(all_snps[[i1]]$se_before)
        H_only$error_before[i1]=mean(all_snps[[i1]]$p_before<0.05)
        H_only$beta_cml[i1]=mean(all_snps[[i1]]$beta_after)
        H_only$sd_cml[i1]=sd(all_snps[[i1]]$beta_after)
        H_only$se_cml1[i1]=mean(all_snps[[i1]]$se_after1)
        H_only$se_cml2[i1]=mean(all_snps[[i1]]$se_after2)
        H_only$se_cml3[i1]=mean(all_snps[[i1]]$se_after3)
        H_only$error_cml1[i1]=mean(all_snps[[i1]]$p_after1<0.05)
        H_only$error_cml2[i1]=mean(all_snps[[i1]]$p_after2<0.05)
        H_only$error_cml3[i1]=mean(all_snps[[i1]]$p_after3<0.05)
        H_only$beta_egger[i1]=mean(all_snps[[i1]]$beta_egger)
        H_only$sd_egger[i1]=sd(all_snps[[i1]]$beta_egger)
        H_only$se_egger1[i1]=mean(all_snps[[i1]]$se_egger1)
        H_only$se_egger2[i1]=mean(all_snps[[i1]]$se_egger2)
        H_only$error_egger1[i1]=mean(all_snps[[i1]]$p_egger1<0.05)
        H_only$error_egger2[i1]=mean(all_snps[[i1]]$p_egger2<0.05)
        H_only$beta_DD[i1]=mean(all_snps[[i1]]$beta_DD)
        H_only$sd_DD[i1]=sd(all_snps[[i1]]$beta_DD)
        H_only$se_DD1[i1]=mean(all_snps[[i1]]$se_DD1)
        H_only$se_DD2[i1]=mean(all_snps[[i1]]$se_DD2)
        H_only$error_DD1[i1]=mean(all_snps[[i1]]$p_DD1<0.05)
        H_only$error_DD2[i1]=mean(all_snps[[i1]]$p_DD2<0.05)
        H_only$beta_SH[i1]=mean(all_snps[[i1]]$beta_SH)
        H_only$sd_SH[i1]=sd(all_snps[[i1]]$beta_SH)
        H_only$se_SH1[i1]=mean(all_snps[[i1]]$se_SH1)
        H_only$se_SH2[i1]=mean(all_snps[[i1]]$se_SH2)
        H_only$error_SH1[i1]=mean(all_snps[[i1]]$p_SH1<0.05)
        H_only$error_SH2[i1]=mean(all_snps[[i1]]$p_SH2<0.05)
        H_only$beta_lasso[i1]=mean(all_snps[[i1]]$beta_lasso)
        H_only$sd_lasso[i1]=sd(all_snps[[i1]]$beta_lasso)
        H_only$se_lasso1[i1]=mean(all_snps[[i1]]$se_lasso1)
        H_only$se_lasso2[i1]=mean(all_snps[[i1]]$se_lasso2)
        H_only$error_lasso1[i1]=mean(all_snps[[i1]]$p_lasso1<0.05)
        H_only$error_lasso2[i1]=mean(all_snps[[i1]]$p_lasso2<0.05)
        H_only$beta_median[i1]=mean(all_snps[[i1]]$beta_median)
        H_only$sd_median[i1]=sd(all_snps[[i1]]$beta_median)
        H_only$se_median1[i1]=mean(all_snps[[i1]]$se_median1)
        H_only$se_median2[i1]=mean(all_snps[[i1]]$se_median2)
        H_only$error_median1[i1]=mean(all_snps[[i1]]$p_median1<0.05)
        H_only$error_median2[i1]=mean(all_snps[[i1]]$p_median2<0.05)
        H_only$beta_ivw[i1]=mean(all_snps[[i1]]$beta_ivw)
        H_only$sd_ivw[i1]=sd(all_snps[[i1]]$beta_ivw)
        H_only$se_ivw1[i1]=mean(all_snps[[i1]]$se_ivw1)
        H_only$se_ivw2[i1]=mean(all_snps[[i1]]$se_ivw2)
        H_only$error_ivw1[i1]=mean(all_snps[[i1]]$p_ivw1<0.05)
        H_only$error_ivw2[i1]=mean(all_snps[[i1]]$p_ivw2<0.05)
      }
      
      H_and_Y=data.frame(beta_true=rep(0,50),beta_before=rep(0,50),sd_before=rep(0,50),se_before=rep(0,50),power_before=rep(0,50),
                         beta_cml=rep(0,50),sd_cml=rep(0,50),se_cml1=rep(0,50),se_cml2=rep(0,50),se_cml3=rep(0,50),power_cml1=rep(0,50),power_cml2=rep(0,50),power_cml3=rep(0,50),
                         beta_egger=rep(0,50),sd_egger=rep(0,50),se_egger1=rep(0,50),se_egger2=rep(0,50),power_egger1=rep(0,50),power_egger2=rep(0,50),
                         beta_DD=rep(0,50),sd_DD=rep(0,50),se_DD1=rep(0,50),se_DD2=rep(0,50),power_DD1=rep(0,50),power_DD2=rep(0,50),
                         beta_SH=rep(0,50),sd_SH=rep(0,50),se_SH1=rep(0,50),se_SH2=rep(0,50),power_SH1=rep(0,50),power_SH2=rep(0,50),
                         beta_ivw=rep(0,50),sd_ivw=rep(0,50),se_ivw1=rep(0,50),se_ivw2=rep(0,50),power_ivw1=rep(0,50),power_ivw2=rep(0,50),
                         beta_lasso=rep(0,50),sd_lasso=rep(0,50),se_lasso1=rep(0,50),se_lasso2=rep(0,50),power_lasso1=rep(0,50),power_lasso2=rep(0,50))
      
      
      for(i1 in 101:150){
        H_and_Y$beta_true[i1-100]=mean(all_snps[[i1]]$beta_true)
        H_and_Y$beta_before[i1-100]=mean(all_snps[[i1]]$beta_before)
        H_and_Y$sd_before[i1-100]=sd(all_snps[[i1]]$beta_before)
        H_and_Y$se_before[i1-100]=mean(all_snps[[i1]]$se_before)
        H_and_Y$power_before[i1-100]=mean(all_snps[[i1]]$p_before<0.05)
        H_and_Y$beta_cml[i1-100]=mean(all_snps[[i1]]$beta_after)
        H_and_Y$sd_cml[i1-100]=sd(all_snps[[i1]]$beta_after)
        H_and_Y$se_cml1[i1-100]=mean(all_snps[[i1]]$se_after1)
        H_and_Y$se_cml2[i1-100]=mean(all_snps[[i1]]$se_after2)
        H_and_Y$se_cml3[i1-100]=mean(all_snps[[i1]]$se_after3)
        H_and_Y$power_cml1[i1-100]=mean(all_snps[[i1]]$p_after1<0.05)
        H_and_Y$power_cml2[i1-100]=mean(all_snps[[i1]]$p_after2<0.05)
        H_and_Y$power_cml3[i1-100]=mean(all_snps[[i1]]$p_after3<0.05)
        H_and_Y$beta_egger[i1-100]=mean(all_snps[[i1]]$beta_egger)
        H_and_Y$sd_egger[i1-100]=sd(all_snps[[i1]]$beta_egger)
        H_and_Y$se_egger1[i1-100]=mean(all_snps[[i1]]$se_egger1)
        H_and_Y$se_egger2[i1-100]=mean(all_snps[[i1]]$se_egger2)
        H_and_Y$power_egger1[i1-100]=mean(all_snps[[i1]]$p_egger1<0.05)
        H_and_Y$power_egger2[i1-100]=mean(all_snps[[i1]]$p_egger2<0.05)
        H_and_Y$beta_DD[i1-100]=mean(all_snps[[i1]]$beta_DD)
        H_and_Y$sd_DD[i1-100]=sd(all_snps[[i1]]$beta_DD)
        H_and_Y$se_DD1[i1-100]=mean(all_snps[[i1]]$se_DD1)
        H_and_Y$se_DD2[i1-100]=mean(all_snps[[i1]]$se_DD2)
        H_and_Y$power_DD1[i1-100]=mean(all_snps[[i1]]$p_DD1<0.05)
        H_and_Y$power_DD2[i1-100]=mean(all_snps[[i1]]$p_DD2<0.05)
        H_and_Y$beta_SH[i1-100]=mean(all_snps[[i1]]$beta_SH)
        H_and_Y$sd_SH[i1-100]=sd(all_snps[[i1]]$beta_SH)
        H_and_Y$se_SH1[i1-100]=mean(all_snps[[i1]]$se_SH1)
        H_and_Y$se_SH2[i1-100]=mean(all_snps[[i1]]$se_SH2)
        H_and_Y$power_SH1[i1-100]=mean(all_snps[[i1]]$p_SH1<0.05)
        H_and_Y$power_SH2[i1-100]=mean(all_snps[[i1]]$p_SH2<0.05)
        H_and_Y$beta_lasso[i1-100]=mean(all_snps[[i1]]$beta_lasso)
        H_and_Y$sd_lasso[i1-100]=sd(all_snps[[i1]]$beta_lasso)
        H_and_Y$se_lasso1[i1-100]=mean(all_snps[[i1]]$se_lasso1)
        H_and_Y$se_lasso2[i1-100]=mean(all_snps[[i1]]$se_lasso2)
        H_and_Y$power_lasso1[i1-100]=mean(all_snps[[i1]]$p_lasso1<0.05)
        H_and_Y$power_lasso2[i1-100]=mean(all_snps[[i1]]$p_lasso2<0.05)
        H_and_Y$beta_median[i1-100]=mean(all_snps[[i1]]$beta_median)
        H_and_Y$sd_median[i1-100]=sd(all_snps[[i1]]$beta_median)
        H_and_Y$se_median1[i1-100]=mean(all_snps[[i1]]$se_median1)
        H_and_Y$se_median2[i1-100]=mean(all_snps[[i1]]$se_median2)
        H_and_Y$power_median1[i1-100]=mean(all_snps[[i1]]$p_median1<0.05)
        H_and_Y$power_median2[i1-100]=mean(all_snps[[i1]]$p_median2<0.05)
        H_and_Y$beta_ivw[i1-100]=mean(all_snps[[i1]]$beta_ivw)
        H_and_Y$sd_ivw[i1-100]=sd(all_snps[[i1]]$beta_ivw)
        H_and_Y$se_ivw1[i1-100]=mean(all_snps[[i1]]$se_ivw1)
        H_and_Y$se_ivw2[i1-100]=mean(all_snps[[i1]]$se_ivw2)
        H_and_Y$power_ivw1[i1-100]=mean(all_snps[[i1]]$p_ivw1<0.05)
        H_and_Y$power_ivw2[i1-100]=mean(all_snps[[i1]]$p_ivw2<0.05)
      }
      
      if(10>0){
        beta_true=numeric()
        beta_before=numeric()
        beta_after=numeric()
        beta_egger=numeric()
        beta_DD=numeric()
        beta_SH=numeric()
        beta_ivw=numeric()
        beta_median=numeric()
        beta_lasso=numeric()
        
        
        sd_before=numeric()
        sd_after=numeric()
        sd_egger=numeric()
        sd_DD=numeric()
        sd_SH=numeric()
        sd_ivw=numeric()
        sd_median=numeric()
        sd_lasso=numeric()
        
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
        x0_DD=numeric()
        y0_DD=numeric()
        x1_DD=numeric()
        y1_DD=numeric()
        x0_SH=numeric()
        y0_SH=numeric()
        x1_SH=numeric()
        y1_SH=numeric()
        x0_median=numeric()
        y0_median=numeric()
        x1_median=numeric()
        y1_median=numeric()
        x0_lasso=numeric()
        y0_lasso=numeric()
        x1_lasso=numeric()
        y1_lasso=numeric()
        x0_ivw=numeric()
        y0_ivw=numeric()
        x1_ivw=numeric()
        y1_ivw=numeric()
        
      }
      
      
      
      for(i in 1:1000){
        beta_before[i]=mean(all_snps[[i]]$beta_before)
        beta_after[i]=mean(all_snps[[i]]$beta_after)
        sd_before[i]=mean(all_snps[[i]]$se_before)
        sd_after[i]=mean(all_snps[[i]]$se_after3)
        beta_true[i]=mean(all_snps[[i]]$beta_true)
        beta_egger[i]=mean(all_snps[[i]]$beta_egger)
        sd_egger[i]=mean(all_snps[[i]]$se_egger2)
        beta_DD[i]=mean(all_snps[[i]]$beta_DD)
        sd_DD[i]=mean(all_snps[[i]]$se_DD1)
        beta_SH[i]=mean(all_snps[[i]]$beta_SH)
        sd_SH[i]=mean(all_snps[[i]]$se_SH1)
        beta_median[i]=mean(all_snps[[i]]$beta_median)
        sd_median[i]=mean(all_snps[[i]]$se_median1)
        beta_lasso[i]=mean(all_snps[[i]]$beta_lasso)
        sd_lasso[i]=mean(all_snps[[i]]$se_lasso1)
        beta_ivw[i]=mean(all_snps[[i]]$beta_ivw)
        sd_ivw[i]=mean(all_snps[[i]]$se_ivw1)
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
        x0_DD[i]=beta_true[i]
        y0_DD[i]=beta_DD[i]-sd_DD[i]
        x1_DD[i]=beta_true[i]
        y1_DD[i]=beta_DD[i]+sd_DD[i]
        x0_SH[i]=beta_true[i]
        y0_SH[i]=beta_SH[i]-sd_SH[i]
        x1_SH[i]=beta_true[i]
        y1_SH[i]=beta_SH[i]+sd_SH[i]
        x0_ivw[i]=beta_true[i]
        y0_ivw[i]=beta_ivw[i]-sd_ivw[i]
        x1_ivw[i]=beta_true[i]
        y1_ivw[i]=beta_ivw[i]+sd_ivw[i]
        x0_median[i]=beta_true[i]
        y0_median[i]=beta_median[i]-sd_median[i]
        x1_median[i]=beta_true[i]
        y1_median[i]=beta_median[i]+sd_median[i]
        x0_lasso[i]=beta_true[i]
        y0_lasso[i]=beta_lasso[i]-sd_lasso[i]
        x1_lasso[i]=beta_true[i]
        y1_lasso[i]=beta_lasso[i]+sd_lasso[i]
      }
      if(10>0){
        jpeg(filename = plotname[1],width=500,height = 500,res=100)
        plot(beta_true,beta_before,xlim=c(min(beta_true),max(beta_true)),ylim=c(min(y0_before,y0_after,y0_egger,y0_DD,y0_SH,y0_ivw,y0_lasso,y0_median),max(y1_before,y1_after,y1_egger,y1_DD,y1_SH,y1_ivw,y1_lasso,y1_median)),main="SNP effects before bias correction",type="n",xlab = "true effect",ylab = "estimated effect")
        
        legend("topleft",pch=c(15,17),col=c("black","darkorchid1"),legend=c("H only","H and Y"),cex=c(0.7,0.7))
        abline(a=0,b=1)
        for(i in 1:50){
          segments(x0=x0_before[i],y0=y0_before[i],x1=x1_before[i],y1=y1_before[i],lwd=0.5,col="gray")
        }
        for(i in 101:150){
          segments(x0=x0_before[i],y0=y0_before[i],x1=x1_before[i],y1=y1_before[i],lwd=0.5,col="gray")
        }
        points(beta_true[1:50],beta_before[1:50],col=1,pch=15,cex=0.7)
        points(beta_true[101:150],beta_before[101:150],col="darkorchid1",pch=17,cex=0.7)
        dev.off()
        
      }
      if(10>0){
        jpeg(filename = plotname[2],width=500,height = 500,res=100)
        
        plot(beta_true,beta_after,xlim=c(min(beta_true),max(beta_true)),ylim=c(min(y0_before,y0_after,y0_egger,y0_DD,y0_SH,y0_ivw,y0_lasso,y0_median),max(y1_before,y1_after,y1_egger,y1_DD,y1_SH,y1_ivw,y1_lasso,y1_median)),main="SNP effects after MVMR-cML bias correction",type="n",xlab = "true effect",ylab = "estimated effect")
        
        legend("topleft",pch=c(15,17),col=c("black","darkorchid1"),legend=c("H only","H and Y"),cex=c(0.7,0.7))
        abline(a=0,b=1)
        for(i in 1:50){
          segments(x0=x0_after[i],y0=y0_after[i],x1=x1_after[i],y1=y1_after[i],lwd=0.5,col="gray")
        }
        for(i in 101:150){
          segments(x0=x0_after[i],y0=y0_after[i],x1=x1_after[i],y1=y1_after[i],lwd=0.5,col="gray")
        }
        
        
        points(beta_true[1:50],beta_after[1:50],col=1,pch=15,cex=0.7)
        points(beta_true[101:150],beta_after[101:150],col="darkorchid1",pch=17,cex=0.7)
        dev.off()
      }
      if(10>0){
        jpeg(filename = plotname[3],width=500,height = 500,res=100)
        
        plot(beta_true,beta_egger,xlim=c(min(beta_true),max(beta_true)),ylim=c(min(y0_before,y0_after,y0_egger,y0_DD,y0_SH,y0_ivw,y0_lasso,y0_median),max(y1_before,y1_after,y1_egger,y1_DD,y1_SH,y1_ivw,y1_lasso,y1_median)),main="SNP effects after MV-Egger bias correction",type="n",xlab = "true effect",ylab = "estimated effect")
        
        legend("topleft",pch=c(15,17),col=c("black","darkorchid1"),legend=c("H only","H and Y"),cex=c(0.7,0.7))
        abline(a=0,b=1)
        for(i in 1:50){
          segments(x0=x0_egger[i],y0=y0_egger[i],x1=x1_egger[i],y1=y1_egger[i],lwd=0.5,col="gray")
        }
        for(i in 101:150){
          segments(x0=x0_egger[i],y0=y0_egger[i],x1=x1_egger[i],y1=y1_egger[i],lwd=0.5,col="gray")
        }
        
        
        points(beta_true[1:50],beta_egger[1:50],col=1,pch=15,cex=0.7)
        points(beta_true[101:150],beta_egger[101:150],col="darkorchid1",pch=17,cex=0.7)
        dev.off()
        
      }
      if(10>0){
        jpeg(filename = plotname[4],width=500,height = 500,res=100)
        
        plot(beta_true,beta_DD,xlim=c(min(beta_true),max(beta_true)),ylim=c(min(y0_before,y0_after,y0_egger,y0_DD,y0_SH,y0_ivw,y0_lasso,y0_median),max(y1_before,y1_after,y1_egger,y1_DD,y1_SH,y1_ivw,y1_lasso,y1_median)),main="SNP effects after DHO bias correction",type="n",xlab = "true effect",ylab = "estimated effect")
        
        legend("topleft",pch=c(15,17),col=c("black","darkorchid1"),legend=c("H only","H and Y"),cex=c(0.7,0.7))
        abline(a=0,b=1)
        for(i in 1:50){
          segments(x0=x0_DD[i],y0=y0_DD[i],x1=x1_DD[i],y1=y1_DD[i],lwd=0.5,col="gray")
        }
        for(i in 101:150){
          segments(x0=x0_DD[i],y0=y0_DD[i],x1=x1_DD[i],y1=y1_DD[i],lwd=0.5,col="gray")
        }
        
        
        points(beta_true[1:50],beta_DD[1:50],col=1,pch=15,cex=0.7)
        points(beta_true[101:150],beta_DD[101:150],col="darkorchid1",pch=17,cex=0.7)
        dev.off()
      }
      if(10>0){
        jpeg(filename = plotname[5],width=500,height = 500,res=100)
        
        plot(beta_true,beta_SH,xlim=c(min(beta_true),max(beta_true)),ylim=c(min(y0_before,y0_after,y0_egger,y0_DD,y0_SH,y0_ivw,y0_lasso,y0_median),max(y1_before,y1_after,y1_egger,y1_DD,y1_SH,y1_ivw,y1_lasso,y1_median)),main="SNP effects after SH bias correction",type="n",xlab = "true effect",ylab = "estimated effect")
        
        legend("topleft",pch=c(15,17),col=c("black","darkorchid1"),legend=c("H only","H and Y"),cex=c(0.7,0.7))
        abline(a=0,b=1)
        for(i in 1:50){
          segments(x0=x0_SH[i],y0=y0_SH[i],x1=x1_SH[i],y1=y1_SH[i],lwd=0.5,col="gray")
        }
        for(i in 101:150){
          segments(x0=x0_SH[i],y0=y0_SH[i],x1=x1_SH[i],y1=y1_SH[i],lwd=0.5,col="gray")
        }
        
        
        points(beta_true[1:50],beta_SH[1:50],col=1,pch=15,cex=0.7)
        points(beta_true[101:150],beta_SH[101:150],col="darkorchid1",pch=17,cex=0.7)
        dev.off()
      }
      if(10>0){
        jpeg(filename = plotname[6],width=500,height = 500,res=100)
        
        plot(beta_true,beta_ivw,xlim=c(min(beta_true),max(beta_true)),ylim=c(min(y0_before,y0_after,y0_egger,y0_DD,y0_ivw,y0_ivw,y0_lasso,y0_median),max(y1_before,y1_after,y1_egger,y1_DD,y1_ivw,y1_ivw,y1_lasso,y1_median)),main="SNP effects after MVMR-IVW bias correction",type="n",xlab = "true effect",ylab = "estimated effect")
        
        legend("topleft",pch=c(15,17),col=c("black","darkorchid1"),legend=c("H only","H and Y"),cex=c(0.7,0.7))
        abline(a=0,b=1)
        for(i in 1:50){
          segments(x0=x0_ivw[i],y0=y0_ivw[i],x1=x1_ivw[i],y1=y1_ivw[i],lwd=0.5,col="gray")
        }
        for(i in 101:150){
          segments(x0=x0_ivw[i],y0=y0_ivw[i],x1=x1_ivw[i],y1=y1_ivw[i],lwd=0.5,col="gray")
        }
        
        
        points(beta_true[1:50],beta_ivw[1:50],col=1,pch=15,cex=0.7)
        points(beta_true[101:150],beta_ivw[101:150],col="darkorchid1",pch=17,cex=0.7)
        dev.off()
      }
      if(10>0){
        jpeg(filename = plotname[7],width=500,height = 500,res=100)
        
        plot(beta_true,beta_lasso,xlim=c(min(beta_true),max(beta_true)),ylim=c(min(y0_before,y0_after,y0_egger,y0_DD,y0_lasso,y0_lasso,y0_lasso,y0_lasso),max(y1_before,y1_after,y1_egger,y1_DD,y1_lasso,y1_lasso,y1_lasso,y1_lasso)),main="SNP effects after MVMR-LASSO bias correction",type="n",xlab = "true effect",ylab = "estimated effect")
        
        legend("topleft",pch=c(15,17),col=c("black","darkorchid1"),legend=c("H only","H and Y"),cex=c(0.7,0.7))
        abline(a=0,b=1)
        for(i in 1:50){
          segments(x0=x0_lasso[i],y0=y0_lasso[i],x1=x1_lasso[i],y1=y1_lasso[i],lwd=0.5,col="gray")
        }
        for(i in 101:150){
          segments(x0=x0_lasso[i],y0=y0_lasso[i],x1=x1_lasso[i],y1=y1_lasso[i],lwd=0.5,col="gray")
        }
        
        
        points(beta_true[1:50],beta_lasso[1:50],col=1,pch=15,cex=0.7)
        points(beta_true[101:150],beta_lasso[101:150],col="darkorchid1",pch=17,cex=0.7)
        dev.off()
        
      }
      if(10>0){
        jpeg(filename = plotname[8],width=500,height = 500,res=100)
        
        plot(beta_true,beta_median,xlim=c(min(beta_true),max(beta_true)),ylim=c(min(y0_before,y0_after,y0_egger,y0_DD,y0_median,y0_median,y0_median,y0_median),max(y1_before,y1_after,y1_egger,y1_DD,y1_median,y1_median,y1_median,y1_median)),main="SNP effects after MVMR-Median bias correction",type="n",xlab = "true effect",ylab = "estimated effect")
        
        legend("topleft",pch=c(15,17),col=c("black","darkorchid1"),legend=c("H only","H and Y"),cex=c(0.7,0.7))
        abline(a=0,b=1)
        for(i in 1:50){
          segments(x0=x0_median[i],y0=y0_median[i],x1=x1_median[i],y1=y1_median[i],lwd=0.5,col="gray")
        }
        for(i in 101:150){
          segments(x0=x0_median[i],y0=y0_median[i],x1=x1_median[i],y1=y1_median[i],lwd=0.5,col="gray")
        }
        
        
        points(beta_true[1:50],beta_median[1:50],col=1,pch=15,cex=0.7)
        points(beta_true[101:150],beta_median[101:150],col="darkorchid1",pch=17,cex=0.7)
        dev.off()
        
      }
      
      
      
      point_est=list(H_only,H_and_Y)
      write.table(result1,file = filename1,sep="\t",col.names = TRUE,quote = FALSE,row.names = FALSE)
      write.table(slope_estimated,file = filename2,sep="\t",col.names = TRUE,quote = FALSE,row.names = FALSE)
      save(point_est,file=filename3)
      
    }
  }
  print(d)
}

