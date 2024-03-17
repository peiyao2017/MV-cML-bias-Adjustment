

for(d in c(2 )){ 
  for(rho in c(0 )){
    if(10>0){
      N=20
      setwd(paste("D:/art sim/New/real_snps_sim/extra_simulation/30invalid/","d=",d,sep=""))
      plotname=numeric()
      name1_1=numeric()
      plotname[1]=paste("D:/art sim/New/real_snps_sim/extra_simulation/30invalid/result/","d=",d,"/","d=",d,"betaXYnot0rho",rho,"_before.jpeg",sep="")
      plotname[2]=paste("D:/art sim/New/real_snps_sim/extra_simulation/30invalid/result/","d=",d,"/","d=",d,"betaXYnot0rho",rho,"_after.jpeg",sep="")
      plotname[3]=paste("D:/art sim/New/real_snps_sim/extra_simulation/30invalid/result/","d=",d,"/","d=",d,"betaXYnot0rho",rho,"_egger.jpeg",sep="")
      
      plotname[4]=paste("D:/art sim/New/real_snps_sim/extra_simulation/30invalid/result/","d=",d,"/","d=",d,"betaXYnot0rho",rho,"_ivw.jpeg",sep="")
      plotname[5]=paste("D:/art sim/New/real_snps_sim/extra_simulation/30invalid/result/","d=",d,"/","d=",d,"betaXYnot0rho",rho,"_lasso.jpeg",sep="")
      plotname[6]=paste("D:/art sim/New/real_snps_sim/extra_simulation/30invalid/result/","d=",d,"/","d=",d,"betaXYnot0rho",rho,"_median.jpeg",sep="")
      filename1=paste("D:/art sim/New/real_snps_sim/extra_simulation/30invalid/result/","d=",d,"/","d=",d,"betaXYnot0rho",rho,".txt",sep="")
      filename2=paste("D:/art sim/New/real_snps_sim/extra_simulation/30invalid/result/","d=",d,"/","d=",d,"betaXYnot0rho",rho,"Slope.txt",sep="")
      filename3=paste("D:/art sim/New/real_snps_sim/extra_simulation/30invalid/result/","d=",d,"/","d=",d,"betaXYnot0rho",rho,"point_effect_estimates.RData",sep="")
      
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
      
      null_snps=list()
      for(i in 1:900){
        null_snps[[i]]=all_snps[[i+100]]
      }
      
      all_snps_no_y=list()
      for(i in c(1:950)){
        if(i<=50){
          all_snps_no_y[[i]]=all_snps[[i]]
        }
        if(i>50){
          all_snps_no_y[[i]]=all_snps[[i+50]]
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
      
      
      all_snps_y=list()
      for(i in 1:50){
        all_snps_y[[i]]=all_snps[[50+i]]
      }
      
      if(10>0){
        a1=matrix(0,nrow=1000,ncol=length(null_snps))
        
        for(i in 1: length(null_snps)){
          a1[,i]=null_snps[[i]]$p_before<0.05
        }
        a11=mean(rowMeans(a1))
        sda11=sd(rowMeans(a1))
        
        a2=matrix(0,nrow=1000,ncol=length(all_snps_H_no_y))
        
        for(i in 1: length(all_snps_H_no_y)){
          a2[,i]=all_snps_H_no_y[[i]]$p_before<0.05
        }
        a21=mean(rowMeans(a2))
        sda21=sd(rowMeans(a2))
        
        
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
        
        a5=matrix(0,nrow=1000,ncol=length(all_snps_y_only))
        for(i in 1: length(all_snps_y_only)){
          a5[,i]=all_snps_y_only[[i]]$p_before<0.05
        }
        a51=mean(rowMeans(a5))
        sda51=sd(rowMeans(a5))
        
        
        
        
        
        
        
        
        
        
        
        power1=numeric()
        power2=numeric()
        
        for(i in 1:length(all_snps_y)){
          power1[i]=mean(all_snps_y[[i]]$p_before<0.05)
          power2[i]=mean(all_snps_y[[i]]$p_after3<0.05)
        }
        change=power2-power1
        if(sum(change>0)>0){
          a61=power1[which.max(change)]
          sda61=sd(all_snps_y[[which.max(change)]]$p_before<0.05)
        }
        if(sum(change>0)==0){
          a61=NA
          sda61=NA
        }
        
        if(sum(change<0)>0){
          a71=power1[which.min(change)]
          sda71=sd(all_snps_y[[which.min(change)]]$p_before<0.05)
        }
        if(sum(change<0)==0){
          a71=NA
          sda71=NA
        }
        
      }
      
      if(10>0){
        b1=matrix(0,nrow=1000,ncol=length(null_snps))
        
        for(i in 1: length(null_snps)){
          b1[,i]=null_snps[[i]]$p_after1<0.05
        }
        cml11=mean(rowMeans(b1))
        sdcml11=sd(rowMeans(b1))
        
        
        b2=matrix(0,nrow=1000,ncol=length(null_snps))
        
        for(i in 1: length(null_snps)){
          b2[,i]=null_snps[[i]]$p_after2<0.05
        }
        cml12=mean(rowMeans(b2))
        sdcml12=sd(rowMeans(b2))
        
        b3=matrix(0,nrow=1000,ncol=length(null_snps))
        
        for(i in 1: length(null_snps)){
          b3[,i]=null_snps[[i]]$p_after3<0.05
        }
        cml13=mean(rowMeans(b3))
        sdcml13=sd(rowMeans(b3))
        
        
        
        
        
        
        b4=matrix(0,nrow=1000,ncol=length(all_snps_H_no_y))
        
        for(i in 1: length(all_snps_H_no_y)){
          b4[,i]= all_snps_H_no_y[[i]]$p_after1<0.05 
          
        }
        cml21=mean(rowMeans(b4))
        sdcml21=sd(rowMeans(b4))
        
        
        b5=matrix(0,nrow=1000,ncol=length(all_snps_H_no_y))
        
        for(i in 1: length(all_snps_H_no_y)){
          b5[,i]= all_snps_H_no_y[[i]]$p_after2<0.05 
          
        }
        cml22=mean(rowMeans(b5))
        sdcml22=sd(rowMeans(b5))
        
        
        b6=matrix(0,nrow=1000,ncol=length(all_snps_H_no_y))
        
        for(i in 1: length(all_snps_H_no_y)){
          b6[,i]= all_snps_H_no_y[[i]]$p_after3<0.05 
          
        }
        cml23=mean(rowMeans(b6))
        sdcml23=sd(rowMeans(b6))
        
        
        
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
        mean(all_snps_no_y[[2]]$p_after2<0.05)
        cml32=b[which.max(a)]
        sdcml32=sd(all_snps_no_y[[2]]$p_after2<0.05)
        
        
        
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
        
        
        b10=matrix(0,nrow=1000,ncol=length(all_snps_y_only))
        for(i in 1: length(all_snps_y_only)){
          b10[,i]=all_snps_y_only[[i]]$p_after1<0.05
        }
        cml51=mean(rowMeans(b10))
        sdcml51=sd(rowMeans(b10))
        
        b11=matrix(0,nrow=1000,ncol=length(all_snps_y_only))
        for(i in 1: length(all_snps_y_only)){
          b11[,i]=all_snps_y_only[[i]]$p_after2<0.05
        }
        cml52=mean(rowMeans(b11))
        sdcml52=sd(rowMeans(b11))
        
        b12=matrix(0,nrow=1000,ncol=length(all_snps_y_only))
        for(i in 1: length(all_snps_y_only)){
          b12[,i]=all_snps_y_only[[i]]$p_after3<0.05
        }
        cml53=mean(rowMeans(b12))
        sdcml53=sd(rowMeans(b12))
        
        
        
        
        
        
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
          cml61=power21[which.max(change1)]
          sdcml61=sd(all_snps_y[[which.max(change1)]]$p_after1<0.05)
        }
        if(sum(change1>0)==0){
          cml61=NA
          sdcml61=NA
        }
        if(sum(change2>0)>0){
          cml62=power22[which.max(change2)]
          sdcml62=sd(all_snps_y[[which.max(change2)]]$p_after2<0.05)
        }
        if(sum(change2>0)==0){
          cml62=NA
          sdcml62=NA
        }
        if(sum(change3>0)>0){
          cml63=power23[which.max(change3)]
          sdcml63=sd(all_snps_y[[which.max(change3)]]$p_after3<0.05)
        }
        if(sum(change3>0)==0){
          cml63=NA
          sdcml63=NA
        }
        
        
        
        if(sum(change1<0)>0){
          cml71=power21[which.min(change1)]
          sdcml71=sd(all_snps_y[[which.min(change1)]]$p_after1<0.05)
        }
        if(sum(change1<0)==0){
          cml71=NA
          sdcml71=NA
        }
        
        if(sum(change2<0)>0){
          cml72=power22[which.min(change2)]
          sdcml72=sd(all_snps_y[[which.min(change2)]]$p_after2<0.05)
        }
        if(sum(change2<0)==0){
          cml72=NA
          sdcml72=NA
        }
        
        
        if(sum(change3<0)>0){
          cml73=power23[which.min(change3)]
          sdcml73=sd(all_snps_y[[which.min(change3)]]$p_after3<0.05)
        }
        if(sum(change3<0)==0){
          cml73=NA
          sdcml73=NA
        }
      }
      
      if(10>0){
        b1=matrix(0,nrow=1000,ncol=length(null_snps))
        
        for(i in 1: length(null_snps)){
          b1[,i]=null_snps[[i]]$p_egger1<0.05
        }
        egger11=mean(rowMeans(b1))
        sdegger11=sd(rowMeans(b1))
        
        
        b2=matrix(0,nrow=1000,ncol=length(null_snps))
        
        for(i in 1: length(null_snps)){
          b2[,i]=null_snps[[i]]$p_egger2<0.05
        }
        egger12=mean(rowMeans(b2))
        sdegger12=sd(rowMeans(b2))
        
        
        
        
        
        
        
        
        b4=matrix(0,nrow=1000,ncol=length(all_snps_H_no_y))
        
        for(i in 1: length(all_snps_H_no_y)){
          b4[,i]= all_snps_H_no_y[[i]]$p_egger1<0.05 
          
        }
        egger21=mean(rowMeans(b4))
        sdegger21=sd(rowMeans(b4))
        
        
        b5=matrix(0,nrow=1000,ncol=length(all_snps_H_no_y))
        
        for(i in 1: length(all_snps_H_no_y)){
          b5[,i]= all_snps_H_no_y[[i]]$p_egger2<0.05 
          
        }
        egger22=mean(rowMeans(b5))
        sdegger22=sd(rowMeans(b5))
        
        
        
        
        
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
        mean(all_snps_no_y[[2]]$p_egger2<0.05)
        egger32=b[which.max(a)]
        sdegger32=sd(all_snps_no_y[[2]]$p_egger2<0.05)
        
        
        
        
        
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
        
        
        
        
        
        b10=matrix(0,nrow=1000,ncol=length(all_snps_y_only))
        for(i in 1: length(all_snps_y_only)){
          b10[,i]=all_snps_y_only[[i]]$p_egger1<0.05
        }
        egger51=mean(rowMeans(b10))
        sdegger51=sd(rowMeans(b10))
        
        b11=matrix(0,nrow=1000,ncol=length(all_snps_y_only))
        for(i in 1: length(all_snps_y_only)){
          b11[,i]=all_snps_y_only[[i]]$p_egger2<0.05
        }
        egger52=mean(rowMeans(b11))
        sdegger52=sd(rowMeans(b11))
        
        
        
        
        
        
        
        
        
        
        power1=numeric()
        power21=numeric()
        power22=numeric()
        
        
        
        for(i in 1:length(all_snps_y)){
          power1[i]=mean(all_snps_y[[i]]$p_before<0.05)
          power21[i]=mean(all_snps_y[[i]]$p_egger1<0.05)
          power22[i]=mean(all_snps_y[[i]]$p_egger2<0.05)
          
        }
        change1=power21-power1
        change2=power22-power1
        
        if(sum(change1>0)>0){
          egger61=power21[which.max(change1)]
          sdegger61=sd(all_snps_y[[which.max(change1)]]$p_egger1<0.05)
        }
        if(sum(change1>0)==0){
          egger61=NA
          sdegger61=NA
        }
        if(sum(change2>0)>0){
          egger62=power22[which.max(change2)]
          sdegger62=sd(all_snps_y[[which.max(change2)]]$p_egger2<0.05)
        }
        if(sum(change2>0)==0){
          egger62=NA
          sdegger62=NA
        }
        
        
        
        
        if(sum(change1<0)>0){
          egger71=power21[which.min(change1)]
          sdegger71=sd(all_snps_y[[which.min(change1)]]$p_egger1<0.05)
        }
        if(sum(change1<0)==0){
          egger71=NA
          sdegger71=NA
        }
        
        if(sum(change2<0)>0){
          egger72=power22[which.min(change2)]
          sdegger72=sd(all_snps_y[[which.min(change2)]]$p_egger2<0.05)
        }
        if(sum(change2<0)==0){
          egger72=NA
          sdegger72=NA
        }
        
        
        
      }
      
      if(10>0){
        b1=matrix(0,nrow=1000,ncol=length(null_snps))
        
        for(i in 1: length(null_snps)){
          b1[,i]=null_snps[[i]]$p_ivw1<0.05
        }
        ivw11=mean(rowMeans(b1))
        sdivw11=sd(rowMeans(b1))
        
        
        b2=matrix(0,nrow=1000,ncol=length(null_snps))
        
        for(i in 1: length(null_snps)){
          b2[,i]=null_snps[[i]]$p_ivw2<0.05
        }
        ivw12=mean(rowMeans(b2))
        sdivw12=sd(rowMeans(b2))
        
        
        
        
        
        
        
        
        b4=matrix(0,nrow=1000,ncol=length(all_snps_H_no_y))
        
        for(i in 1: length(all_snps_H_no_y)){
          b4[,i]= all_snps_H_no_y[[i]]$p_ivw1<0.05 
          
        }
        ivw21=mean(rowMeans(b4))
        sdivw21=sd(rowMeans(b4))
        
        
        b5=matrix(0,nrow=1000,ncol=length(all_snps_H_no_y))
        
        for(i in 1: length(all_snps_H_no_y)){
          b5[,i]= all_snps_H_no_y[[i]]$p_ivw2<0.05 
          
        }
        ivw22=mean(rowMeans(b5))
        sdivw22=sd(rowMeans(b5))
        
        
        
        
        
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
        mean(all_snps_no_y[[2]]$p_ivw2<0.05)
        ivw32=b[which.max(a)]
        sdivw32=sd(all_snps_no_y[[2]]$p_ivw2<0.05)
        
        
        
        
        
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
        
        
        
        
        
        b10=matrix(0,nrow=1000,ncol=length(all_snps_y_only))
        for(i in 1: length(all_snps_y_only)){
          b10[,i]=all_snps_y_only[[i]]$p_ivw1<0.05
        }
        ivw51=mean(rowMeans(b10))
        sdivw51=sd(rowMeans(b10))
        
        b11=matrix(0,nrow=1000,ncol=length(all_snps_y_only))
        for(i in 1: length(all_snps_y_only)){
          b11[,i]=all_snps_y_only[[i]]$p_ivw2<0.05
        }
        ivw52=mean(rowMeans(b11))
        sdivw52=sd(rowMeans(b11))
        
        
        
        
        
        
        
        
        
        power1=numeric()
        power21=numeric()
        power22=numeric()
        
        
        
        for(i in 1:length(all_snps_y)){
          power1[i]=mean(all_snps_y[[i]]$p_before<0.05)
          power21[i]=mean(all_snps_y[[i]]$p_ivw1<0.05)
          power22[i]=mean(all_snps_y[[i]]$p_ivw2<0.05)
          
        }
        change1=power21-power1
        change2=power22-power1
        
        if(sum(change1>0)>0){
          ivw61=power21[which.max(change1)]
          sdivw61=sd(all_snps_y[[which.max(change1)]]$p_ivw1<0.05)
        }
        if(sum(change1>0)==0){
          ivw61=NA
          sdivw61=NA
        }
        if(sum(change2>0)>0){
          ivw62=power22[which.max(change2)]
          sdivw62=sd(all_snps_y[[which.max(change2)]]$p_ivw2<0.05)
        }
        if(sum(change2>0)==0){
          ivw62=NA
          sdivw62=NA
        }
        
        
        
        
        if(sum(change1<0)>0){
          ivw71=power21[which.min(change1)]
          sdivw71=sd(all_snps_y[[which.min(change1)]]$p_ivw1<0.05)
        }
        if(sum(change1<0)==0){
          ivw71=NA
          sdivw71=NA
        }
        
        if(sum(change2<0)>0){
          ivw72=power22[which.min(change2)]
          sdivw72=sd(all_snps_y[[which.min(change2)]]$p_ivw2<0.05)
        }
        if(sum(change2<0)==0){
          ivw72=NA
          sdivw72=NA
        }
        
        
        
      }
      
      if(10>0){
        b1=matrix(0,nrow=1000,ncol=length(null_snps))
        
        for(i in 1: length(null_snps)){
          b1[,i]=null_snps[[i]]$p_median1<0.05
        }
        median11=mean(rowMeans(b1))
        sdmedian11=sd(rowMeans(b1))
        
        
        b2=matrix(0,nrow=1000,ncol=length(null_snps))
        
        for(i in 1: length(null_snps)){
          b2[,i]=null_snps[[i]]$p_median2<0.05
        }
        median12=mean(rowMeans(b2))
        sdmedian12=sd(rowMeans(b2))
        
        
        
        
        
        
        
        
        b4=matrix(0,nrow=1000,ncol=length(all_snps_H_no_y))
        
        for(i in 1: length(all_snps_H_no_y)){
          b4[,i]= all_snps_H_no_y[[i]]$p_median1<0.05 
          
        }
        median21=mean(rowMeans(b4))
        sdmedian21=sd(rowMeans(b4))
        
        
        b5=matrix(0,nrow=1000,ncol=length(all_snps_H_no_y))
        
        for(i in 1: length(all_snps_H_no_y)){
          b5[,i]= all_snps_H_no_y[[i]]$p_median2<0.05 
          
        }
        median22=mean(rowMeans(b5))
        sdmedian22=sd(rowMeans(b5))
        
        
        
        
        
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
        mean(all_snps_no_y[[2]]$p_median2<0.05)
        median32=b[which.max(a)]
        sdmedian32=sd(all_snps_no_y[[2]]$p_median2<0.05)
        
        
        
        
        
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
        
        
        
        
        
        b10=matrix(0,nrow=1000,ncol=length(all_snps_y_only))
        for(i in 1: length(all_snps_y_only)){
          b10[,i]=all_snps_y_only[[i]]$p_median1<0.05
        }
        median51=mean(rowMeans(b10))
        sdmedian51=sd(rowMeans(b10))
        
        b11=matrix(0,nrow=1000,ncol=length(all_snps_y_only))
        for(i in 1: length(all_snps_y_only)){
          b11[,i]=all_snps_y_only[[i]]$p_median2<0.05
        }
        median52=mean(rowMeans(b11))
        sdmedian52=sd(rowMeans(b11))
        
        
        
        
        
        
        
        
        
        
        
        power1=numeric()
        power21=numeric()
        power22=numeric()
        
        
        
        for(i in 1:length(all_snps_y)){
          power1[i]=mean(all_snps_y[[i]]$p_before<0.05)
          power21[i]=mean(all_snps_y[[i]]$p_median1<0.05)
          power22[i]=mean(all_snps_y[[i]]$p_median2<0.05)
          
        }
        change1=power21-power1
        change2=power22-power1
        
        if(sum(change1>0)>0){
          median61=power21[which.max(change1)]
          sdmedian61=sd(all_snps_y[[which.max(change1)]]$p_median1<0.05)
        }
        if(sum(change1>0)==0){
          median61=NA
          sdmedian61=NA
        }
        if(sum(change2>0)>0){
          median62=power22[which.max(change2)]
          sdmedian62=sd(all_snps_y[[which.max(change2)]]$p_median2<0.05)
        }
        if(sum(change2>0)==0){
          median62=NA
          sdmedian62=NA
        }
        
        
        
        
        if(sum(change1<0)>0){
          median71=power21[which.min(change1)]
          sdmedian71=sd(all_snps_y[[which.min(change1)]]$p_median1<0.05)
        }
        if(sum(change1<0)==0){
          median71=NA
          sdmedian71=NA
        }
        
        if(sum(change2<0)>0){
          median72=power22[which.min(change2)]
          sdmedian72=sd(all_snps_y[[which.min(change2)]]$p_median2<0.05)
        }
        if(sum(change2<0)==0){
          median72=NA
          sdmedian72=NA
        }
        
        
        
      }
      
      if(10>0){
        b1=matrix(0,nrow=1000,ncol=length(null_snps))
        
        for(i in 1: length(null_snps)){
          b1[,i]=null_snps[[i]]$p_lasso1<0.05
        }
        lasso11=mean(rowMeans(b1))
        sdlasso11=sd(rowMeans(b1))
        
        
        b2=matrix(0,nrow=1000,ncol=length(null_snps))
        
        for(i in 1: length(null_snps)){
          b2[,i]=null_snps[[i]]$p_lasso2<0.05
        }
        lasso12=mean(rowMeans(b2))
        sdlasso12=sd(rowMeans(b2))
        
        
        
        
        
        
        
        
        b4=matrix(0,nrow=1000,ncol=length(all_snps_H_no_y))
        
        for(i in 1: length(all_snps_H_no_y)){
          b4[,i]= all_snps_H_no_y[[i]]$p_lasso1<0.05 
          
        }
        lasso21=mean(rowMeans(b4))
        sdlasso21=sd(rowMeans(b4))
        
        
        b5=matrix(0,nrow=1000,ncol=length(all_snps_H_no_y))
        
        for(i in 1: length(all_snps_H_no_y)){
          b5[,i]= all_snps_H_no_y[[i]]$p_lasso2<0.05 
          
        }
        lasso22=mean(rowMeans(b5))
        sdlasso22=sd(rowMeans(b5))
        
        
        
        
        
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
        mean(all_snps_no_y[[2]]$p_lasso2<0.05)
        lasso32=b[which.max(a)]
        sdlasso32=sd(all_snps_no_y[[2]]$p_lasso2<0.05)
        
        
        
        
        
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
        
        
        
        
        
        b10=matrix(0,nrow=1000,ncol=length(all_snps_y_only))
        for(i in 1: length(all_snps_y_only)){
          b10[,i]=all_snps_y_only[[i]]$p_lasso1<0.05
        }
        lasso51=mean(rowMeans(b10))
        sdlasso51=sd(rowMeans(b10))
        
        b11=matrix(0,nrow=1000,ncol=length(all_snps_y_only))
        for(i in 1: length(all_snps_y_only)){
          b11[,i]=all_snps_y_only[[i]]$p_lasso2<0.05
        }
        lasso52=mean(rowMeans(b11))
        sdlasso52=sd(rowMeans(b11))
        
        
        
        
        
        
        
        
        
        
        power1=numeric()
        power21=numeric()
        power22=numeric()
        
        
        
        for(i in 1:length(all_snps_y)){
          power1[i]=mean(all_snps_y[[i]]$p_before<0.05)
          power21[i]=mean(all_snps_y[[i]]$p_lasso1<0.05)
          power22[i]=mean(all_snps_y[[i]]$p_lasso2<0.05)
          
        }
        change1=power21-power1
        change2=power22-power1
        
        if(sum(change1>0)>0){
          lasso61=power21[which.max(change1)]
          sdlasso61=sd(all_snps_y[[which.max(change1)]]$p_lasso1<0.05)
        }
        if(sum(change1>0)==0){
          lasso61=NA
          sdlasso61=NA
        }
        if(sum(change2>0)>0){
          lasso62=power22[which.max(change2)]
          sdlasso62=sd(all_snps_y[[which.max(change2)]]$p_lasso2<0.05)
        }
        if(sum(change2>0)==0){
          lasso62=NA
          sdlasso62=NA
        }
        
        
        
        
        if(sum(change1<0)>0){
          lasso71=power21[which.min(change1)]
          sdlasso71=sd(all_snps_y[[which.min(change1)]]$p_lasso1<0.05)
        }
        if(sum(change1<0)==0){
          lasso71=NA
          sdlasso71=NA
        }
        
        if(sum(change2<0)>0){
          lasso72=power22[which.min(change2)]
          sdlasso72=sd(all_snps_y[[which.min(change2)]]$p_lasso2<0.05)
        }
        if(sum(change2<0)==0){
          lasso72=NA
          sdlasso72=NA
        }
        
        
        
      }
      
      
      
      if(10>0){
        no=round(c(a11,a21,a31,a41,a51,a61,a71  ),3)
        sdno=round(c(sda11,sda21,sda31,sda41,sda51,sda61,sda71  ),3)
        cml1=round(c(cml11,cml21,cml31,cml41,cml51,cml61,cml71 ),3)
        sdcml1=round(c(sdcml11,sdcml21,sdcml31,sdcml41,sdcml51,sdcml61,sdcml71  ),3)
        cml2=round(c(cml12,cml22,cml32,cml42,cml52,cml62,cml72  ),3)
        sdcml2=round(c(sdcml12,sdcml22,sdcml32,sdcml42,sdcml52,sdcml62,sdcml72  ),3)
        cml3=round(c(cml13,cml23,cml33,cml43,cml53,cml63,cml73  ),3)
        sdcml3=round(c(sdcml13,sdcml23,sdcml33,sdcml43,sdcml53,sdcml63,sdcml73  ),3)
        
        egger1=round(c(egger11,egger21,egger31,egger41,egger51,egger61,egger71  ),3)
        sdegger1=round(c(sdegger11,sdegger21,sdegger31,sdegger41,sdegger51,sdegger61,sdegger71  ),3)
        egger2=round(c(egger12,egger22,egger32,egger42,egger52,egger62,egger72 ),3)
        sdegger2=round(c(sdegger12,sdegger22,sdegger32,sdegger42,sdegger52,sdegger62,sdegger72 ),3)
        
        
        
        median1=round(c(median11,median21,median31,median41,median51,median61,median71 ),3)
        sdmedian1=round(c(sdmedian11,sdmedian21,sdmedian31,sdmedian41,sdmedian51,sdmedian61,sdmedian71 ),3)
        median2=round(c(median12,median22,median32,median42,median52,median62,median72 ),3)
        sdmedian2=round(c(sdmedian12,sdmedian22,sdmedian32,sdmedian42,sdmedian52,sdmedian62,sdmedian72 ),3)
        
        ivw1=round(c(ivw11,ivw21,ivw31,ivw41,ivw51,ivw61,ivw71 ),3)
        sdivw1=round(c(sdivw11,sdivw21,sdivw31,sdivw41,sdivw51,sdivw61,sdivw71 ),3)
        ivw2=round(c(ivw12,ivw22,ivw32,ivw42,ivw52,ivw62,ivw72 ),3)
        sdivw2=round(c(sdivw12,sdivw22,sdivw32,sdivw42,sdivw52,sdivw62,sdivw72 ),3)
        
        
        lasso1=round(c(lasso11,lasso21,lasso31,lasso41,lasso51,lasso61,lasso71 ),3)
        sdlasso1=round(c(sdlasso11,sdlasso21,sdlasso31,sdlasso41,sdlasso51,sdlasso61,sdlasso71 ),3)
        lasso2=round(c(lasso12,lasso22,lasso32,lasso42,lasso52,lasso62,lasso72 ),3)
        sdlasso2=round(c(sdlasso12,sdlasso22,sdlasso32,sdlasso42,sdlasso52,sdlasso62,sdlasso72 ),3)
        result1=data.frame(no,sdno,cml1,sdcml1,cml2,sdcml2,cml3,sdcml3,egger1,sdegger1,egger2,sdegger2,median1,sdmedian1,median2,sdmedian2,ivw1,sdivw1,ivw2,sdivw2,lasso1,sdlasso1,lasso2,sdlasso2)
      }
      
      
      
      
      
      
      
      
      
      if(10>0){
        H_only=data.frame(beta_true=rep(0,50),beta_before=rep(0,50),sd_before=rep(0,50),se_before=rep(0,50),error_before=rep(0,50),
                          beta_cml=rep(0,50),sd_cml=rep(0,50),se_cml1=rep(0,50),se_cml2=rep(0,50),se_cml3=rep(0,50),error_cml1=rep(0,50),error_cml2=rep(0,50),error_cml3=rep(0,50),
                          beta_egger=rep(0,50),sd_egger=rep(0,50),se_egger1=rep(0,50),se_egger2=rep(0,50),error_egger1=rep(0,50),error_egger2=rep(0,50),
                          beta_ivw=rep(0,50),sd_ivw=rep(0,50),se_ivw1=rep(0,50),se_ivw2=rep(0,50),error_ivw1=rep(0,50),error_ivw2=rep(0,50),
                          beta_lasso=rep(0,50),sd_lasso=rep(0,50),se_lasso1=rep(0,50),se_lasso2=rep(0,50),error_lasso1=rep(0,50),error_lasso2=rep(0,50),
                          beta_median=rep(0,50),sd_median=rep(0,50),se_median1=rep(0,50),se_median2=rep(0,50),error_median1=rep(0,50),error_median2=rep(0,50))
        
        
        
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
        
        Y_only=data.frame(beta_true=rep(0,50),beta_before=rep(0,50),sd_before=rep(0,50),se_before=rep(0,50),power_before=rep(0,50),
                          beta_cml=rep(0,50),sd_cml=rep(0,50),se_cml1=rep(0,50),se_cml2=rep(0,50),se_cml3=rep(0,50),power_cml1=rep(0,50),power_cml2=rep(0,50),power_cml3=rep(0,50),
                          beta_egger=rep(0,50),sd_egger=rep(0,50),se_egger1=rep(0,50),se_egger2=rep(0,50),power_egger1=rep(0,50),power_egger2=rep(0,50),
                          beta_ivw=rep(0,50),sd_ivw=rep(0,50),se_ivw1=rep(0,50),se_ivw2=rep(0,50),power_ivw1=rep(0,50),power_ivw2=rep(0,50),
                          beta_lasso=rep(0,50),sd_lasso=rep(0,50),se_lasso1=rep(0,50),se_lasso2=rep(0,50),power_lasso1=rep(0,50),power_lasso2=rep(0,50),
                          beta_median=rep(0,50),sd_median=rep(0,50),se_median1=rep(0,50),se_median2=rep(0,50),power_median1=rep(0,50),power_median2=rep(0,50))
        
        
        
        for(i1 in 1:50){
          Y_only$beta_true[i1]=mean(all_snps[[i1+50]]$beta_true)
          Y_only$beta_before[i1]=mean(all_snps[[i1+50]]$beta_before)
          Y_only$sd_before[i1]=sd(all_snps[[i1+50]]$beta_before)
          Y_only$se_before[i1]=mean(all_snps[[i1+50]]$se_before)
          Y_only$power_before[i1]=mean(all_snps[[i1+50]]$p_before<0.05)
          Y_only$beta_cml[i1]=mean(all_snps[[i1+50]]$beta_after)
          Y_only$sd_cml[i1]=sd(all_snps[[i1+50]]$beta_after)
          Y_only$se_cml1[i1]=mean(all_snps[[i1+50]]$se_after1)
          Y_only$se_cml2[i1]=mean(all_snps[[i1+50]]$se_after2)
          Y_only$se_cml3[i1]=mean(all_snps[[i1+50]]$se_after3)
          Y_only$power_cml1[i1]=mean(all_snps[[i1+50]]$p_after1<0.05)
          Y_only$power_cml2[i1]=mean(all_snps[[i1+50]]$p_after2<0.05)
          Y_only$power_cml3[i1]=mean(all_snps[[i1+50]]$p_after3<0.05)
          Y_only$beta_egger[i1]=mean(all_snps[[i1+50]]$beta_egger)
          Y_only$sd_egger[i1]=sd(all_snps[[i1+50]]$beta_egger)
          Y_only$se_egger1[i1]=mean(all_snps[[i1+50]]$se_egger1)
          Y_only$se_egger2[i1]=mean(all_snps[[i1+50]]$se_egger2)
          Y_only$power_egger1[i1]=mean(all_snps[[i1+50]]$p_egger1<0.05)
          Y_only$power_egger2[i1]=mean(all_snps[[i1+50]]$p_egger2<0.05)
          
          Y_only$beta_lasso[i1]=mean(all_snps[[i1+50]]$beta_lasso)
          Y_only$sd_lasso[i1]=sd(all_snps[[i1+50]]$beta_lasso)
          Y_only$se_lasso1[i1]=mean(all_snps[[i1+50]]$se_lasso1)
          Y_only$se_lasso2[i1]=mean(all_snps[[i1+50]]$se_lasso2)
          Y_only$power_lasso1[i1]=mean(all_snps[[i1+50]]$p_lasso1<0.05)
          Y_only$power_lasso2[i1]=mean(all_snps[[i1+50]]$p_lasso2<0.05)
          Y_only$beta_median[i1]=mean(all_snps[[i1+50]]$beta_median)
          Y_only$sd_median[i1]=sd(all_snps[[i1+50]]$beta_median)
          Y_only$se_median1[i1]=mean(all_snps[[i1+50]]$se_median1)
          Y_only$se_median2[i1]=mean(all_snps[[i1+50]]$se_median2)
          Y_only$power_median1[i1]=mean(all_snps[[i1+50]]$p_median1<0.05)
          Y_only$power_median2[i1]=mean(all_snps[[i1+50]]$p_median2<0.05)
          Y_only$beta_ivw[i1]=mean(all_snps[[i1+50]]$beta_ivw)
          Y_only$sd_ivw[i1]=sd(all_snps[[i1+50]]$beta_ivw)
          Y_only$se_ivw1[i1]=mean(all_snps[[i1+50]]$se_ivw1)
          Y_only$se_ivw2[i1]=mean(all_snps[[i1+50]]$se_ivw2)
          Y_only$power_ivw1[i1]=mean(all_snps[[i1+50]]$p_ivw1<0.05)
          Y_only$power_ivw2[i1]=mean(all_snps[[i1+50]]$p_ivw2<0.05)
        }
        
        Null_snps=data.frame(beta_true=rep(0,900),beta_before=rep(0,900),sd_before=rep(0,900),se_before=rep(0,900),error_before=rep(0,900),
                             beta_cml=rep(0,900),sd_cml=rep(0,900),se_cml1=rep(0,900),se_cml2=rep(0,900),se_cml3=rep(0,900),error_cml1=rep(0,900),error_cml2=rep(0,900),error_cml3=rep(0,900),
                             beta_egger=rep(0,900),sd_egger=rep(0,900),se_egger1=rep(0,900),se_egger2=rep(0,900),error_egger1=rep(0,900),error_egger2=rep(0,900),
                             beta_ivw=rep(0,900),sd_ivw=rep(0,900),se_ivw1=rep(0,900),se_ivw2=rep(0,900),error_ivw1=rep(0,900),error_ivw2=rep(0,900),
                             beta_lasso=rep(0,900),sd_lasso=rep(0,900),se_lasso1=rep(0,900),se_lasso2=rep(0,900),error_lasso1=rep(0,900),error_lasso2=rep(0,900),
                             beta_median=rep(0,900),sd_median=rep(0,900),se_median1=rep(0,900),se_median2=rep(0,900),error_median1=rep(0,900),error_median2=rep(0,900))
        
        
        
        for(i1 in 1:900){
          Null_snps$beta_true[i1]=mean(all_snps[[i1+100]]$beta_true)
          Null_snps$beta_before[i1]=mean(all_snps[[i1+100]]$beta_before)
          Null_snps$sd_before[i1]=sd(all_snps[[i1+100]]$beta_before)
          Null_snps$se_before[i1]=mean(all_snps[[i1+100]]$se_before)
          Null_snps$error_before[i1]=mean(all_snps[[i1+100]]$p_before<0.05)
          Null_snps$beta_cml[i1]=mean(all_snps[[i1+100]]$beta_after)
          Null_snps$sd_cml[i1]=sd(all_snps[[i1+100]]$beta_after)
          Null_snps$se_cml1[i1]=mean(all_snps[[i1+100]]$se_after1)
          Null_snps$se_cml2[i1]=mean(all_snps[[i1+100]]$se_after2)
          Null_snps$se_cml3[i1]=mean(all_snps[[i1+100]]$se_after3)
          Null_snps$error_cml1[i1]=mean(all_snps[[i1+100]]$p_after1<0.05)
          Null_snps$error_cml2[i1]=mean(all_snps[[i1+100]]$p_after2<0.05)
          Null_snps$error_cml3[i1]=mean(all_snps[[i1+100]]$p_after3<0.05)
          Null_snps$beta_egger[i1]=mean(all_snps[[i1+100]]$beta_egger)
          Null_snps$sd_egger[i1]=sd(all_snps[[i1+100]]$beta_egger)
          Null_snps$se_egger1[i1]=mean(all_snps[[i1+100]]$se_egger1)
          Null_snps$se_egger2[i1]=mean(all_snps[[i1+100]]$se_egger2)
          Null_snps$error_egger1[i1]=mean(all_snps[[i1+100]]$p_egger1<0.05)
          Null_snps$error_egger2[i1]=mean(all_snps[[i1+100]]$p_egger2<0.05)
          
          Null_snps$beta_lasso[i1]=mean(all_snps[[i1+100]]$beta_lasso)
          Null_snps$sd_lasso[i1]=sd(all_snps[[i1+100]]$beta_lasso)
          Null_snps$se_lasso1[i1]=mean(all_snps[[i1+100]]$se_lasso1)
          Null_snps$se_lasso2[i1]=mean(all_snps[[i1+100]]$se_lasso2)
          Null_snps$error_lasso1[i1]=mean(all_snps[[i1+100]]$p_lasso1<0.05)
          Null_snps$error_lasso2[i1]=mean(all_snps[[i1+100]]$p_lasso2<0.05)
          Null_snps$beta_median[i1]=mean(all_snps[[i1+100]]$beta_median)
          Null_snps$sd_median[i1]=sd(all_snps[[i1+100]]$beta_median)
          Null_snps$se_median1[i1]=mean(all_snps[[i1+100]]$se_median1)
          Null_snps$se_median2[i1]=mean(all_snps[[i1+100]]$se_median2)
          Null_snps$error_median1[i1]=mean(all_snps[[i1+100]]$p_median1<0.05)
          Null_snps$error_median2[i1]=mean(all_snps[[i1+100]]$p_median2<0.05)
          Null_snps$beta_ivw[i1]=mean(all_snps[[i1+100]]$beta_ivw)
          Null_snps$sd_ivw[i1]=sd(all_snps[[i1+100]]$beta_ivw)
          Null_snps$se_ivw1[i1]=mean(all_snps[[i1+100]]$se_ivw1)
          Null_snps$se_ivw2[i1]=mean(all_snps[[i1+100]]$se_ivw2)
          Null_snps$error_ivw1[i1]=mean(all_snps[[i1+100]]$p_ivw1<0.05)
          Null_snps$error_ivw2[i1]=mean(all_snps[[i1+100]]$p_ivw2<0.05)
        }
        
        b_true_name=numeric()
        b_egger_name=numeric()
        b_cml_name=numeric()
        b_ivw_name=numeric()
        b_lasso_name=numeric()
        b_median_name=numeric()
        se_b_true_name=numeric()
        se_b_egger_name=numeric()
        se_b_cml_name=numeric()
        se_b_ivw_name=numeric()
        se_b_lasso_name=numeric()
        se_b_median_name=numeric()
        for(i1 in 1:d){
          b_true_name[i1]=paste("true_slop",i1,sep="")
          b_egger_name[i1]=paste("egger_slop",i1,sep="")
          b_cml_name[i1]=paste("cml_slop",i1,sep="")
          b_ivw_name[i1]=paste("ivw_slop",i1,sep="")
          b_lasso_name[i1]=paste("lasso_slop",i1,sep="")
          b_median_name[i1]=paste("median_slop",i1,sep="")
          se_b_true_name[i1]=paste("true_slop_se",i1,sep="")
          se_b_egger_name[i1]=paste("egger_slop_se",i1,sep="")
          se_b_cml_name[i1]=paste("cml_slop_se",i1,sep="")
          se_b_ivw_name[i1]=paste("ivw_slop_se",i1,sep="")
          se_b_lasso_name[i1]=paste("lasso_slop_se",i1,sep="")
          se_b_median_name[i1]=paste("median_slop_se",i1,sep="")
        }
        
        slope_estimated=data.frame(b_cml=sapply(all_snps[[1]][b_cml_name],mean),b_cml_sd=sapply(all_snps[[1]][b_cml_name],sd),b_cml_se= sapply(all_snps[[1]][se_b_cml_name],mean ),
                                   b_egger=sapply(all_snps[[1]][b_egger_name],mean),b_egger_sd=sapply(all_snps[[1]][b_egger_name],sd),b_egger_se= sapply(all_snps[[1]][se_b_egger_name],mean ),
                                   b_ivw=sapply(all_snps[[1]][b_ivw_name],mean),b_ivw_sd=sapply(all_snps[[1]][b_ivw_name],sd),b_ivw_se= sapply(all_snps[[1]][se_b_ivw_name],mean) ,
                                   b_lasso=sapply(all_snps[[1]][b_lasso_name],mean),b_lasso_sd=sapply(all_snps[[1]][b_lasso_name],sd),b_lasso_se= sapply(all_snps[[1]][se_b_lasso_name],mean ),
                                   b_median=sapply(all_snps[[1]][b_median_name],mean),b_median_sd=sapply(all_snps[[1]][b_median_name],sd),b_median_se= sapply(all_snps[[1]][se_b_median_name],mean ), 
                                   b_true=sapply(all_snps[[1]][b_true_name],mean),b_true_sd=sapply(all_snps[[1]][b_true_name],sd),b_true_se= sapply(all_snps[[1]][se_b_true_name],mean ))
        slope_estimated=round(slope_estimated,digits=3)
        
      }
      
      
      
      
      
      
      
      
      if(10>0){
        beta_true=numeric()
        beta_before=numeric()
        beta_after=numeric()
        beta_egger=numeric()
        
        beta_ivw=numeric()
        beta_median=numeric()
        beta_lasso=numeric()
        
        
        sd_before=numeric()
        sd_after=numeric()
        sd_egger=numeric()
        
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
        plot(beta_true,beta_before,xlim=c(min(beta_true),max(beta_true)),ylim=c(min(y0_before,y0_after,y0_egger,y0_ivw,y0_lasso,y0_median),max(y1_before,y1_after,y1_egger,y1_ivw,y1_lasso,y1_median)),main="SNP effects before bias correction",type="n",xlab = "true effect",ylab = "estimated effect")
        
        
        abline(a=0,b=1)
        for(i in 50:100){
          segments(x0=x0_before[i],y0=y0_before[i],x1=x1_before[i],y1=y1_before[i],lwd=0.5,col="gray")
        }
        
        
        points(beta_true[50:100],beta_before[50:100],col="darkorchid1",pch=17,cex=0.7)
        dev.off()
        
      }
      if(10>0){
        jpeg(filename = plotname[2],width=500,height = 500,res=100)
        
        plot(beta_true,beta_after,xlim=c(min(beta_true),max(beta_true)),ylim=c(min(y0_before,y0_after,y0_egger,y0_ivw,y0_lasso,y0_median),max(y1_before,y1_after,y1_egger,y1_ivw,y1_lasso,y1_median)),main="SNP effects after MVMR-cML bias correction",type="n",xlab = "true effect",ylab = "estimated effect")
        
        
        abline(a=0,b=1)
        
        for(i in 50:100){
          segments(x0=x0_after[i],y0=y0_after[i],x1=x1_after[i],y1=y1_after[i],lwd=0.5,col="gray")
        }
        
        
        
        points(beta_true[50:100],beta_after[50:100],col="darkorchid1",pch=17,cex=0.7)
        dev.off()
      }
      if(10>0){
        jpeg(filename = plotname[3],width=500,height = 500,res=100)
        
        plot(beta_true,beta_egger,xlim=c(min(beta_true),max(beta_true)),ylim=c(min(y0_before,y0_after,y0_egger,y0_ivw,y0_lasso,y0_median),max(y1_before,y1_after,y1_egger,y1_ivw,y1_lasso,y1_median)),main="SNP effects after MV-Egger bias correction",type="n",xlab = "true effect",ylab = "estimated effect")
        
        
        abline(a=0,b=1)
        
        for(i in 50:100){
          segments(x0=x0_egger[i],y0=y0_egger[i],x1=x1_egger[i],y1=y1_egger[i],lwd=0.5,col="gray")
        }
        
        
        
        points(beta_true[50:100],beta_egger[50:100],col="darkorchid1",pch=17,cex=0.7)
        dev.off()
        
      }
      
      
      if(10>0){
        jpeg(filename = plotname[4],width=500,height = 500,res=100)
        
        plot(beta_true,beta_ivw,xlim=c(min(beta_true),max(beta_true)),ylim=c(min(y0_before,y0_after,y0_egger,y0_ivw,y0_ivw,y0_lasso,y0_median),max(y1_before,y1_after,y1_egger,y1_ivw,y1_ivw,y1_lasso,y1_median)),main="SNP effects after MVMR-IVW bias correction",type="n",xlab = "true effect",ylab = "estimated effect")
        
        
        abline(a=0,b=1)
        
        for(i in 50:100){
          segments(x0=x0_ivw[i],y0=y0_ivw[i],x1=x1_ivw[i],y1=y1_ivw[i],lwd=0.5,col="gray")
        }
        
        
        
        points(beta_true[50:100],beta_ivw[50:100],col="darkorchid1",pch=17,cex=0.7)
        dev.off()
      }
      if(10>0){
        jpeg(filename = plotname[5],width=500,height = 500,res=100)
        
        plot(beta_true,beta_lasso,xlim=c(min(beta_true),max(beta_true)),ylim=c(min(y0_before,y0_after,y0_egger,y0_lasso,y0_lasso,y0_lasso,y0_lasso),max(y1_before,y1_after,y1_egger,y1_lasso,y1_lasso,y1_lasso,y1_lasso)),main="SNP effects after MVMR-LASSO bias correction",type="n",xlab = "true effect",ylab = "estimated effect")
        
        
        abline(a=0,b=1)
        
        for(i in 50:100){
          segments(x0=x0_lasso[i],y0=y0_lasso[i],x1=x1_lasso[i],y1=y1_lasso[i],lwd=0.5,col="gray")
        }
        
        
        
        points(beta_true[50:100],beta_lasso[50:100],col="darkorchid1",pch=17,cex=0.7)
        dev.off()
        
      }
      if(10>0){
        jpeg(filename = plotname[6],width=500,height = 500,res=100)
        
        plot(beta_true,beta_median,xlim=c(min(beta_true),max(beta_true)),ylim=c(min(y0_before,y0_after,y0_egger,y0_median,y0_median,y0_median,y0_median),max(y1_before,y1_after,y1_egger,y1_median,y1_median,y1_median,y1_median)),main="SNP effects after MVMR-Median bias correction",type="n",xlab = "true effect",ylab = "estimated effect")
        
        
        abline(a=0,b=1)
        
        for(i in 50:100){
          segments(x0=x0_median[i],y0=y0_median[i],x1=x1_median[i],y1=y1_median[i],lwd=0.5,col="gray")
        }
        
        
        
        points(beta_true[50:100],beta_median[50:100],col="darkorchid1",pch=17,cex=0.7)
        dev.off()
        
      }
      
      point_est=list(H_only,Y_only,Null )
      write.table(result1,file = filename1,sep="\t",col.names = TRUE,quote = FALSE,row.names = FALSE)
      write.table(slope_estimated,file = filename2,sep="\t",col.names = TRUE,quote = FALSE,row.names = FALSE)
      save(point_est,file=filename3)
      
    }
  }
  print(d)
}


