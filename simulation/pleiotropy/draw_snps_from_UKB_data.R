library(BEDMatrix)
setwd("/panfs/jay/groups/20/panwei/wan01299/collider_bias/BMI/MV_Egger/new/real_snp_sim/pleiotropy/")
nIV=30
Nsnps=1000
Nsnpsz=3000
invalid=0.5
Gmatrixtotal=BEDMatrix("/panfs/jay/groups/20/panwei/wan01299/collider_bias/BMI/MV_Egger/new/new_gwas/110K_QCed0.1.bed",simple_names = TRUE)
Indsnps=read.table("/panfs/jay/groups/20/panwei/wan01299/collider_bias/BMI/MV_Egger/new/new_gwas/110K_QCed0.001.bim",header = FALSE)
snps=read.table("/panfs/jay/groups/20/panwei/wan01299/collider_bias/BMI/MV_Egger/new/new_gwas/110K_QCed0.1.bim",header = FALSE)
pos=read.table("/panfs/jay/groups/20/panwei/wan01299/collider_bias/BMI/MV_Egger/new/real_snp_sim/pleiotropy/fourier_ls-all.bed",sep="\t",header=TRUE)


used_snps=snps[snps$V1==1,]

used_pos=pos[pos$chr=="chr1 ",]
block_list=list()
k=0
for(i in 1:nrow(used_pos)){
  start1=used_pos$start[i]
  stop1=used_pos$stop[i]
  a=numeric()
  for(j in 1:nrow(used_snps)){
    if(used_snps$V4[j]>=start1&used_snps$V4[j]<=stop1){
      a=c(a,used_snps$V2[j])
    }
    
  }
  if(length(a)>0){
    k=k+1
    block_list[[k]]=a
  }
}

a1=c(block_list[[1]],block_list[[2]],block_list[[3]],block_list[[4]], block_list[[5]],block_list[[6]],block_list[[7]],block_list[[8]],block_list[[9]], block_list[[10]])
a2=c(block_list[[length(block_list)]],block_list[[length(block_list)-1]],block_list[[length(block_list)-2]],block_list[[length(block_list)-3]],block_list[[length(block_list)-4]],block_list[[length(block_list)-5]],block_list[[length(block_list)-6]],block_list[[length(block_list)-7]],block_list[[length(block_list)-8]],block_list[[length(block_list)-9]])

non_null_snps=sample(a1,size = 0.15*Nsnps-nIV,replace = FALSE)
IVID=sample( intersect(Indsnps$V2, a1[!a1%in%non_null_snps])  ,size=nIV,replace = FALSE)

null_snps=sample(a2,size = 0.85*Nsnps,replace = FALSE)

Null_z=sample(intersect(Indsnps$V2[Indsnps$V1!=1  ],snps$V2),size=Nsnpsz,replace = FALSE)
 

vp=1-invalid
X_only=c(non_null_snps[1:(Nsnps*0.05-nIV*vp)],IVID[1:(nIV*vp)])
n1=(Nsnps*0.05-nIV*vp)+1
n2=(Nsnps*0.05-nIV*vp)+Nsnps*0.05
Y_only=non_null_snps[n1:n2]
n3=n2+1
n4=n2+Nsnps*0.05-nIV*(1-vp)
n5=nIV*vp+1
n6=nIV
X_and_Y=c(non_null_snps[n3:n4],IVID[n5:n6])
 

simID=c(X_only,Y_only,X_and_Y,null_snps,Null_z )
simIV=IVID
Gmatrix=Gmatrixtotal[,colnames(Gmatrixtotal)%in%simID]
for(i in 1:ncol(Gmatrix)){
  Gmatrix[,i][is.na( Gmatrix[,i])]=mean(Gmatrix[,i],na.rm=TRUE)
}

Geno=list(Gmatrix,simID,simIV)
save(Geno,file = "Geno50.RDATA")


library(BEDMatrix)
setwd("/panfs/jay/groups/20/panwei/wan01299/collider_bias/BMI/MV_Egger/new/real_snp_sim/pleiotropy/")
nIV=30
Nsnps=1000
Nsnpsz=3000
invalid=0.3
Gmatrixtotal=BEDMatrix("/panfs/jay/groups/20/panwei/wan01299/collider_bias/BMI/MV_Egger/new/new_gwas/110K_QCed0.1.bed",simple_names = TRUE)
Indsnps=read.table("/panfs/jay/groups/20/panwei/wan01299/collider_bias/BMI/MV_Egger/new/new_gwas/110K_QCed0.001.bim",header = FALSE)
snps=read.table("/panfs/jay/groups/20/panwei/wan01299/collider_bias/BMI/MV_Egger/new/new_gwas/110K_QCed0.1.bim",header = FALSE)
pos=read.table("/panfs/jay/groups/20/panwei/wan01299/collider_bias/BMI/MV_Egger/new/real_snp_sim/pleiotropy/fourier_ls-all.bed",sep="\t",header=TRUE)


used_snps=snps[snps$V1==1,]

used_pos=pos[pos$chr=="chr1 ",]
block_list=list()
k=0
for(i in 1:nrow(used_pos)){
  start1=used_pos$start[i]
  stop1=used_pos$stop[i]
  a=numeric()
  for(j in 1:nrow(used_snps)){
    if(used_snps$V4[j]>=start1&used_snps$V4[j]<=stop1){
      a=c(a,used_snps$V2[j])
    }
    
  }
  if(length(a)>0){
    k=k+1
    block_list[[k]]=a
  }
}

a1=c(block_list[[1]],block_list[[2]],block_list[[3]],block_list[[4]], block_list[[5]],block_list[[6]],block_list[[7]],block_list[[8]],block_list[[9]], block_list[[10]])
a2=c(block_list[[length(block_list)]],block_list[[length(block_list)-1]],block_list[[length(block_list)-2]],block_list[[length(block_list)-3]],block_list[[length(block_list)-4]],block_list[[length(block_list)-5]],block_list[[length(block_list)-6]],block_list[[length(block_list)-7]],block_list[[length(block_list)-8]],block_list[[length(block_list)-9]])

non_null_snps=sample(a1,size = 0.15*Nsnps-nIV,replace = FALSE)
IVID=sample( intersect(Indsnps$V2, a1[!a1%in%non_null_snps])  ,size=nIV,replace = FALSE)

null_snps=sample(a2,size = 0.85*Nsnps,replace = FALSE)

Null_z=sample(intersect(Indsnps$V2[Indsnps$V1!=1  ],snps$V2),size=Nsnpsz,replace = FALSE)
 

vp=1-invalid
X_only=c(non_null_snps[1:(Nsnps*0.05-nIV*vp)],IVID[1:(nIV*vp)])
n1=(Nsnps*0.05-nIV*vp)+1
n2=(Nsnps*0.05-nIV*vp)+Nsnps*0.05
Y_only=non_null_snps[n1:n2]
n3=n2+1
n4=n2+Nsnps*0.05-nIV*(1-vp)
n5=nIV*vp+1
n6=nIV
X_and_Y=c(non_null_snps[n3:n4],IVID[n5:n6])
 

simID=c(X_only,Y_only,X_and_Y,null_snps,Null_z )
simIV=IVID
Gmatrix=Gmatrixtotal[,colnames(Gmatrixtotal)%in%simID]
for(i in 1:ncol(Gmatrix)){
  Gmatrix[,i][is.na( Gmatrix[,i])]=mean(Gmatrix[,i],na.rm=TRUE)
}

Geno=list(Gmatrix,simID,simIV)
save(Geno,file = "Geno30.RDATA")





library(BEDMatrix)
setwd("/panfs/jay/groups/20/panwei/wan01299/collider_bias/BMI/MV_Egger/new/real_snp_sim/pleiotropy/")
nIV=30
Nsnps=1000
Nsnpsz=3000
invalid=0 
Gmatrixtotal=BEDMatrix("/panfs/jay/groups/20/panwei/wan01299/collider_bias/BMI/MV_Egger/new/new_gwas/110K_QCed0.1.bed",simple_names = TRUE)
Indsnps=read.table("/panfs/jay/groups/20/panwei/wan01299/collider_bias/BMI/MV_Egger/new/new_gwas/110K_QCed0.001.bim",header = FALSE)
snps=read.table("/panfs/jay/groups/20/panwei/wan01299/collider_bias/BMI/MV_Egger/new/new_gwas/110K_QCed0.1.bim",header = FALSE)
pos=read.table("/panfs/jay/groups/20/panwei/wan01299/collider_bias/BMI/MV_Egger/new/real_snp_sim/pleiotropy/fourier_ls-all.bed",sep="\t",header=TRUE)


used_snps=snps[snps$V1==1,]

used_pos=pos[pos$chr=="chr1 ",]
block_list=list()
k=0
for(i in 1:nrow(used_pos)){
  start1=used_pos$start[i]
  stop1=used_pos$stop[i]
  a=numeric()
  for(j in 1:nrow(used_snps)){
    if(used_snps$V4[j]>=start1&used_snps$V4[j]<=stop1){
      a=c(a,used_snps$V2[j])
    }
    
  }
  if(length(a)>0){
    k=k+1
    block_list[[k]]=a
  }
}

a1=c(block_list[[1]],block_list[[2]],block_list[[3]],block_list[[4]], block_list[[5]],block_list[[6]],block_list[[7]],block_list[[8]],block_list[[9]], block_list[[10]])
a2=c(block_list[[length(block_list)]],block_list[[length(block_list)-1]],block_list[[length(block_list)-2]],block_list[[length(block_list)-3]],block_list[[length(block_list)-4]],block_list[[length(block_list)-5]],block_list[[length(block_list)-6]],block_list[[length(block_list)-7]],block_list[[length(block_list)-8]],block_list[[length(block_list)-9]])

non_null_snps=sample(a1,size = 0.15*Nsnps-nIV,replace = FALSE)
IVID=sample( intersect(Indsnps$V2, a1[!a1%in%non_null_snps])  ,size=nIV,replace = FALSE)

null_snps=sample(a2,size = 0.85*Nsnps,replace = FALSE)

Null_z=sample(intersect(Indsnps$V2[Indsnps$V1!=1  ],snps$V2),size=Nsnpsz,replace = FALSE)
 

vp=1-invalid
X_only=c(non_null_snps[1:(Nsnps*0.05-nIV*vp)],IVID[1:(nIV*vp)])
n1=(Nsnps*0.05-nIV*vp)+1
n2=(Nsnps*0.05-nIV*vp)+Nsnps*0.05
Y_only=non_null_snps[n1:n2]
n3=n2+1
n4=n2+Nsnps*0.05-nIV*(1-vp)
 
X_and_Y=c(non_null_snps[n3:n4])
 

simID=c(X_only,Y_only,X_and_Y,null_snps,Null_z )
simIV=IVID
Gmatrix=Gmatrixtotal[,colnames(Gmatrixtotal)%in%simID]
for(i in 1:ncol(Gmatrix)){
  Gmatrix[,i][is.na( Gmatrix[,i])]=mean(Gmatrix[,i],na.rm=TRUE)
}

Geno=list(Gmatrix,simID,simIV)
save(Geno,file = "Geno0.RDATA")

