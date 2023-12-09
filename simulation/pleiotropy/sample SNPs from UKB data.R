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



 
