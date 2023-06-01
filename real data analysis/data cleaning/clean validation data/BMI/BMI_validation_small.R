library(R.utils)
setwd("C:/Users/wangp/Downloads/BMI")
data=gunzip("C:/Users/wangp/Downloads/BMI/raw UKB validation/21001_irnt.gwas.imputed_v3.both_sexes.tsv.bgz", "BMI.tsv")
data=fread("BMI.tsv",sep="\t")
data=data[data$pval<=5e-8,]
vari=data$variant
chr=numeric()
BP=numeric()
for(i in 1:length(vari)){
  chr[i]=substr(vari[i],start=1,stop=unlist(gregexpr(':',vari[i] ))[1]-1)
  BP[i]=substr(vari[i],start=unlist(gregexpr(':',vari[i] ))[1]+1,stop=unlist(gregexpr(':',vari[i] ))[2]-1)
  print(i)
}
chr=chr[chr!="X"]
BP=BP[1:length(chr)]
BP=na.omit(BP)
chr=na.omit(chr)
BP=as.numeric(BP)
chr=as.numeric(chr)
BMI_validation_small=data.frame(BP=BP,chr=chr)
write.table(BMI_validation_small,file="BMI_validation_small.txt",sep="\t",row.names = FALSE,col.names = TRUE,quote=FALSE)