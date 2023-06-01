setwd("C:/Users/wangp/Downloads/BMI")
data=read.table(file="Meta-analysis_Locke_et_al+UKBiobank_2018_top_941_from_COJO_analysis_UPDATED.txt",sep="\t",header = TRUE)
data=data[data$P<=5e-8,]
data=na.omit(data)
BMI_validation_large=data.frame(BP=data$POS,chr=data$CHR)
write.table(BMI_validation_large,file="BMI_validation_large.txt",sep="\t",row.names = FALSE,col.names = TRUE,quote=FALSE)