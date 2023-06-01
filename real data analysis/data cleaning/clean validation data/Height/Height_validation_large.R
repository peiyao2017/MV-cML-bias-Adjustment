setwd("C:/Users/wangp/Downloads/Height")
data=read.table("Meta-analysis_Wood_et_al+UKBiobank_2018_top_3290_from_COJO_analysis.txt",header = TRUE,sep="\t")
data=data[data$P<=5e-8,]
data=na.omit(data)
Height_validation_large=data.frame(BP=data$POS,chr=data$CHR)
write.table(Height_validation_large,file="Height_validation_large.txt",sep="\t",row.names = FALSE,col.names = TRUE,quote=FALSE)