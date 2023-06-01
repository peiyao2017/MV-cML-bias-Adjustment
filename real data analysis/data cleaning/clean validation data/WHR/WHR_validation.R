library(data.table)
setwd("C:/Users/wangp/Downloads/WHR")

data=fread("whr.giant-ukbb.meta-analysis.combined.23May2018.txt",select = c("CHR","POS","P"),sep=" ")
head(data)
data=data[data$P<=5e-8,]
WHR_validation_large=data.frame(BP=data$POS,chr=data$CHR)
write.table(WHR_validation_large,file="WHR_validation",sep="\t",row.names = FALSE,col.names = TRUE,quote=FALSE)