
install.packages("VennDiagram")               

library("VennDiagram")     


setwd("D:/make vein plots/BMI")
if(1>0){
  filename1=numeric()
  filename2=numeric()
  for(i in 1:10){
    filename1[i]=paste("Res_ldblockM1_second_turn",i,".RData",sep="")
    filename2[i]=paste("BMI_M1_",i,".jpeg",sep="")
  }
  for(i in 1:10){
    load(filename1[i])
    a1=Res$a1
    a2=Res$a2
    a3=Res$a3
    a4=Res$a4
    jpeg(filename = filename2[i],width=500,height=500,res=100)
    
    draw.quad.venn(area1=length(a1), area2=length(a2), area3=length(a3), 
                   area4 =length(a4), n12=sum(a1%in%a2), n23=sum(a2%in%a3), n13=sum(a1%in%a3), 
                   n14= sum(a1%in%a4),n24=sum(a2%in%a4), n34=sum(a3%in%a4), n123=length(intersect(intersect(a1,a2),a3)), 
                   n124=length(intersect(intersect(a1,a2),a4)), n234=length(intersect(intersect(a2,a3),a4)), n134=length(intersect(intersect(a1,a3),a4)), n1234=length(intersect(intersect(a1,a2),intersect(a3,a4))), 
                   category=c("before","after","UKB validation","other validation"),
                   col="Green",fill=c("Yellow","Pink","Blue","Orange"),lty="dashed")
    dev.off()
  }
  
  filename3=numeric()
  filename4=numeric()
  for(i in 1:10){
    filename3[i]=paste("Res_ldblockM2_second_turn",i,".RData",sep="")
    filename4[i]=paste("BMI_M2_",i,".jpeg",sep="")
  }
  for(i in 1:10){
    load(filename3[i])
    a1=Res$a1
    a2=Res$a2
    a3=Res$a3
    a4=Res$a4
    jpeg(filename = filename4[i],width=500,height=500,res=100)
    
    draw.quad.venn(area1=length(a1), area2=length(a2), area3=length(a3), 
                   area4 =length(a4), n12=sum(a1%in%a2), n23=sum(a2%in%a3), n13=sum(a1%in%a3), 
                   n14= sum(a1%in%a4),n24=sum(a2%in%a4), n34=sum(a3%in%a4), n123=length(intersect(intersect(a1,a2),a3)), 
                   n124=length(intersect(intersect(a1,a2),a4)), n234=length(intersect(intersect(a2,a3),a4)), n134=length(intersect(intersect(a1,a3),a4)), n1234=length(intersect(intersect(a1,a2),intersect(a3,a4))), 
                   category=c("M2","cML-M1","UKB validation","other validation"),
                   col="Green",fill=c("Yellow","Pink","Blue","Orange"),lty="dashed")
    dev.off()
  }
  
  
}


setwd("D:/make vein plots/DBP")
if(1>0){
  filename1=numeric()
  filename2=numeric()
  for(i in 1:10){
    filename1[i]=paste("Res_ldblockM1_second_turn",i,".RData",sep="")
    filename2[i]=paste("DBP_M1_",i,".jpeg",sep="")
  }
  for(i in 1:10){
    load(filename1[i])
    a1=Res$a1
    a2=Res$a2
    a3=Res$a3
    a4=Res$a4
    jpeg(filename = filename2[i],width=500,height=500,res=100)
    
    draw.quad.venn(area1=length(a1), area2=length(a2), area3=length(a3), 
                   area4 =length(a4), n12=sum(a1%in%a2), n23=sum(a2%in%a3), n13=sum(a1%in%a3), 
                   n14= sum(a1%in%a4),n24=sum(a2%in%a4), n34=sum(a3%in%a4), n123=length(intersect(intersect(a1,a2),a3)), 
                   n124=length(intersect(intersect(a1,a2),a4)), n234=length(intersect(intersect(a2,a3),a4)), n134=length(intersect(intersect(a1,a3),a4)), n1234=length(intersect(intersect(a1,a2),intersect(a3,a4))), 
                   category=c("before","after","UKB validation","other validation"),
                   col="Green",fill=c("Yellow","Pink","Blue","Orange"),lty="dashed")
    dev.off()
  }
  
  filename3=numeric()
  filename4=numeric()
  for(i in 1:10){
    filename3[i]=paste("Res_ldblockM2_second_turn",i,".RData",sep="")
    filename4[i]=paste("DBP_M2_",i,".jpeg",sep="")
  }
  for(i in 1:10){
    load(filename3[i])
    a1=Res$a1
    a2=Res$a2
    a3=Res$a3
    a4=Res$a4
    jpeg(filename = filename4[i],width=500,height=500,res=100)
    
    draw.quad.venn(area1=length(a1), area2=length(a2), area3=length(a3), 
                   area4 =length(a4), n12=sum(a1%in%a2), n23=sum(a2%in%a3), n13=sum(a1%in%a3), 
                   n14= sum(a1%in%a4),n24=sum(a2%in%a4), n34=sum(a3%in%a4), n123=length(intersect(intersect(a1,a2),a3)), 
                   n124=length(intersect(intersect(a1,a2),a4)), n234=length(intersect(intersect(a2,a3),a4)), n134=length(intersect(intersect(a1,a3),a4)), n1234=length(intersect(intersect(a1,a2),intersect(a3,a4))), 
                   category=c("M2","cML-M1","UKB validation","other validation"),
                   col="Green",fill=c("Yellow","Pink","Blue","Orange"),lty="dashed")
    dev.off()
  }
  
  
}


setwd("D:/make vein plots/SBP")
if(1>0){
  filename1=numeric()
  filename2=numeric()
  for(i in 1:10){
    filename1[i]=paste("Res_ldblockM1_second_turn",i,".RData",sep="")
    filename2[i]=paste("SBP_M1_",i,".jpeg",sep="")
  }
  for(i in 1:10){
    load(filename1[i])
    a1=Res$a1
    a2=Res$a2
    a3=Res$a3
    a4=Res$a4
    jpeg(filename = filename2[i],width=500,height=500,res=100)
    
    draw.quad.venn(area1=length(a1), area2=length(a2), area3=length(a3), 
                   area4 =length(a4), n12=sum(a1%in%a2), n23=sum(a2%in%a3), n13=sum(a1%in%a3), 
                   n14= sum(a1%in%a4),n24=sum(a2%in%a4), n34=sum(a3%in%a4), n123=length(intersect(intersect(a1,a2),a3)), 
                   n124=length(intersect(intersect(a1,a2),a4)), n234=length(intersect(intersect(a2,a3),a4)), n134=length(intersect(intersect(a1,a3),a4)), n1234=length(intersect(intersect(a1,a2),intersect(a3,a4))), 
                   category=c("before","after","UKB validation","other validation"),
                   col="Green",fill=c("Yellow","Pink","Blue","Orange"),lty="dashed")
    dev.off()
  }
  
  filename3=numeric()
  filename4=numeric()
  for(i in 1:10){
    filename3[i]=paste("Res_ldblockM2_second_turn",i,".RData",sep="")
    filename4[i]=paste("SBP_M2_",i,".jpeg",sep="")
  }
  for(i in 1:10){
    load(filename3[i])
    a1=Res$a1
    a2=Res$a2
    a3=Res$a3
    a4=Res$a4
    jpeg(filename = filename4[i],width=500,height=500,res=100)
    
    draw.quad.venn(area1=length(a1), area2=length(a2), area3=length(a3), 
                   area4 =length(a4), n12=sum(a1%in%a2), n23=sum(a2%in%a3), n13=sum(a1%in%a3), 
                   n14= sum(a1%in%a4),n24=sum(a2%in%a4), n34=sum(a3%in%a4), n123=length(intersect(intersect(a1,a2),a3)), 
                   n124=length(intersect(intersect(a1,a2),a4)), n234=length(intersect(intersect(a2,a3),a4)), n134=length(intersect(intersect(a1,a3),a4)), n1234=length(intersect(intersect(a1,a2),intersect(a3,a4))), 
                   category=c("M2","cML-M1","UKB validation","other validation"),
                   col="Green",fill=c("Yellow","Pink","Blue","Orange"),lty="dashed")
    dev.off()
  }
  
  
}

setwd("D:/make vein plots/Height")
if(1>0){
  filename1=numeric()
  filename2=numeric()
  for(i in 1:10){
    filename1[i]=paste("Res_ldblockM1_second_turn",i,".RData",sep="")
    filename2[i]=paste("Height_M1_",i,".jpeg",sep="")
  }
  for(i in 1:10){
    load(filename1[i])
    a1=Res$a1
    a2=Res$a2
    a3=Res$a3
    a4=Res$a4
    jpeg(filename = filename2[i],width=500,height=500,res=100)
    
    draw.quad.venn(area1=length(a1), area2=length(a2), area3=length(a3), 
                   area4 =length(a4), n12=sum(a1%in%a2), n23=sum(a2%in%a3), n13=sum(a1%in%a3), 
                   n14= sum(a1%in%a4),n24=sum(a2%in%a4), n34=sum(a3%in%a4), n123=length(intersect(intersect(a1,a2),a3)), 
                   n124=length(intersect(intersect(a1,a2),a4)), n234=length(intersect(intersect(a2,a3),a4)), n134=length(intersect(intersect(a1,a3),a4)), n1234=length(intersect(intersect(a1,a2),intersect(a3,a4))), 
                   category=c("before","after","UKB validation","other validation"),
                   col="Green",fill=c("Yellow","Pink","Blue","Orange"),lty="dashed")
    dev.off()
  }
  
  filename3=numeric()
  filename4=numeric()
  for(i in 1:10){
    filename3[i]=paste("Res_ldblockM2_second_turn",i,".RData",sep="")
    filename4[i]=paste("Height_M2_",i,".jpeg",sep="")
  }
  for(i in 1:10){
    load(filename3[i])
    a1=Res$a1
    a2=Res$a2
    a3=Res$a3
    a4=Res$a4
    jpeg(filename = filename4[i],width=500,height=500,res=100)
    
    draw.quad.venn(area1=length(a1), area2=length(a2), area3=length(a3), 
                   area4 =length(a4), n12=sum(a1%in%a2), n23=sum(a2%in%a3), n13=sum(a1%in%a3), 
                   n14= sum(a1%in%a4),n24=sum(a2%in%a4), n34=sum(a3%in%a4), n123=length(intersect(intersect(a1,a2),a3)), 
                   n124=length(intersect(intersect(a1,a2),a4)), n234=length(intersect(intersect(a2,a3),a4)), n134=length(intersect(intersect(a1,a3),a4)), n1234=length(intersect(intersect(a1,a2),intersect(a3,a4))), 
                   category=c("M2","cML-M1","UKB validation","other validation"),
                   col="Green",fill=c("Yellow","Pink","Blue","Orange"),lty="dashed")
    dev.off()
  }
  
  
}

setwd("D:/make vein plots/WHR")
if(1>0){
  filename1=numeric()
  filename2=numeric()
  for(i in 1:10){
    filename1[i]=paste("Res_ldblockM1_second_turn",i,".RData",sep="")
    filename2[i]=paste("WHR_M1_",i,".jpeg",sep="")
  }
  for(i in 1:10){
    load(filename1[i])
    a1=Res$a1
    a2=Res$a2
    a3=Res$a3
    a4=Res$a4
    jpeg(filename = filename2[i],width=500,height=500,res=100)
    
    draw.quad.venn(area1=length(a1), area2=length(a2), area3=length(a3), 
                   area4 =length(a4), n12=sum(a1%in%a2), n23=sum(a2%in%a3), n13=sum(a1%in%a3), 
                   n14= sum(a1%in%a4),n24=sum(a2%in%a4), n34=sum(a3%in%a4), n123=length(intersect(intersect(a1,a2),a3)), 
                   n124=length(intersect(intersect(a1,a2),a4)), n234=length(intersect(intersect(a2,a3),a4)), n134=length(intersect(intersect(a1,a3),a4)), n1234=length(intersect(intersect(a1,a2),intersect(a3,a4))), 
                   category=c("before","after","UKB validation","other validation"),
                   col="Green",fill=c("Yellow","Pink","Blue","Orange"),lty="dashed")
    dev.off()
  }
  
  filename3=numeric()
  filename4=numeric()
  for(i in 1:10){
    filename3[i]=paste("Res_ldblockM2_second_turn",i,".RData",sep="")
    filename4[i]=paste("WHR_M2_",i,".jpeg",sep="")
  }
  for(i in 1:10){
    load(filename3[i])
    a1=Res$a1
    a2=Res$a2
    a3=Res$a3
    a4=Res$a4
    jpeg(filename = filename4[i],width=500,height=500,res=100)
    
    draw.quad.venn(area1=length(a1), area2=length(a2), area3=length(a3), 
                   area4 =length(a4), n12=sum(a1%in%a2), n23=sum(a2%in%a3), n13=sum(a1%in%a3), 
                   n14= sum(a1%in%a4),n24=sum(a2%in%a4), n34=sum(a3%in%a4), n123=length(intersect(intersect(a1,a2),a3)), 
                   n124=length(intersect(intersect(a1,a2),a4)), n234=length(intersect(intersect(a2,a3),a4)), n134=length(intersect(intersect(a1,a3),a4)), n1234=length(intersect(intersect(a1,a2),intersect(a3,a4))), 
                   category=c("M2","cML-M1","UKB validation","other validation"),
                   col="Green",fill=c("Yellow","Pink","Blue","Orange"),lty="dashed")
    dev.off()
  }
  
  
}


setwd("D:/make vein plots/WC")
if(1>0){
  filename1=numeric()
  filename2=numeric()
  for(i in 1:10){
    filename1[i]=paste("Res_ldblockM1_second_turn",i,".RData",sep="")
    filename2[i]=paste("WC_M1_",i,".jpeg",sep="")
  }
  for(i in 1:10){
    load(filename1[i])
    a1=Res$a1
    a2=Res$a2
    a3=Res$a3
    a4=Res$a4
    jpeg(filename = filename2[i],width=500,height=500,res=100)
    
    draw.quad.venn(area1=length(a1), area2=length(a2), area3=length(a3), 
                   area4 =length(a4), n12=sum(a1%in%a2), n23=sum(a2%in%a3), n13=sum(a1%in%a3), 
                   n14= sum(a1%in%a4),n24=sum(a2%in%a4), n34=sum(a3%in%a4), n123=length(intersect(intersect(a1,a2),a3)), 
                   n124=length(intersect(intersect(a1,a2),a4)), n234=length(intersect(intersect(a2,a3),a4)), n134=length(intersect(intersect(a1,a3),a4)), n1234=length(intersect(intersect(a1,a2),intersect(a3,a4))), 
                   category=c("before","after","UKB validation","other validation"),
                   col="Green",fill=c("Yellow","Pink","Blue","Orange"),lty="dashed")
    dev.off()
  }
  
  filename3=numeric()
  filename4=numeric()
  for(i in 1:10){
    filename3[i]=paste("Res_ldblockM2_second_turn",i,".RData",sep="")
    filename4[i]=paste("WC_M2_",i,".jpeg",sep="")
  }
  for(i in 1:10){
    load(filename3[i])
    a1=Res$a1
    a2=Res$a2
    a3=Res$a3
    a4=Res$a4
    jpeg(filename = filename4[i],width=500,height=500,res=100)
    
    draw.quad.venn(area1=length(a1), area2=length(a2), area3=length(a3), 
                   area4 =length(a4), n12=sum(a1%in%a2), n23=sum(a2%in%a3), n13=sum(a1%in%a3), 
                   n14= sum(a1%in%a4),n24=sum(a2%in%a4), n34=sum(a3%in%a4), n123=length(intersect(intersect(a1,a2),a3)), 
                   n124=length(intersect(intersect(a1,a2),a4)), n234=length(intersect(intersect(a2,a3),a4)), n134=length(intersect(intersect(a1,a3),a4)), n1234=length(intersect(intersect(a1,a2),intersect(a3,a4))), 
                   category=c("M2","cML-M1","UKB validation","other validation"),
                   col="Green",fill=c("Yellow","Pink","Blue","Orange"),lty="dashed")
    dev.off()
  }
  
  
}

setwd("D:/make vein plots/HC")
if(1>0){
  filename1=numeric()
  filename2=numeric()
  for(i in 1:10){
    filename1[i]=paste("Res_ldblockM1_second_turn",i,".RData",sep="")
    filename2[i]=paste("HC_M1_",i,".jpeg",sep="")
  }
  for(i in 1:10){
    load(filename1[i])
    a1=Res$a1
    a2=Res$a2
    a3=Res$a3
    a4=Res$a4
    jpeg(filename = filename2[i],width=500,height=500,res=100)
    
    draw.quad.venn(area1=length(a1), area2=length(a2), area3=length(a3), 
                   area4 =length(a4), n12=sum(a1%in%a2), n23=sum(a2%in%a3), n13=sum(a1%in%a3), 
                   n14= sum(a1%in%a4),n24=sum(a2%in%a4), n34=sum(a3%in%a4), n123=length(intersect(intersect(a1,a2),a3)), 
                   n124=length(intersect(intersect(a1,a2),a4)), n234=length(intersect(intersect(a2,a3),a4)), n134=length(intersect(intersect(a1,a3),a4)), n1234=length(intersect(intersect(a1,a2),intersect(a3,a4))), 
                   category=c("before","after","UKB validation","other validation"),
                   col="Green",fill=c("Yellow","Pink","Blue","Orange"),lty="dashed")
    dev.off()
  }
  
  filename3=numeric()
  filename4=numeric()
  for(i in 1:10){
    filename3[i]=paste("Res_ldblockM2_second_turn",i,".RData",sep="")
    filename4[i]=paste("HC_M2_",i,".jpeg",sep="")
  }
  for(i in 1:10){
    load(filename3[i])
    a1=Res$a1
    a2=Res$a2
    a3=Res$a3
    a4=Res$a4
    jpeg(filename = filename4[i],width=500,height=500,res=100)
    
    draw.quad.venn(area1=length(a1), area2=length(a2), area3=length(a3), 
                   area4 =length(a4), n12=sum(a1%in%a2), n23=sum(a2%in%a3), n13=sum(a1%in%a3), 
                   n14= sum(a1%in%a4),n24=sum(a2%in%a4), n34=sum(a3%in%a4), n123=length(intersect(intersect(a1,a2),a3)), 
                   n124=length(intersect(intersect(a1,a2),a4)), n234=length(intersect(intersect(a2,a3),a4)), n134=length(intersect(intersect(a1,a3),a4)), n1234=length(intersect(intersect(a1,a2),intersect(a3,a4))), 
                   category=c("M2","cML-M1","UKB validation","other validation"),
                   col="Green",fill=c("Yellow","Pink","Blue","Orange"),lty="dashed")
    dev.off()
  }
  
  
}
