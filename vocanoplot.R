### definition of volcanoplot function
#
#firt colunm is analytes, second columnis p value, third column is ratio
#
library(ggplot2)
library(ggrepel)
library(openxlsx)
library(dplyr)
library(magrittr)


volcanoplotfdr<-function(pfc,title="volcanoplot",fccutoff=2, label="labelsig",labelsize=3,pointsize=1, pcutoff=0.05, fdrcutoff=2){
  names(pfc)[1:4]<-c("analytes", "pvalue","ratio", "fdr")
  
  #pfc <- pfc %>% drop_na()
  
  lowestp<-min(pfc$pvalue[pfc$pvalue>0])
  pfc$pvalue[pfc$pvalue==0]<-min(lowestp,1e-30 )
 

  pfc[, "colgroup"] <- "black"
  pfdrcut <- pfc$pvalue  < pcutoff & pfc$fdr < fdrcutoff
  pfdrratio_p <-  pfdrcut & pfc$ratio >=fccutoff
  pfdrratio_n <-  pfdrcut & pfc$ratio <=1/fccutoff
  
  pfc[pfdrcut, "colgroup"] <- "orange"
  pfc[pfdrratio_p, "colgroup"] <- "red" 
  pfc[pfdrratio_n, "colgroup"] <- "green" 
  pfc$colgroup<-as.factor(pfc$colgroup)
  
  changed <- c(total= nrow(pfc),increase=nrow(pfc[pfc$colgroup=="red",]), decrease=nrow(pfc[pfc$colgroup=="green",]))
  
  changed <- paste(changed, collapse=":")
  
  
  p<-ggplot(data=pfc, aes(log(ratio,2), -log(pvalue,10)))+
    theme_bw()+ 
    geom_point(size=pointsize, alpha=ifelse(pfc$colgroup=="black", 1, 1), color=pfc$colgroup)+
    geom_hline(yintercept = -log10(pcutoff),lty=3,lwd=0.5)+
    geom_vline(xintercept =log2(fccutoff) ,lty=3,lwd=0.5)+
    geom_vline(xintercept =log2(1/fccutoff),lty=3,lwd=0.5)+
    labs(x="log2(Ratio)",y="-log10(p value)")+
    ggtitle(paste(title,changed))+
    theme(plot.title = element_text(color="black", face="bold", size=12, hjust=0))+
    theme(axis.title = element_text(color="black", face="bold", size=12))+
    theme(plot.margin=unit(c(2,2,2,2),"line"))
  if (label=="labelsig") {
    p + geom_text_repel(aes(label=substr(analytes,1,25)),
                        segment.size = 0.06,
                        min.segment.length = 0,
                        max.overlaps=100,
                        data = subset(pfc,colgroup=="red"|colgroup=="green"), size=labelsize)
    
  } else if(label=="nolabel")
  { 
    p
  }else if(label=="labelall"){
    
    p+geom_text_repel(aes(label=analytes), size=labelsize)
    
  }
}
