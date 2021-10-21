#
#Olink Human Plasma data analysis
#
#library(outliers)
#library(beeswarm)
library(ggplot2)
library(magrittr)
library(mixOmics)
library(ggrepel)
# library(survminer)
# library(survival)
library(purrr)
library(tidyr)
library(dplyr)
library(openxlsx)
library(tibble)
library(stringr)
library(ggpubr)
rm(list=ls())
# removeprotein <- read.csv("C:/zhijuncao/cardiotoxicity/SomaLogic Returned Data/11-29-2016/soma 5 protein and calibrators remove.csv", sep=",", header=TRUE)
# somaid <- read.csv("C:/zhijuncao/cardiotoxicity/SomaLogic Returned Data/11-29-2016/somaid1.csv", sep=",", header=TRUE)


ratio <- function(x){x/x[1]}
###somascan
soma <- read.csv("C:/zhijuncao/cardiotoxicity/SomaLogic Returned Data/04142018/sample_replace_duplicate_with_mean_04142018_subject33labchanged_time.csv")
somaid <- read.csv("C:/zhijuncao/cardiotoxicity/SomaLogic Returned Data/11-29-2016/somaid1.csv", sep=",", header=TRUE)


soma_ratio <- soma %>% filter(SubjectID!=15,SubjectID!=44,SubjectID!=45) %>% group_by(SubjectID) %>% mutate_at(8:1312, ratio) %>% data.frame()
write.xlsx(soma_ratio, "C:/zhijuncao/cardiotoxicity/DOX_Olink/12162019/sample_replace_duplicate_with_mean_04142018_subject33labchanged_time_ratio.xlsx")

##### olink bridge 2018 and 2019 data
olinkid <- read.xlsx("C:/zhijuncao/cardiotoxicity/DOX_Olink/20190281/20190281_Yu_NPX.xlsx", sheet = "OlinkID")
toIntensity <- function(x){2**x} # olink data is log2 scale, transformed back to NXP unit

### olink 2018 set
olink2018 <- read.xlsx("C:/zhijuncao/cardiotoxicity/DOX_Olink/20180576_Yu_NPX_LOD_for Statistic analysis.xlsx", sheet ="NPX Data 2018")
olink2018_intensity <- olink2018 %>% mutate_if(is.numeric, toIntensity)
olink2018meta <- read.xlsx("C:/zhijuncao/cardiotoxicity/DOX_Olink/20180576_Yu_NPX_LOD_for Statistic analysis.xlsx", sheet ="metasample")
olink2018meta <- olink2018meta %>% mutate_if(is.numeric, as.character)
olink2018_meta_intensity <- inner_join(olink2018meta, olink2018_intensity)
olink2018_intensity_ratio <- olink2018_meta_intensity %>% group_by(SubjectID) %>% mutate_if(is.numeric, ratio) %>% data.frame()

write.xlsx(list(olink2018_meta_intensity=olink2018_meta_intensity,olink2018_intensity_ratio=olink2018_intensity_ratio), 
                 "C:/zhijuncao/cardiotoxicity/DOX_Olink/12162019/olink2018.xlsx")

### olink 2019 set
olink2019 <- read.xlsx("C:/zhijuncao/cardiotoxicity/DOX_Olink/20190281/20190281_Yu_NPX.xlsx", sheet = "NPX Data R")
olink2019_intensity <- olink2019 %>% mutate_if(is.numeric, toIntensity)
olink2019meta <- read.xlsx("C:/zhijuncao/cardiotoxicity/DOX_Olink/20190714_Yu_sample_metadata.xlsx", sheet = "metasample")
olink2019meta <- olink2019meta %>% mutate_if(is.numeric, as.character) %>% filter(SubjectID!="44")
olink2019_meta_intensity <- inner_join(olink2019meta, olink2019_intensity)
olink2019_intensity_ratio <- olink2019_meta_intensity %>% group_by(SubjectID) %>% mutate_if(is.numeric, ratio) %>% data.frame()
write.xlsx(list(olink2019_meta_intensity=olink2019_meta_intensity,olink2019_intensity_ratio=olink2019_intensity_ratio), 
           "C:/zhijuncao/cardiotoxicity/DOX_Olink/12162019/olink2019.xlsx")

olink2019_meta_intensity1 <- olink2019_meta_intensity %>% 
  filter(!(SubjectID %in% c("84", "85","86", "87", "88", "89", "90")))%>% 
  select(-c(OID00563:OID00654))
olink2019_intensity_ratio1 <- olink2019_intensity_ratio %>% 
                filter(!(SubjectID %in% c("84", "85","86", "87", "88", "89", "90"))) %>% 
                select(-c(OID00563:OID00654))
write.xlsx(list(olink2019_meta_intensity=olink2019_meta_intensity1,olink2019_intensity_ratio=olink2019_intensity_ratio1), 
           "C:/zhijuncao/cardiotoxicity/DOX_Olink/12162019/olink2019_to2panelsSubject83.xlsx")

#####################################################################################################################
##########              function:     wilcox and Welch's t test
#####################################################################################################################
ttest <- function(sample_mean1=allsample,Firsti=9, Lasti=284){
  
  resultfoldp<-data.frame("analytes"=character(),
                          "A.s.t.p"=numeric(),"A.s.w.p"=numeric(),"A.s.fc"=numeric(),"A.s.fcmd"=numeric(),
                          "B.s.t.p"=numeric(),"B.s.w.p"=numeric(),"B.s.fc"=numeric(),"B.s.fcmd"=numeric(),
                          "C.s.t.p"=numeric(),"C.s.w.p"=numeric(),"C.s.fc"=numeric(),"C.s.fcmd"=numeric(),
                          "N.b.t.p"=numeric(),"N.b.w.p"=numeric(),"N.b.fc"=numeric(),"N.b.fcmd"=numeric(),
                          "N.c.t.p"=numeric(),"N.c.w.p"=numeric(),"N.c.fc"=numeric(),"N.c.fcmd"=numeric(),
                          "S.b.t.p"=numeric(),"S.b.w.p"=numeric(),"S.b.fc"=numeric(),"S.b.fcmd"=numeric(),
                          "S.c.t.p"=numeric(),"S.c.w.p"=numeric(),"S.c.fc"=numeric(),"S.c.fcmd"=numeric(),
                          stringsAsFactors=FALSE)
  for (i in Firsti:Lasti){
    #i=6
    A.n<-sample_mean1[sample_mean1$Status=="Normal"&sample_mean1$TimePoint=="A",i]
    
    A.s<-sample_mean1[sample_mean1$Status=="Abnormal"&sample_mean1$TimePoint=="A",i]
    
    
    A.s.fc<-mean(A.s, na.rm = TRUE)/mean(A.n, na.rm = TRUE)
    A.s.fcmd<-median(A.s, na.rm = TRUE)/median(A.n, na.rm = TRUE)
    
    A.s.t<-t.test(log(A.s,2), log(A.n,2),paired=FALSE,var.qual=FALSE)$p.value
    A.s.w<-wilcox.test(log(A.s,2), log(A.n,2),paired=FALSE,var.qual=FALSE)$p.value
    
    B.n<-sample_mean1[sample_mean1$Status=="Normal"&sample_mean1$TimePoint=="B",i]
    
    B.s<-sample_mean1[sample_mean1$Status=="Abnormal"&sample_mean1$TimePoint=="B",i]
    
    B.s.fc<-mean(B.s, na.rm = TRUE)/mean(B.n, na.rm = TRUE)
    B.s.fcmd<-median(B.s, na.rm = TRUE)/median(B.n, na.rm = TRUE)
    
    if (length(B.s)>=2){
      B.s.t<-t.test(log(B.s,2), log(B.n,2),paired=FALSE,var.qual=FALSE)$p.value
      B.s.w<-wilcox.test(log(B.s,2), log(B.n,2),paired=FALSE,var.qual=FALSE)$p.value
    }    else {B.s.t=NA
    B.s.w=NA}
    
    C.n<-sample_mean1[sample_mean1$Status=="Normal"&sample_mean1$TimePoint=="C",i]
    
    C.s<-sample_mean1[sample_mean1$Status=="Abnormal"&sample_mean1$TimePoint=="C",i]
    
    
    C.s.fc<-mean(C.s, na.rm = TRUE)/mean(C.n, na.rm = TRUE)
    C.s.fcmd<-median(C.s, na.rm = TRUE)/median(C.n, na.rm = TRUE)
    
    C.s.t<-t.test(log(C.s,2), log(C.n,2),paired=FALSE,var.qual=FALSE)$p.value
    C.s.w<-wilcox.test(log(C.s,2), log(C.n,2),paired=FALSE,var.qual=FALSE)$p.value
    
    
    
    N.b.fc<-mean(B.n, na.rm = TRUE)/mean(A.n, na.rm = TRUE)
    N.b.fcmd<-median(B.n, na.rm = TRUE)/median(A.n, na.rm = TRUE)
    
    N.c.fc<-mean(C.n, na.rm = TRUE)/mean(A.n, na.rm = TRUE)
    N.c.fcmd<-median(C.n, na.rm = TRUE)/median(A.n, na.rm = TRUE)
    
    N.b.t<-t.test(log(B.n,2), log(A.n,2),paired=FALSE,var.qual=FALSE)$p.value
    N.b.w<-wilcox.test(log(B.n,2), log(A.n,2),paired=FALSE,var.qual=FALSE)$p.value
    N.c.t<-t.test(log(C.n,2), log(A.n,2),paired=FALSE,var.qual=FALSE)$p.value
    N.c.w<-wilcox.test(log(C.n,2), log(A.n,2),paired=FALSE,var.qual=FALSE)$p.value
    
    S.b.fc<-mean(B.s, na.rm = TRUE)/mean(A.s, na.rm = TRUE)
    S.b.fcmd<-median(B.s, na.rm = TRUE)/median(A.s, na.rm = TRUE)
    
    S.c.fc<-mean(C.s, na.rm = TRUE)/mean(A.s, na.rm = TRUE)
    S.c.fcmd<-median(C.s, na.rm = TRUE)/median(A.s, na.rm = TRUE)
    
    if(length(B.s)>=2){
      S.b.t<-t.test(log(B.s,2), log(A.s,2),paired=FALSE,var.qual=FALSE)$p.value
      S.b.w<-wilcox.test(log(B.s,2), log(A.s,2),paired=FALSE,var.qual=FALSE)$p.value
    }      else {S.b.t=NA
    S.b.w=NA}
    S.c.t<-t.test(log(C.s,2), log(A.s,2),paired=FALSE,var.qual=FALSE)$p.value
    S.c.w<-wilcox.test(log(C.s,2), log(A.s,2),paired=FALSE,var.qual=FALSE)$p.value
    
    
    
    resultfoldp[i-Firsti+1,]<-c(names(sample_mean1)[i],
                                A.s.t,A.s.w, A.s.fc, A.s.fcmd,
                                B.s.t,B.s.w, B.s.fc, B.s.fcmd,
                                C.s.t, C.s.w, C.s.fc, C.s.fcmd,
                                N.b.t, N.b.w, N.b.fc, N.b.fcmd,
                                N.c.t,N.c.w, N.c.fc, N.c.fcmd,
                                S.b.t, S.b.w, S.b.fc, S.b.fcmd,
                                S.c.t,S.c.w, S.c.fc, S.c.fcmd)
    
    
  }
  
  resultfoldp[,-1] <- sapply(resultfoldp[,-1], as.numeric)
  
  
  resultfoldp[,"A.s.t.fdr"]<-p.adjust(resultfoldp$A.s.t.p,"BH")
  
  resultfoldp[,"B.s.t.fdr"]<-p.adjust(resultfoldp$B.s.t.p,"BH")
  
  resultfoldp[,"C.s.t.fdr"]<-p.adjust(resultfoldp$C.s.t.p,"BH")
  
  resultfoldp[,"N.b.t.fdr"]<-p.adjust(resultfoldp$N.b.t.p,"BH")
  resultfoldp[,"N.c.t.fdr"]<-p.adjust(resultfoldp$N.c.t.p,"BH")
  
  
  resultfoldp[,"S.b.t.fdr"]<-p.adjust(resultfoldp$S.b.t.p,"BH")
  resultfoldp[,"S.c.t.fdr"]<-p.adjust(resultfoldp$S.c.t.p,"BH")
  
  ### fdr with wilcox
  resultfoldp[,"A.s.w.fdr"]<-p.adjust(resultfoldp$A.s.w.p,"BH")
  
  resultfoldp[,"B.s.w.fdr"]<-p.adjust(resultfoldp$B.s.w.p,"BH")
  
  resultfoldp[,"C.s.w.fdr"]<-p.adjust(resultfoldp$C.s.w.p,"BH")
  
  resultfoldp[,"N.b.w.fdr"]<-p.adjust(resultfoldp$N.b.w.p,"BH")
  resultfoldp[,"N.c.w.fdr"]<-p.adjust(resultfoldp$N.c.w.p,"BH")
  
  
  resultfoldp[,"S.b.w.fdr"]<-p.adjust(resultfoldp$S.b.w.p,"BH")
  resultfoldp[,"S.c.w.fdr"]<-p.adjust(resultfoldp$S.c.w.p,"BH")
  

  resultfoldp1 <- resultfoldp %>% select(c(analytes, 
                                           A.s.w.p, A.s.w.fdr, A.s.fcmd, B.s.w.p, B.s.w.fdr, B.s.fcmd, C.s.w.p, C.s.w.fdr, C.s.fcmd, 
                                           N.b.w.p, N.b.w.fdr, N.b.fcmd, N.c.w.p, N.c.w.fdr, N.c.fcmd, 
                                           S.b.w.p, S.b.w.fdr, S.b.fcmd, S.c.w.p, S.c.w.fdr, S.c.fcmd,
                                           
                                           A.s.t.p, A.s.t.fdr, A.s.fc,  B.s.t.p, B.s.t.fdr,  B.s.fc, C.s.t.p, C.s.t.fdr, C.s.fc,
                                           N.b.t.p, N.b.t.fdr, N.b.fc,  N.c.t.p, N.c.t.fdr, N.c.fc, 
                                           S.b.t.p, S.b.t.fdr, S.b.fc,  S.c.t.p, S.c.t.fdr, S.c.fc) )
  return (resultfoldp1)
  
}
###
ttest_ratio <- function(sample_mean1=allsample,Firsti=9, Lasti=284){
  
  resultfoldp<-data.frame("analytes"=character(),
                          "A.s.t.p"=numeric(),"A.s.w.p"=numeric(),"A.s.fc"=numeric(),"A.s.fcmd"=numeric(),
                          "B.s.t.p"=numeric(),"B.s.w.p"=numeric(),"B.s.fc"=numeric(),"B.s.fcmd"=numeric(),
                          "C.s.t.p"=numeric(),"C.s.w.p"=numeric(),"C.s.fc"=numeric(),"C.s.fcmd"=numeric(),
                          "N.b.t.p"=numeric(),"N.b.w.p"=numeric(),"N.b.fc"=numeric(),"N.b.fcmd"=numeric(),
                          "N.c.t.p"=numeric(),"N.c.w.p"=numeric(),"N.c.fc"=numeric(),"N.c.fcmd"=numeric(),
                          "S.b.t.p"=numeric(),"S.b.w.p"=numeric(),"S.b.fc"=numeric(),"S.b.fcmd"=numeric(),
                          "S.c.t.p"=numeric(),"S.c.w.p"=numeric(),"S.c.fc"=numeric(),"S.c.fcmd"=numeric(),
                          stringsAsFactors=FALSE)
  for (i in Firsti:Lasti){
    #i=6
    A.n<-sample_mean1[sample_mean1$Status=="Normal"&sample_mean1$TimePoint=="A",i]
    
    A.s<-sample_mean1[sample_mean1$Status=="Abnormal"&sample_mean1$TimePoint=="A",i]
    
    
    A.s.fc<-mean(A.s, na.rm = TRUE)/mean(A.n, na.rm = TRUE)
    A.s.fcmd<-median(A.s, na.rm = TRUE)/median(A.n, na.rm = TRUE)
    
    A.s.t<-t.test(log(A.s,2), log(A.n,2),paired=FALSE,var.qual=FALSE)$p.value
    A.s.w<-wilcox.test(log(A.s,2), log(A.n,2),paired=FALSE,var.qual=FALSE)$p.value
    
    B.n<-sample_mean1[sample_mean1$Status=="Normal"&sample_mean1$TimePoint=="B",i]
    
    B.s<-sample_mean1[sample_mean1$Status=="Abnormal"&sample_mean1$TimePoint=="B",i]
    
    B.s.fc<-mean(B.s, na.rm = TRUE)/mean(B.n, na.rm = TRUE)
    B.s.fcmd<-median(B.s, na.rm = TRUE)/median(B.n, na.rm = TRUE)
    
    if (length(B.s)>=2){
      B.s.t<-t.test(log(B.s,2), log(B.n,2),paired=FALSE,var.qual=FALSE)$p.value
      B.s.w<-wilcox.test(log(B.s,2), log(B.n,2),paired=FALSE,var.qual=FALSE)$p.value
    }    else {B.s.t=NA
    B.s.w=NA}
    
    C.n<-sample_mean1[sample_mean1$Status=="Normal"&sample_mean1$TimePoint=="C",i]
    
    C.s<-sample_mean1[sample_mean1$Status=="Abnormal"&sample_mean1$TimePoint=="C",i]
    
    
    C.s.fc<-mean(C.s, na.rm = TRUE)/mean(C.n, na.rm = TRUE)
    C.s.fcmd<-median(C.s, na.rm = TRUE)/median(C.n, na.rm = TRUE)
    
    C.s.t<-t.test(log(C.s,2), log(C.n,2),paired=FALSE,var.qual=FALSE)$p.value
    C.s.w<-wilcox.test(log(C.s,2), log(C.n,2),paired=FALSE,var.qual=FALSE)$p.value
    
    
    
    N.b.fc<-mean(B.n, na.rm = TRUE)/mean(A.n, na.rm = TRUE)
    N.b.fcmd<-median(B.n, na.rm = TRUE)/median(A.n, na.rm = TRUE)
    
    N.c.fc<-mean(C.n, na.rm = TRUE)/mean(A.n, na.rm = TRUE)
    N.c.fcmd<-median(C.n, na.rm = TRUE)/median(A.n, na.rm = TRUE)
    
    N.b.t<-t.test(log(B.n,2), mu=0,var.qual=FALSE)$p.value
    N.b.w<-wilcox.test(log(B.n,2), mu=0,var.qual=FALSE)$p.value
    N.c.t<-t.test(log(C.n,2), mu=0,var.qual=FALSE)$p.value
    N.c.w<-wilcox.test(log(C.n,2), mu=0,var.qual=FALSE)$p.value
    
    S.b.fc<-mean(B.s, na.rm = TRUE)/mean(A.s, na.rm = TRUE)
    S.b.fcmd<-median(B.s, na.rm = TRUE)/median(A.s, na.rm = TRUE)
    
    S.c.fc<-mean(C.s, na.rm = TRUE)/mean(A.s, na.rm = TRUE)
    S.c.fcmd<-median(C.s, na.rm = TRUE)/median(A.s, na.rm = TRUE)
    
    if(length(B.s)>=2){
      S.b.t<-t.test(log(B.s,2), mu=0,var.qual=FALSE)$p.value
      S.b.w<-wilcox.test(log(B.s,2), mu=0,var.qual=FALSE)$p.value
    }      else {S.b.t=NA
    S.b.w=NA}
    S.c.t<-t.test(log(C.s,2), mu=0,var.qual=FALSE)$p.value
    S.c.w<-wilcox.test(log(C.s,2), mu=0,var.qual=FALSE)$p.value
    
    
    
    resultfoldp[i-Firsti+1,]<-c(names(sample_mean1)[i],
                                A.s.t,A.s.w, A.s.fc, A.s.fcmd,
                                B.s.t,B.s.w, B.s.fc, B.s.fcmd,
                                C.s.t, C.s.w, C.s.fc, C.s.fcmd,
                                N.b.t, N.b.w, N.b.fc, N.b.fcmd,
                                N.c.t,N.c.w, N.c.fc, N.c.fcmd,
                                S.b.t, S.b.w, S.b.fc, S.b.fcmd,
                                S.c.t,S.c.w, S.c.fc, S.c.fcmd)
    
    
  }
  
  resultfoldp[,-1] <- sapply(resultfoldp[,-1], as.numeric)
  
  
  resultfoldp[,"A.s.t.fdr"]<-p.adjust(resultfoldp$A.s.t.p,"BH")
  
  resultfoldp[,"B.s.t.fdr"]<-p.adjust(resultfoldp$B.s.t.p,"BH")
  
  resultfoldp[,"C.s.t.fdr"]<-p.adjust(resultfoldp$C.s.t.p,"BH")
  
  resultfoldp[,"N.b.t.fdr"]<-p.adjust(resultfoldp$N.b.t.p,"BH")
  resultfoldp[,"N.c.t.fdr"]<-p.adjust(resultfoldp$N.c.t.p,"BH")
  
  
  resultfoldp[,"S.b.t.fdr"]<-p.adjust(resultfoldp$S.b.t.p,"BH")
  resultfoldp[,"S.c.t.fdr"]<-p.adjust(resultfoldp$S.c.t.p,"BH")
  
  ### fdr with wilcox
  resultfoldp[,"A.s.w.fdr"]<-p.adjust(resultfoldp$A.s.w.p,"BH")
  
  resultfoldp[,"B.s.w.fdr"]<-p.adjust(resultfoldp$B.s.w.p,"BH")
  
  resultfoldp[,"C.s.w.fdr"]<-p.adjust(resultfoldp$C.s.w.p,"BH")
  
  resultfoldp[,"N.b.w.fdr"]<-p.adjust(resultfoldp$N.b.w.p,"BH")
  resultfoldp[,"N.c.w.fdr"]<-p.adjust(resultfoldp$N.c.w.p,"BH")
  
  
  resultfoldp[,"S.b.w.fdr"]<-p.adjust(resultfoldp$S.b.w.p,"BH")
  resultfoldp[,"S.c.w.fdr"]<-p.adjust(resultfoldp$S.c.w.p,"BH")
  
  
  resultfoldp1 <- resultfoldp %>% select(c(analytes, 
                                           A.s.w.p, A.s.w.fdr, A.s.fcmd, B.s.w.p, B.s.w.fdr, B.s.fcmd, C.s.w.p, C.s.w.fdr, C.s.fcmd, 
                                           N.b.w.p, N.b.w.fdr, N.b.fcmd, N.c.w.p, N.c.w.fdr, N.c.fcmd, 
                                           S.b.w.p, S.b.w.fdr, S.b.fcmd, S.c.w.p, S.c.w.fdr, S.c.fcmd,
                                           
                                           A.s.t.p, A.s.t.fdr, A.s.fc,  B.s.t.p, B.s.t.fdr,  B.s.fc, C.s.t.p, C.s.t.fdr, C.s.fc,
                                           N.b.t.p, N.b.t.fdr, N.b.fc,  N.c.t.p, N.c.t.fdr, N.c.fc, 
                                           S.b.t.p, S.b.t.fdr, S.b.fc,  S.c.t.p, S.c.t.fdr, S.c.fc) )
  return (resultfoldp1)
  
}




###################################end function######################################################################


olink2018_meta_intensity_re <- ttest(sample_mean1=olink2018_meta_intensity, Firsti = 8, Lasti = 99)
olink2018_intensity_ratio_re <- ttest_ratio(sample_mean1=olink2018_intensity_ratio, Firsti = 8, Lasti = 99)
olink2018_2019

olink2019_meta_intensity_re <- ttest(sample_mean1=olink2019_meta_intensity1, Firsti = 8, Lasti = 191)
olink2019_intensity_ratio_re <- ttest_ratio(sample_mean1=olink2019_intensity_ratio1, Firsti = 8, Lasti = 191)

olink2018_2019_intensity_re <- rbind(olink2018_meta_intensity_re,olink2019_meta_intensity_re)
olink2018_2019_intensity_ratio_re <- rbind(olink2018_intensity_ratio_re,olink2019_intensity_ratio_re)

olink2018_2019_intensity_re1 <- inner_join(olinkid,olink2018_2019_intensity_re, by=c("OlinkID"="analytes"))
olink2018_2019_intensity_ratio_re1 <- inner_join(olinkid,olink2018_2019_intensity_ratio_re, by=c("OlinkID"="analytes"))

olink2018_2019_intensity_ratio_re1$A.s.w.p <- olink2018_2019_intensity_re1$A.s.w.p
olink2018_2019_intensity_ratio_re1$A.s.w.fdr <- olink2018_2019_intensity_re1$A.s.w.fdr
olink2018_2019_intensity_ratio_re1$A.s.fcmd <- olink2018_2019_intensity_re1$A.s.fcmd

olink2018_2019_intensity_ratio_re1$A.s.t.p <- olink2018_2019_intensity_re1$A.s.t.p
olink2018_2019_intensity_ratio_re1$A.s.t.fdr <- olink2018_2019_intensity_re1$A.s.t.fdr
olink2018_2019_intensity_ratio_re1$A.s.fc <- olink2018_2019_intensity_re1$A.s.fc


soma_re <- ttest(sample_mean1=soma, Firsti = 8, Lasti = 1312)
soma_ratio_re <- ttest_ratio(sample_mean1=soma_ratio, Firsti = 8, Lasti = 1312)

soma_re1 <- inner_join(somaid,soma_re, by=c("SomaId"="analytes"))
soma_ratio_re1 <- inner_join(somaid, soma_ratio_re, by=c("SomaId"="analytes"))

soma_ratio_re1$A.s.w.p <- soma_re1$A.s.w.p
soma_ratio_re1$A.s.w.fdr <- soma_re1$A.s.w.fdr
soma_ratio_re1$A.s.fcmd <- soma_re1$A.s.fcmd

soma_ratio_re1$A.s.t.p <- soma_re1$A.s.t.p
soma_ratio_re1$A.s.t.fdr <- soma_re1$A.s.t.fdr
soma_ratio_re1$A.s.fc <- soma_re1$A.s.fc

openxlsx::write.xlsx(list(olink2018_2019_intensity=olink2018_2019_intensity_re1,
                olink2018_2019_intensity_ratio=olink2018_2019_intensity_ratio_re1,
                soma=soma_re1,
                soma_ratio=soma_ratio_re1),
           "C:/zhijuncao/cardiotoxicity/DOX_Olink/12192019/olinke_soma_result_12192019.xlsx")


dim(olink2018_2019_intensity_re1)
dim(olink2018_2019_intensity_ratio_re1)
dim(soma_re1)
dim(soma_ratio_re1)
####volcanoplot
vocanplot_wilcox <- function(data=data, protein="Assay",title=title, plotname=plotname){
  pdf(plotname, width = 8, height = 8)
  
  p1 <- volcanoplotfdr(data[, c(protein, "A.s.w.p", "A.s.fcmd", "A.s.w.fdr")], title=paste("Intensity","A.s"),fccutoff=1.2, label="labelsig",labelsize=3, pcutoff=0.05, fdrcutoff=2)
  p2 <- volcanoplotfdr(data[, c(protein, "B.s.w.p", "B.s.fcmd", "B.s.w.fdr")], title=paste(title,"B.s"),fccutoff=1.2, label="labelsig",labelsize=3, pcutoff=0.05, fdrcutoff=2)
  p3 <- volcanoplotfdr(data[, c(protein, "C.s.w.p", "C.s.fcmd", "C.s.w.fdr")], title=paste(title,"C.s"),fccutoff=1.2, label="labelsig",labelsize=3, pcutoff=0.05, fdrcutoff=2)
  
  p4 <- volcanoplotfdr(data[, c(protein, "N.b.w.p", "N.b.fcmd", "N.b.w.fdr")], title=paste(title,"N.b"),fccutoff=1.2, label="labelsig",labelsize=3, pcutoff=0.05, fdrcutoff=2)
  p5 <- volcanoplotfdr(data[, c(protein, "N.c.w.p", "N.c.fcmd", "N.c.w.fdr")], title=paste(title,"N.c"),fccutoff=1.2, label="labelsig",labelsize=3, pcutoff=0.05, fdrcutoff=2)
  p6 <- volcanoplotfdr(data[, c(protein, "S.b.w.p", "S.b.fcmd", "S.b.w.fdr")], title=paste(title,"S.b"),fccutoff=1.2, label="labelsig",labelsize=3, pcutoff=0.05, fdrcutoff=2)
  p7 <- volcanoplotfdr(data[, c(protein, "S.c.w.p", "S.c.fcmd", "S.c.w.fdr")], title=paste(title,"S.c"),fccutoff=1.2, label="labelsig",labelsize=3, pcutoff=0.05, fdrcutoff=2)
  print (p1)
  print (p2)
  print (p3)
  print (p4)
  print (p5)
  print (p6)
  print (p7)

  dev.off()
}


vocanplot_ttest <- function(data=data, protein="Assay",title=title, plotname=plotname){
  pdf(plotname, width = 8, height = 8)

p1 <- volcanoplotfdr(data[, c(protein, "A.s.t.p", "A.s.fc", "A.s.t.fdr")], title=paste("Intensity","A.s"),fccutoff=1.2, label="labelsig",labelsize=3, pcutoff=0.05, fdrcutoff=2)
p2 <-volcanoplotfdr(data[, c(protein, "B.s.t.p", "B.s.fc", "B.s.t.fdr")], title=paste(title,"B.s"),fccutoff=1.2, label="labelsig",labelsize=3, pcutoff=0.05, fdrcutoff=2)
p3 <-volcanoplotfdr(data[, c(protein, "C.s.t.p", "C.s.fc", "C.s.t.fdr")], title=paste(title,"C.s"),fccutoff=1.2, label="labelsig",labelsize=3, pcutoff=0.05, fdrcutoff=2)

p4 <-volcanoplotfdr(data[, c(protein, "N.b.t.p", "N.b.fc", "N.b.t.fdr")], title=paste(title,"N.b"),fccutoff=1.2, label="labelsig",labelsize=3, pcutoff=0.05, fdrcutoff=2)
p5 <-volcanoplotfdr(data[, c(protein, "N.c.t.p", "N.c.fc", "N.c.t.fdr")], title=paste(title,"N.c"),fccutoff=1.2, label="labelsig",labelsize=3, pcutoff=0.05, fdrcutoff=2)
p6 <-volcanoplotfdr(data[, c(protein, "S.b.t.p", "S.b.fc", "S.b.t.fdr")], title=paste(title,"S.b"),fccutoff=1.2, label="labelsig",labelsize=3, pcutoff=0.05, fdrcutoff=2)
p7 <-volcanoplotfdr(data[, c(protein, "S.c.t.p", "S.c.fc", "S.c.t.fdr")], title=paste(title,"S.c"),fccutoff=1.2, label="labelsig",labelsize=3, pcutoff=0.05, fdrcutoff=2)
print (p1)
print (p2)
print (p3)
print (p4)
print (p5)
print (p6)
print (p7)
dev.off()
}

###soma
data <- soma_re1
title <- "Intensity"
protein <- "Target"
plotname <- "C:/zhijuncao/cardiotoxicity/DOX_Olink/12192019/soma_intensity_wilcox_vocanoplot.pdf"
vocanplot_wilcox(data=data, protein="Target",title=title, plotname=plotname)
plotname <- "C:/zhijuncao/cardiotoxicity/DOX_Olink/12192019/soma_intensity_ttest_vocanoplot.pdf"
vocanplot_ttest(data=data, protein="Target",title=title, plotname=plotname)

data <- soma_ratio_re1
title <- "Ratio"
protein <- "Target"
plotname <- "C:/zhijuncao/cardiotoxicity/DOX_Olink/12192019/soma_ratio_wilcox_vocanoplot.pdf"
vocanplot_wilcox(data=data, protein="Target",title=title, plotname=plotname)
plotname <- "C:/zhijuncao/cardiotoxicity/DOX_Olink/12192019/soma_ratio_ttest_vocanoplot.pdf"
vocanplot_ttest(data=data, protein="Target",title=title, plotname=plotname)

###olink
data <- olink2018_2019_intensity_re1
title <- "Intensity"
protein <- "Assay"
plotname <- "C:/zhijuncao/cardiotoxicity/DOX_Olink/12192019/olink2018_2019_intensity_wilcox_vocanoplot.pdf"
vocanplot_wilcox(data=data, protein="Assay",title=title, plotname=plotname)
plotname <- "C:/zhijuncao/cardiotoxicity/DOX_Olink/12192019/olink2018_2019_intensity_ttest_vocanoplot.pdf"
vocanplot_ttest(data=data, protein="Assay",title=title, plotname=plotname)

data <- olink2018_2019_intensity_ratio_re1
title <- "Ratio"
protein <- "Assay"
plotname <- "C:/zhijuncao/cardiotoxicity/DOX_Olink/12192019/olink2018_2019_ratio_wilcox_vocanoplot.pdf"
vocanplot_wilcox(data=data, protein="Assay",title=title, plotname=plotname)
plotname <- "C:/zhijuncao/cardiotoxicity/DOX_Olink/12192019/olink2018_2019_ratio_ttest_vocanoplot.pdf"
vocanplot_ttest(data=data, protein="Assay",title=title, plotname=plotname)



###extract significanlty change proteins' functions
ttest_sig <- function(data=data, p=0.05, fc=1.2){
  all_re_sig <- data %>% filter(
    ((A.s.t.p<0.05) & (A.s.fc>=fc | A.s.fc<=1/fc))|
    ((B.s.t.p<0.05) & (B.s.fc>=fc | B.s.fc<=1/fc))|
    ((C.s.t.p<0.05) & (C.s.fc>=fc | C.s.fc<=1/fc))|
    
    ((N.b.t.p<0.05) & (N.b.fc>=fc | N.b.fc<=1/fc))|
    ((N.c.t.p<0.05) & (N.c.fc>=fc | N.c.fc<=1/fc))|
    
    ((S.b.t.p<0.05) & (S.b.fc>=fc | S.b.fc<=1/fc))|
    ((S.c.t.p<0.05) & (S.c.fc>=fc | S.c.fc<=1/fc))
)
  return (all_re_sig)
}

wilcoxtest_sig <- function(data=data, p=0.05, fc=1.2){
  all_re_sig <- data %>% filter(
      ((A.s.w.p<0.05) & (A.s.fcmd>=fc | A.s.fcmd<=1/fc))|
      ((B.s.w.p<0.05) & (B.s.fcmd>=fc | B.s.fcmd<=1/fc))|
      ((C.s.w.p<0.05) & (C.s.fcmd>=fc | C.s.fcmd<=1/fc))|
      
      ((N.b.w.p<0.05) & (N.b.fcmd>=fc | N.b.fcmd<=1/fc))|
      ((N.c.w.p<0.05) & (N.c.fcmd>=fc | N.c.fcmd<=1/fc))|
      
      ((S.b.w.p<0.05) & (S.b.fcmd>=fc | S.b.fcmd<=1/fc))|
      ((S.c.w.p<0.05) & (S.c.fcmd>=fc | S.c.fcmd<=1/fc))
  )
  return (all_re_sig)
}

olink2018_2019_intensity_wilcoxsig <- wilcoxtest_sig(data=olink2018_2019_intensity_re1, p=0.05, fc=1.2)
olink2018_2019_intensity_ratio_wilcoxsig <- wilcoxtest_sig(data=olink2018_2019_intensity_ratio_re1, p=0.05, fc=1.2)

olink2018_2019_intensity_tsig <- ttest_sig(data=olink2018_2019_intensity_re1, p=0.05, fc=1.2)
olink2018_2019_intensity_ratio_tsig <- ttest_sig(data=olink2018_2019_intensity_ratio_re1, p=0.05, fc=1.2)

soma_wilcoxsig <- wilcoxtest_sig(data=soma_re1, p=0.05, fc=1.2)
soma_ratio_wilcoxsig <- wilcoxtest_sig(data=soma_ratio_re1, p=0.05, fc=1.2)

soma_tsig <- ttest_sig(data=soma_re1, p=0.05, fc=1.2)
soma_ratio_tsig <- ttest_sig(data=soma_ratio_re1, p=0.05, fc=1.2)

openxlsx::write.xlsx(list(olink_intensity_wilcoxsig=olink2018_2019_intensity_wilcoxsig,
                olink_ratio_tsig=olink2018_2019_intensity_ratio_tsig,
                soma_wilcoxsig=soma_wilcoxsig,
                soma_ratio_wilcoxsig=soma_ratio_wilcoxsig,
                
                olink_intensity_tsig=olink2018_2019_intensity_tsig,
                olink_ratio_tsig=olink2018_2019_intensity_ratio_tsig,
                soma_tsig=soma_tsig,
                soma_ratio_tsig=soma_ratio_tsig),
           "C:/zhijuncao/cardiotoxicity/DOX_Olink/12192019/olinke_soma_result_12192019_p05fc1d2.xlsx")


