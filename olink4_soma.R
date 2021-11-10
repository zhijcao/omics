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
library(future)
library(furrr)
library(tidyr)
library(dplyr)
library(openxlsx)
library(tibble)
library(stringr)
library(ggpubr)
rm(list=ls())


somaid <- read.csv("C:/zhijuncao/cardiotoxicity/SomaLogic Returned Data/11-29-2016/somaid1.csv", sep=",", header=TRUE)
sample_ratio <- read.csv("C:/zhijuncao/cardiotoxicity/SomaLogic Returned Data/04142018/06212018/BCsample_ratio_time.csv")
sample <- read.csv("C:/zhijuncao/cardiotoxicity/SomaLogic Returned Data/04142018/sample_replace_duplicate_with_mean_04142018_subject33labchanged_time.csv")
sample <- sample %>% rename(Status=Status1)
sample_ratio <-sample_ratio %>% rename(Status=Status1)



#####################################################################################################################
##########              function:  wilcox and  Welch's t test
#####################################################################################################################
ttest <- function(...) {
  obj<-try(t.test(...), silent=TRUE)
  if (is(obj, "try-error")) return(NA) else return(obj$p.value)
}
wtest <- function(...) {
  obj<-try(wilcox.test(...), silent=TRUE)
  if (is(obj, "try-error")) return(NA) else return(obj$p.value)
}

abnormal_normal <- function(data=sample, timepoint="A", analytes=names(sample)[-c(1:7)]){
  j <- 0
  t.re <- list()
  w.re <- list()
  for (analyte in analytes){
    cat(j, " ")
    j <- j+1
     normal <- data %>% filter(Status=="Normal", TimePoint==timepoint) %>% pluck(analyte)
     abnormal <- data %>% filter(Status=="Abnormal", TimePoint==timepoint)%>% pluck(analyte)
     fc.mean <- mean(abnormal)/mean(normal)
     t.pvalue<-ttest(log(abnormal,2), log(normal,2),paired=FALSE,var.qual=FALSE)
   
     fc.median <- median(abnormal)/median(normal)
      w.pvalue<-wtest(abnormal, normal, paired=FALSE,var.qual=FALSE)
  
      t.re[[analyte]] <- c(fc.mean=fc.mean,t.pvalue=t.pvalue)
      w.re[[analyte]] <- c(fc.median=fc.median, w.pvalue=w.pvalue)
  }
 
  t.re <- t(data.frame(t.re)) %>%data.frame() %>% 
    add_column(t.FDR=p.adjust(.$t.pvalue, method="BH"), .after = 2)
   
  w.re <- t(data.frame(w.re)) %>%data.frame() %>% 
      add_column(w.FDR=p.adjust(.$w.pvalue, method="BH"), .after = 2)

  return (list(t.re=t.re, w.re=w.re))
  
}


toA <- function(data=sample, timepoint="B", status="Normal", analytes=names(sample)[-c(1:7)]){
  j <- 0
  t.re <- list()
  w.re <- list()
  for (analyte in analytes){
     cat(j, " ")
     j <- j+1
     a <- data %>% filter(Status==status, TimePoint=="A") %>% pluck(analyte)
     treat <- data %>% filter(Status==status, TimePoint==timepoint) %>% pluck(analyte)
  
     fc.mean <- mean(treat)/mean(a)
     t.pvalue<-ttest(log(treat,2), log(a,2),paired=FALSE,var.qual=FALSE)
  
     fc.median <- median(treat)/median(a)
     w.pvalue<-wtest(treat, a, paired=FALSE,var.qual=FALSE)
  
     t.re[[analyte]] <- c(fc.mean=fc.mean,t.pvalue=t.pvalue)
     w.re[[analyte]] <- c(fc.median=fc.median, w.pvalue=w.pvalue)
  }
 
  
  t.re <- t(data.frame(t.re)) %>%data.frame() %>% 
    add_column(t.FDR=p.adjust(.$t.pvalue, method="BH"), .after = 2)
  
  w.re <- t(data.frame(w.re)) %>%data.frame() %>% 
    add_column(w.FDR=p.adjust(.$w.pvalue, method="BH"), .after = 2)
  
  return (list(t.re=t.re, w.re=w.re))
}


length(names(sample)[-c(1:7)])
A_abnormal.normal_test <- abnormal_normal(data=sample, timepoint="A", analytes=names(sample)[8:10])
normal_B.A_test <- toA(data=sample, timepoint="B", status="Normal", analytes=names(sample)[8:10])
abnormal_B.A_test <- toA(data=sample, timepoint="B", status="Abnormal", analytes=names(sample)[8:10])

  A_abnormal.normal <- abnormal_normal(data=sample, timepoint="A", analytes=names(sample)[-c(1:7)])
  B_abnormal.normal <- abnormal_normal(data=sample, timepoint="B", analytes=names(sample)[-c(1:7)])
  C_abnormal.normal <- abnormal_normal(data=sample, timepoint="C", analytes=names(sample)[-c(1:7)])
  # D_abnormal.normal <- abnormal_normal(data=sample, timepoint="D", analytes=names(sample)[-c(1:7)])
  # D0_abnormal.normal <- abnormal_normal(data=sample, timepoint="D0", analytes=names(sample)[-c(1:7)])
  # E_abnormal.normal <- abnormal_normal(data=sample, timepoint="E", analytes=names(sample)[-c(1:7)])
  # F_abnormal.normal <- abnormal_normal(data=sample, timepoint="F", analytes=names(sample)[-c(1:7)])

  normal_B.A <- toA(data=sample, timepoint="B", status="Normal", analytes=names(sample)[-c(1:7)])
  normal_C.A <- toA(data=sample, timepoint="C", status="Normal", analytes=names(sample)[-c(1:7)])
  normal_D.A <- toA(data=sample, timepoint="D", status="Normal", analytes=names(sample)[-c(1:7)])
  # normal_D0.A <- toA(data=sample, timepoint="D0", status="Normal", analytes=names(sample)[-c(1:7)])
  # normal_E.A <- toA(data=sample, timepoint="E", status="Normal", analytes=names(sample)[-c(1:7)])
  # normal_F.A <- toA(data=sample, timepoint="F", status="Normal", analytes=names(sample)[-c(1:7)])

  abnormal_B.A <- toA(data=sample, timepoint="B", status="Abnormal", analytes=names(sample)[-c(1:7)])
  abnormal_C.A <- toA(data=sample, timepoint="C", status="Abnormal", analytes=names(sample)[-c(1:7)])
  # abnormal_D.A <- toA(data=sample, timepoint="D", status="Abnormal", analytes=names(sample)[-c(1:7)])
  # abnormal_D0.A <- toA(data=sample, timepoint="D0", status="Abnormal", analytes=names(sample)[-c(1:7)])
  # abnormal_E.A <- toA(data=sample, timepoint="E", status="Abnormal", analytes=names(sample)[-c(1:7)])
  # abnormal_F.A <- toA(data=sample, timepoint="F", status="Abnormal", analytes=names(sample)[-c(1:7)])

  re_all <- function(test="w.re") {
  
    re.all <-list(A_abnormal.normal=A_abnormal.normal[[test]],
              B_abnormal.normal=B_abnormal.normal[[test]],
              C_abnormal.normal=C_abnormal.normal[[test]],
              # D_abnormal.normal=D_abnormal.normal[[test]],
              # D0_abnormal.normal=D0_abnormal.normal[[test]],
              # E_abnormal.normal=E_abnormal.normal[[test]],
              # F_abnormal.normal=F_abnormal.normal[[test]],
              
              normal_B.A=normal_B.A[[test]],
              normal_C.A=normal_C.A[[test]],
              # normal_D.A=normal_D.A[[test]],
              # normal_D0.A=normal_D0.A[[test]],
              # normal_E.A=normal_E.A[[test]],
              # normal_F.A=normal_F.A[[test]],
              # 
              abnormal_B.A=abnormal_B.A[[test]],
              abnormal_C.A=abnormal_C.A[[test]]
              # abnormal_D.A=abnormal_D.A[[test]],
              # abnormal_D0.A=abnormal_D0.A[[test]],
              # abnormal_E.A=abnormal_E.A[[test]],
              # abnormal_F.A=abnormal_F.A[[test]]
              ) 
    return (re.all)
  }
  w.re <-  re_all(test="w.re")
  t.re <-  re_all(test="t.re")
  
  saveRDS( w.re, "C:/zhijuncao/cardiotoxicity/DOX_Olink/somascan_w_test_results.Rds")
  saveRDS(t.re, "C:/zhijuncao/cardiotoxicity/DOX_Olink/somascan_t_test_results.Rds")
  
  w_sig <- function(data, fc=1.2){
     data %>% filter((A_abnormal.normal.w.pvalue<0.05&(A_abnormal.normal.fc.median>=fc|A_abnormal.normal.fc.median<=1/fc))|
                     (B_abnormal.normal.w.pvalue<0.05&(B_abnormal.normal.fc.median>=fc|B_abnormal.normal.fc.median<=1/fc))|
                       (C_abnormal.normal.w.pvalue<0.05&(C_abnormal.normal.fc.median>=fc|C_abnormal.normal.fc.median<=1/fc))|
                       # (D_abnormal.normal.w.pvalue<0.05&(D_abnormal.normal.fc.median>=fc|D_abnormal.normal.fc.median<=1/fc))|
                       # (D0_abnormal.normal.w.pvalue<0.05&(D0_abnormal.normal.fc.median>=fc|D0_abnormal.normal.fc.median<=1/fc))|
                       # (E_abnormal.normal.w.pvalue<0.05&(E_abnormal.normal.fc.median>=fc|E_abnormal.normal.fc.median<=1/fc))|
                       # (F_abnormal.normal.w.pvalue<0.05&(F_abnormal.normal.fc.median>=fc|F_abnormal.normal.fc.median<=1/fc))|
                       
                       (normal_B.A.w.pvalue<0.05&(normal_B.A.fc.median>=fc|normal_B.A.fc.median<=1/fc))|
                       (normal_C.A.w.pvalue<0.05&(normal_C.A.fc.median>=fc|normal_C.A.fc.median<=1/fc))|
                       # (normal_D.A.w.pvalue<0.05&(normal_D.A.fc.median>=fc|normal_D.A.fc.median<=1/fc))|
                       # (normal_D0.A.w.pvalue<0.05&(normal_D0.A.fc.median>=fc|normal_D0.A.fc.median<=1/fc))|
                       # (normal_E.A.w.pvalue<0.05&(normal_E.A.fc.median>=fc|normal_E.A.fc.median<=1/fc))|
                       # (normal_F.A.w.pvalue<0.05&(normal_F.A.fc.median>=fc|normal_F.A.fc.median<=1/fc))|
                       
                       (abnormal_B.A.w.pvalue<0.05&(abnormal_B.A.fc.median>=fc|abnormal_B.A.fc.median<=1/fc))|
                       (abnormal_C.A.w.pvalue<0.05&(abnormal_C.A.fc.median>=fc|abnormal_C.A.fc.median<=1/fc))
                       # (abnormal_D.A.w.pvalue<0.05&(abnormal_D.A.fc.median>=fc|abnormal_D.A.fc.median<=1/fc))|
                       # (abnormal_D0.A.w.pvalue<0.05&(abnormal_D0.A.fc.median>=fc|abnormal_D0.A.fc.median<=1/fc))|
                       # (abnormal_E.A.w.pvalue<0.05&(abnormal_E.A.fc.median>=fc|abnormal_E.A.fc.median<=1/fc))|
                       # (abnormal_F.A.w.pvalue<0.05&(abnormal_F.A.fc.median>=fc|abnormal_F.A.fc.median<=1/fc))
                     )
  }
  
  
  t_sig <- function(data, fc=1.2){
    data %>% filter((A_abnormal.normal.t.pvalue<0.05&(A_abnormal.normal.fc.mean>=fc|A_abnormal.normal.fc.mean<=1/fc))|
                      (B_abnormal.normal.t.pvalue<0.05&(B_abnormal.normal.fc.mean>=fc|B_abnormal.normal.fc.mean<=1/fc))|
                      (C_abnormal.normal.t.pvalue<0.05&(C_abnormal.normal.fc.mean>=fc|C_abnormal.normal.fc.mean<=1/fc))|
                      # (D_abnormal.normal.t.pvalue<0.05&(D_abnormal.normal.fc.mean>=fc|D_abnormal.normal.fc.mean<=1/fc))|
                      # (D0_abnormal.normal.t.pvalue<0.05&(D0_abnormal.normal.fc.mean>=fc|D0_abnormal.normal.fc.mean<=1/fc))|
                      # (E_abnormal.normal.t.pvalue<0.05&(E_abnormal.normal.fc.mean>=fc|E_abnormal.normal.fc.mean<=1/fc))|
                      # (F_abnormal.normal.t.pvalue<0.05&(F_abnormal.normal.fc.mean>=fc|F_abnormal.normal.fc.mean<=1/fc))|
                      
                      (normal_B.A.t.pvalue<0.05&(normal_B.A.fc.mean>=fc|normal_B.A.fc.mean<=1/fc))|
                      (normal_C.A.t.pvalue<0.05&(normal_C.A.fc.mean>=fc|normal_C.A.fc.mean<=1/fc))|
                      # (normal_D.A.t.pvalue<0.05&(normal_D.A.fc.mean>=fc|normal_D.A.fc.mean<=1/fc))|
                      # (normal_D0.A.t.pvalue<0.05&(normal_D0.A.fc.mean>=fc|normal_D0.A.fc.mean<=1/fc))|
                      # (normal_E.A.t.pvalue<0.05&(normal_E.A.fc.mean>=fc|normal_E.A.fc.mean<=1/fc))|
                      # (normal_F.A.t.pvalue<0.05&(normal_F.A.fc.mean>=fc|normal_F.A.fc.mean<=1/fc))|
                      
                      (abnormal_B.A.t.pvalue<0.05&(abnormal_B.A.fc.mean>=fc|abnormal_B.A.fc.mean<=1/fc))|
                      (abnormal_C.A.t.pvalue<0.05&(abnormal_C.A.fc.mean>=fc|abnormal_C.A.fc.mean<=1/fc))
                      # (abnormal_D.A.t.pvalue<0.05&(abnormal_D.A.fc.mean>=fc|abnormal_D.A.fc.mean<=1/fc))|
                      # (abnormal_D0.A.t.pvalue<0.05&(abnormal_D0.A.fc.mean>=fc|abnormal_D0.A.fc.mean<=1/fc))|
                      # (abnormal_E.A.t.pvalue<0.05&(abnormal_E.A.fc.mean>=fc|abnormal_E.A.fc.mean<=1/fc))|
                      # (abnormal_F.A.t.pvalue<0.05&(abnormal_F.A.fc.mean>=fc|abnormal_F.A.fc.mean<=1/fc))
    )
  }
  
  
  w.re_df <- w.re %>% data.frame() %>% rownames_to_column("SomaId") %>% inner_join(somaid, .)
  w.re_df_sig <- w_sig(data=w.re_df, fc=1.2)
  write.xlsx(list(all=w.re_df,sig=w.re_df_sig), "C:/zhijuncao/cardiotoxicity/DOX_Olink/somacan_wilcox_result_all_11102021.xlsx")
  
  w.re_id <-  map(w.re, ~.x%>% rownames_to_column("SomaId") %>% inner_join(somaid,.))
  write.xlsx(w.re_id, "C:/zhijuncao/cardiotoxicity/DOX_Olink/somacan_wilcox_result_separated_11102021.xlsx") 
  
  
  
  t.re_df <- t.re %>% data.frame() %>% rownames_to_column("SomaId") %>% inner_join(somaid, .)
  t.re_df_sig <- t_sig(data=t.re_df, fc=1.2)
  write.xlsx(list(all=t.re_df, sig=t.re_df_sig), "C:/zhijuncao/cardiotoxicity/DOX_Olink/somacan_ttest_result_all_11102021.xlsx") 
  
  t.re_id <-  map(t.re, ~.x%>% rownames_to_column("SomaId") %>% inner_join(somaid,.))
  write.xlsx(t.re_id, "C:/zhijuncao/cardiotoxicity/DOX_Olink/somacan_ttest_result_separated_11102021.xlsx") 

###Ratio
  A_abnormal.normal_ratio <- abnormal_normal(data=sample_ratio, timepoint="A", analytes=names(sample_ratio)[-c(1:7)])
  B_abnormal.normal_ratio <- abnormal_normal(data=sample_ratio, timepoint="B", analytes=names(sample_ratio)[-c(1:7)])
  C_abnormal.normal_ratio <- abnormal_normal(data=sample_ratio, timepoint="C", analytes=names(sample_ratio)[-c(1:7)])
  D_abnormal.normal_ratio <- abnormal_normal(data=sample_ratio, timepoint="D", analytes=names(sample_ratio)[-c(1:7)])
  D0_abnormal.normal_ratio <- abnormal_normal(data=sample_ratio, timepoint="D0", analytes=names(sample_ratio)[-c(1:7)])
  E_abnormal.normal_ratio <- abnormal_normal(data=sample_ratio, timepoint="E", analytes=names(sample_ratio)[-c(1:7)])
  F_abnormal.normal_ratio<- abnormal_normal(data=sample_ratio, timepoint="F", analytes=names(sample_ratio)[-c(1:7)])
  
  normal_B.A_ratio <- toA(data=sample_ratio, timepoint="B", status="Normal", analytes=names(sample_ratio)[-c(1:7)])
  normal_C.A_ratio <- toA(data=sample_ratio, timepoint="C", status="Normal", analytes=names(sample_ratio)[-c(1:7)])
  normal_D.A_ratio<- toA(data=sample_ratio, timepoint="D", status="Normal", analytes=names(sample_ratio)[-c(1:7)])
  normal_D0.A_ratio <- toA(data=sample_ratio, timepoint="D0", status="Normal", analytes=names(sample_ratio)[-c(1:7)])
  normal_E.A_ratio <- toA(data=sample_ratio, timepoint="E", status="Normal", analytes=names(sample_ratio)[-c(1:7)])
  normal_F.A_ratio <- toA(data=sample_ratio, timepoint="F", status="Normal", analytes=names(sample_ratio)[-c(1:7)])
  
  abnormal_B.A_ratio <- toA(data=sample_ratio, timepoint="B", status="Abnormal", analytes=names(sample_ratio)[-c(1:7)])
  abnormal_C.A_ratio <- toA(data=sample_ratio, timepoint="C", status="Abnormal", analytes=names(sample_ratio)[-c(1:7)])
  abnormal_D.A_ratio <- toA(data=sample_ratio, timepoint="D", status="Abnormal", analytes=names(sample_ratio)[-c(1:7)])
  abnormal_D0.A_ratio <- toA(data=sample_ratio, timepoint="D0", status="Abnormal", analytes=names(sample_ratio)[-c(1:7)])
  abnormal_E.A_ratio <- toA(data=sample_ratio, timepoint="E", status="Abnormal", analytes=names(sample_ratio)[-c(1:7)])
  abnormal_F.A_ratio <- toA(data=sample_ratio, timepoint="F", status="Abnormal", analytes=names(sample_ratio)[-c(1:7)])
  
  re_all_ratio <- function(test="w.re") {
    
    re.all <-list(A_abnormal.normal=A_abnormal.normal_ratio[[test]],
                  B_abnormal.normal=B_abnormal.normal_ratio[[test]],
                  C_abnormal.normal=C_abnormal.normal_ratio[[test]],
                  D_abnormal.normal=D_abnormal.normal_ratio[[test]],
                  D0_abnormal.normal=D0_abnormal.normal_ratio[[test]],
                  E_abnormal.normal=E_abnormal.normal_ratio[[test]],
                  F_abnormal.normal=F_abnormal.normal_ratio[[test]],
                  
                  normal_B.A=normal_B.A_ratio[[test]],
                  normal_C.A=normal_C.A_ratio[[test]],
                  normal_D.A=normal_D.A_ratio[[test]],
                  normal_D0.A=normal_D0.A_ratio[[test]],
                  normal_E.A=normal_E.A_ratio[[test]],
                  normal_F.A=normal_F.A_ratio[[test]],
                  
                  abnormal_B.A=abnormal_B.A_ratio[[test]],
                  abnormal_C.A=abnormal_C.A_ratio[[test]],
                  abnormal_D.A=abnormal_D.A_ratio[[test]],
                  abnormal_D0.A=abnormal_D0.A_ratio[[test]],
                  abnormal_E.A=abnormal_E.A_ratio[[test]],
                  abnormal_F.A=abnormal_F.A_ratio[[test]]) 
    return (re.all)
  }
  
  w.re_ratio <-  re_all_ratio(test="w.re")
  t.re_ratio <-  re_all_ratio(test="t.re")
  saveRDS( w.re_ratio, "C:/zhijuncao/cardiotoxicity/DOX_Olink/somascan_w_test_results_ratio.Rds")
  saveRDS(t.re_ratio, "C:/zhijuncao/cardiotoxicity/DOX_Olink/somascan_t_test_results_ratio.Rds")
  
  w.re_df_ratio <- w.re_ratio %>% data.frame() %>% rownames_to_column("SomaId") %>% inner_join(somaid, .)
  w.re_df_sig_ratio <- w_sig(data=w.re_df_ratio, fc=1.2)
  write.xlsx(list(all=w.re_df_ratio,sig=w.re_df_sig_ratio), "C:/zhijuncao/cardiotoxicity/DOX_Olink/somacan_wilcox_result_all_ratio_11102021.xlsx")
  
  w.re_id_ratio <-  map(w.re_ratio, ~.x%>% rownames_to_column("SomaId") %>% inner_join(somaid,.))
  write.xlsx(w.re_id_ratio, "C:/zhijuncao/cardiotoxicity/DOX_Olink/somacan_wilcox_result_separated_ratio_11102021.xlsx") 
  
  t.re_df_ratio <- t.re_ratio %>% data.frame() %>% rownames_to_column("SomaId") %>% inner_join(somaid, .)
  t.re_df_sig_ratio <- t_sig(data=t.re_df_ratio, fc=1.2)
  write.xlsx(list(all=t.re_df_ratio, sig=t.re_df_sig_ratio), "C:/zhijuncao/cardiotoxicity/DOX_Olink/somacan_ttest_result_all_ratio_11102021.xlsx") 
  
  t.re_id_ratio <-  map(t.re_ratio, ~.x%>% rownames_to_column("SomaId") %>% inner_join(somaid,.))
  write.xlsx(t.re_id_ratio, "C:/zhijuncao/cardiotoxicity/DOX_Olink/somacan_ttest_result_separated_ratio_11102021.xlsx") 
  

###volcanoplot
  
pfc <- c("Assay", "w.pvalue", "fc.median","w.FDR")
pdf("C:/zhijuncao/cardiotoxicity/DOX_Olink/somacan_wtest_result_volcanoplot.pdf", width = 6, height = 6)
  map2(w.re_id, names(w.re_id),~volcanoplotfdr(.x[,pfc], title=.y,fccutoff=1.2, label="labelsig",labelsize=3,pointsize=1, pcutoff=0.05, fdrcutoff=2))
dev.off()

pfc <- c("Assay", "t.pvalue", "fc.mean","t.FDR")
pdf("C:/zhijuncao/cardiotoxicity/DOX_Olink/somacan_ttest_result_volcanoplot.pdf", width = 6, height = 6)
map2(t.re_id, names(t.re_id),~volcanoplotfdr(.x[,pfc], title=.y,fccutoff=1.2, label="labelsig",labelsize=3,pointsize=1, pcutoff=0.05, fdrcutoff=2))

dev.off()

##Ratio
pfc <- c("Assay", "w.pvalue", "fc.median","w.FDR")
pdf("C:/zhijuncao/cardiotoxicity/DOX_Olink/somacan_wtest_result_ratio_volcanoplot.pdf", width = 6, height = 6)
map2(w.re_id_ratio, names(w.re_id_ratio),~volcanoplotfdr(.x[,pfc], title=.y,fccutoff=1.2, label="labelsig",labelsize=3,pointsize=1, pcutoff=0.05, fdrcutoff=2))
dev.off()

pfc <- c("Assay", "t.pvalue", "fc.mean","t.FDR")
pdf("C:/zhijuncao/cardiotoxicity/DOX_Olink/somacan_ttest_result_ratio_volcanoplot.pdf", width = 6, height = 6)
map2(t.re_id_ratio, names(t.re_id_ratio),~volcanoplotfdr(.x[,pfc], title=.y,fccutoff=1.2, label="labelsig",labelsize=3,pointsize=1, pcutoff=0.05, fdrcutoff=2))

dev.off()


####point_errobar plot
point_errorbar_plot <-
  function(data = intest05min,
           analytes = names(intest05min)[c(6:7)],
           ids = intest_clean$ids,
           intratio="Intenstity") {
    j <- 0
    for (analyte in analytes) {
      cat(j, "")
      j <- j+1
      title <- paste(analyte, somaid$Assay[somaid$SomaId == analyte], sep=": ")
      p <- ggplot(data, aes(TimePoint, get(analyte))) +
        facet_wrap( ~Status) +
        theme_classic() +
        geom_point(alpha = 0.5) +
        geom_line(aes(group=SubjectID, color=SubjectID))+
        stat_summary(
          fun.data = "mean_sdl",
          fun.args = list(mult = 1),
          geom = "errorbar",
          width = 0.4,
          size=0.5,
          position = position_dodge(0.6)
        ) +
        stat_summary(
          aes(y = get(analyte)),
          fun = "mean",
          geom = "point",
          shape = 4,
          position = position_dodge(0.6)
        ) +
        labs(y = intratio, x = "TimePoint", title = title) +
        theme(
          axis.title = element_text(size = 10),
          axis.text = element_text(size = 10),
          axis.text.x = element_text(angle = 0)
        ) +
        theme(plot.title = element_text(size = 7))+
        theme(legend.key.size = unit(0.3, 'cm'))
      
      print(p)
    }
    
  }
####

#### time point A, B, C
ABC <- sample %>% filter(TimePoint %in% c("A", "B", "C"))
ABC_ratio <- sample_ratio %>% filter(TimePoint %in% c("A", "B", "C"))

pdf("C:/zhijuncao/cardiotoxicity/DOX_Olink/somacan_errorbarplot_ABC.pdf", width = 7, height = 4)
point_errorbar_plot(data = ABC,
                    analytes = names(ABC)[-c(1:7)],
                    ids = somaid,
                    intratio = "Intensity")
dev.off()

pdf("C:/zhijuncao/cardiotoxicity/DOX_Olink/somacan_errorbarplot_ABC_ratio.pdf", width = 7, height = 4)
point_errorbar_plot(data = ABC_ratio,
                    analytes = names(ABC_ratio)[-c(1:7)],
                    ids = somaid,
                    intratio = "Ratio")
dev.off()

#### All time point
pdf("C:/zhijuncao/cardiotoxicity/DOX_Olink/somacan_errorbarplot.pdf", width = 7, height = 4)
point_errorbar_plot(data = sample,
                    analytes = names(sample)[-c(1:7)],
                    ids = somaid)
dev.off()

pdf("C:/zhijuncao/cardiotoxicity/DOX_Olink/somacan_wilcox_sig_errorbarplot.pdf", width = 7, height = 4)
point_errorbar_plot(data = sample,
                    analytes = w.re_df_sig$SomaId,
                    ids = somaid)
dev.off()

sample_ratio


pdf("C:/zhijuncao/cardiotoxicity/DOX_Olink/somacan_wilcox_sig_errorbarplot.pdf", width = 7, height = 4)
point_errorbar_plot(data = sample,
                    analytes = w.re_df_sig$SomaId,
                    ids = somaid)
dev.off()

####Ratio
pdf("C:/zhijuncao/cardiotoxicity/DOX_Olink/somacan_wilcox_sig_ratio_errorbarplot.pdf", width = 7, height = 4)
point_errorbar_plot(data = sample_ratio,
                    analytes = w.re_df_sig_ratio$SomaId,
                    ids = somaid)
dev.off()

sample_ratio


pdf("C:/zhijuncao/cardiotoxicity/DOX_Olink/somacan_ttest_sig_ratio_errorbarplot.pdf", width = 7, height = 4)
point_errorbar_plot(data = sample_ratio,
                    analytes = t.re_df_sig_ratio$SomaId,
                    ids = somaid)
dev.off()











